# gextract(file = "*.parquet"): checked against pyarrow, the reader the files are written for.

# This file re-roots into databases of its own. Put a usable root back for whichever file
# the parallel worker picks up next, which would otherwise inherit this one.
restore_groot_on_exit()

# Path of a python3 that can import pyarrow, or NULL
pyarrow_python <- function() {
    py <- Sys.which("python3")
    if (!nzchar(py)) {
        return(NULL)
    }
    ok <- suppressWarnings(system2(py, c("-I", "-c", shQuote("import pyarrow")), stdout = FALSE, stderr = FALSE))
    if (identical(ok, 0L)) py else NULL
}

skip_without_pyarrow <- function() {
    skip_if(is.null(pyarrow_python()), "python3 with pyarrow is not available")
}

# Python writes floats in hex so that R reads back the exact float32 value, and decodes
# binary to text so that a chrom column written without the STRING annotation still loads
# (and fails the type check rather than the script).
pq_dump_script <- '
import json, sys
import pyarrow.parquet as pq

path, out_tsv, out_json = sys.argv[1:4]

def py(v):
    if isinstance(v, bytes):
        return v.decode()
    if isinstance(v, float):
        return "NaN" if v != v else v.hex()
    return v

f = pq.ParquetFile(path)
md = f.metadata
meta = {"schema": [{"name": x.name, "type": str(x.type)} for x in f.schema_arrow], "row_groups": []}
for i in range(md.num_row_groups):
    rg = md.row_group(i)
    cols = []
    for j in range(rg.num_columns):
        c = rg.column(j)
        col = {"name": c.path_in_schema, "encodings": list(c.encodings)}
        s = c.statistics
        if s is not None and s.has_min_max:
            col["min"] = py(s.min)
            col["max"] = py(s.max)
        if s is not None and s.has_null_count:
            col["null_count"] = s.null_count
        cols.append(col)
    meta["row_groups"].append({"num_rows": rg.num_rows, "columns": cols})
with open(out_json, "w") as fh:
    json.dump(meta, fh)
    fh.write("\\n")

t = f.read()
cols = [t.column(k).to_pylist() for k in range(t.num_columns)]
with open(out_tsv, "w") as fh:
    fh.write("\\t".join(t.column_names) + "\\n")
    for row in zip(*cols):
        fh.write("\\t".join(str(py(v)) for v in row) + "\\n")
'

# Reads a Parquet file with pyarrow: the table (value columns as doubles holding the exact
# float32 values), the arrow schema (name, type) and per-row-group metadata.
read_pq <- function(path) {
    script <- withr::local_tempfile(fileext = ".py")
    out_tsv <- withr::local_tempfile(fileext = ".tsv")
    out_json <- withr::local_tempfile(fileext = ".json")
    writeLines(pq_dump_script, script)
    status <- system2(pyarrow_python(), c("-I", shQuote(script), shQuote(path), shQuote(out_tsv), shQuote(out_json)))
    if (!identical(status, 0L)) {
        stop("pyarrow could not read ", path)
    }
    tab <- utils::read.delim(out_tsv,
        colClasses = "character", check.names = FALSE, quote = "",
        comment.char = "", na.strings = character()
    )
    for (k in names(tab)[-1]) {
        tab[[k]] <- as.numeric(tab[[k]])
    }
    meta <- yaml::read_yaml(out_json, handlers = list(int = as.numeric))
    list(
        table = tab,
        names = vapply(meta$schema, function(x) x$name, ""),
        types = vapply(meta$schema, function(x) x$type, ""),
        row_groups = meta$row_groups
    )
}

rg_column <- function(rg, name) {
    Filter(function(x) identical(x$name, name), rg$columns)[[1]]
}

# A fresh database rooted for the calling test only
local_pq_db <- function(sizes, env = parent.frame()) {
    local_db_state(env)
    root <- file.path(withr::local_tempdir(.local_envir = env), "db")
    create_test_db(root, chrom_sizes = data.frame(chrom = names(sizes), size = unname(sizes)))
    suppressMessages(gdb.init(root))
    invisible(root)
}

bins_of <- function(chrom, size, binsize) {
    start <- seq(0, size - 1, by = binsize)
    data.frame(chrom = chrom, start = start, end = pmin(start + binsize, size))
}

# chr1 = 1000 bp, chr2 = 600 bp; dense track pq (binsize 32) = start / 32 + 0.25, + 1000 on chr2
local_pq_track_db <- function(env = parent.frame()) {
    local_pq_db(c(chr1 = 1000, chr2 = 600), env)
    intervs <- rbind(bins_of("chr1", 1000, 32), bins_of("chr2", 600, 32))
    gtrack.create_dense("pq", "pq", intervs, intervs$start / 32 + 0.25 + ifelse(intervs$chrom == "chr2", 1000, 0), 32)
}

# Rounds to the nearest float32 (as.single() only tags a vector for .C and rounds nothing)
float32 <- function(x) {
    readBin(writeBin(as.numeric(x), raw(), size = 4), "numeric", size = 4, n = length(x))
}

expect_same_rows <- function(got, inmem, value_cols) {
    expect_equal(got$chrom, as.character(inmem$chrom))
    expect_equal(got$start, as.numeric(inmem$start))
    expect_equal(got$end, as.numeric(inmem$end))
    for (k in value_cols) {
        expect_identical(got[[k]], float32(inmem[[k]]), label = k)
    }
}

test_that("parquet output equals in-memory gextract", {
    skip_without_pyarrow()
    local_pq_track_db()
    f <- file.path(withr::local_tempdir(), "x.parquet")
    exprs <- c("pq", "pq * 2", "pq + 0.5")
    cols <- c("a", "a.b c", "test_pq * 2")

    inmem <- gextract(exprs, intervals = gintervals.all(), iterator = 32, colnames = cols)
    gextract(exprs, intervals = gintervals.all(), iterator = 32, colnames = cols, file = f)
    pq <- read_pq(f)

    expect_identical(pq$names, c("chrom", "start", "end", cols))
    expect_identical(pq$types, c("string", "int32", "int32", "float", "float", "float"))
    expect_equal(nrow(pq$table), 32 + 19)
    expect_same_rows(pq$table, inmem, cols)

    expect_length(pq$row_groups, 2)
    for (i in 1:2) {
        chrom <- c("chr1", "chr2")[i]
        st <- rg_column(pq$row_groups[[i]], "chrom")
        expect_identical(st$min, chrom)
        expect_identical(st$max, chrom)
        expect_equal(st$null_count, 0)
        start_st <- rg_column(pq$row_groups[[i]], "start")
        expect_equal(start_st$min, min(inmem$start[inmem$chrom == chrom]))
        expect_equal(start_st$max, max(inmem$start[inmem$chrom == chrom]))
        expect_null(rg_column(pq$row_groups[[i]], "a")$min)
    }

    # the extension is matched case-insensitively
    f2 <- file.path(dirname(f), "X.PARQUET")
    gextract("pq", intervals = gintervals.all(), iterator = 32, file = f2)
    expect_identical(readBin(f2, "raw", 4), charToRaw("PAR1"))
})

test_that("many columns (Thrift lists of 15 elements or more), values rounded to float32", {
    skip_without_pyarrow()
    local_pq_track_db()
    f <- file.path(withr::local_tempdir(), "x.parquet")
    # pq / 3 and the like are not float32 values, so this checks the rounding too
    exprs <- paste0("pq / ", 1:20)

    inmem <- gextract(exprs, intervals = gintervals.all(), iterator = 32)
    gextract(exprs, intervals = gintervals.all(), iterator = 32, file = f)
    pq <- read_pq(f)

    expect_identical(pq$names, c("chrom", "start", "end", exprs))
    expect_false(identical(float32(inmem[["pq / 3"]]), inmem[["pq / 3"]]))
    expect_same_rows(pq$table, inmem, exprs)
})

test_that("row groups split at the cap and at chromosome changes", {
    skip_without_pyarrow()
    local_pq_db(c(chr1 = 1100000, chr2 = 1000))
    intervs <- rbind(bins_of("chr1", 1100000, 1), bins_of("chr2", 1000, 1))
    gtrack.create_dense("pos", "pos", intervs, intervs$start, 1)
    f <- file.path(withr::local_tempdir(), "x.parquet")
    exprs <- c("pos", "pos * 0")

    inmem <- gextract(exprs, intervals = gintervals.all(), iterator = 1)
    gextract(exprs, intervals = gintervals.all(), iterator = 1, file = f)
    pq <- read_pq(f)

    expect_equal(vapply(pq$row_groups, function(rg) rg$num_rows, 0), c(1048576, 51424, 1000))
    expect_false("RLE_DICTIONARY" %in% rg_column(pq$row_groups[[1]], "pos")$encodings)
    for (rg in pq$row_groups) {
        expect_true("RLE_DICTIONARY" %in% rg_column(rg, "pos * 0")$encodings)
    }
    expect_same_rows(pq$table, inmem, exprs)
})

test_that("nan column stays dictionary-encoded", {
    skip_without_pyarrow()
    local_pq_db(c(chr1 = 1100000, chr2 = 1000))
    intervs <- rbind(bins_of("chr1", 1100000, 1000), bins_of("chr2", 1000, 1000))
    intervs$end <- intervs$start + 1
    gtrack.create_dense("sparse_vals", "sparse", intervs, intervs$start, 1, defval = NaN)
    f <- file.path(withr::local_tempdir(), "x.parquet")

    inmem <- gextract("sparse_vals", intervals = gintervals.all(), iterator = 1)
    gextract("sparse_vals", intervals = gintervals.all(), iterator = 1, file = f)
    pq <- read_pq(f)

    for (rg in pq$row_groups) {
        expect_true("RLE_DICTIONARY" %in% rg_column(rg, "sparse_vals")$encodings)
    }
    expect_identical(is.na(pq$table$sparse_vals), is.na(inmem$sparse_vals))
    expect_gt(sum(!is.na(inmem$sparse_vals)), 1000)
    expect_same_rows(pq$table, inmem, "sparse_vals")
})

test_that("a chromosome of 2^31 bp or more gets INT64 coordinates", {
    skip_without_pyarrow()
    local_pq_db(c(chr1 = 10))
    # the size misha reads; no sequence of that length is needed for this
    writeLines("chr1\t2147483648", file.path(.misha$GROOT, "chrom_sizes.txt"))
    suppressMessages(gdb.init(.misha$GROOT))
    gtrack.create_dense("big", "big", gintervals("chr1", 2147480000, 2147483648), 5, 1e6)
    f <- file.path(withr::local_tempdir(), "x.parquet")
    scope <- gintervals("chr1", 2140000000, 2147483648)

    inmem <- gextract("big", intervals = scope, iterator = 1e6)
    gextract("big", intervals = scope, iterator = 1e6, file = f)
    pq <- read_pq(f)

    expect_identical(pq$types, c("string", "int64", "int64", "float"))
    expect_equal(max(pq$table$end), 2147483648)
    expect_equal(rg_column(pq$row_groups[[1]], "end")$max, 2147483648)
    expect_same_rows(pq$table, inmem, "big")
})

test_that("errors leave no file behind", {
    local_pq_track_db()
    dir <- withr::local_tempdir()
    f <- file.path(dir, "x.parquet")

    # (a) duplicate column names, (b) a value column named like a coordinate column
    expect_error(gextract("pq", "pq * 2", intervals = gintervals.all(), iterator = 32, colnames = c("v", "v"), file = f), "v")
    expect_false(file.exists(f))
    expect_error(gextract("pq", intervals = gintervals.all(), iterator = 32, colnames = "start", file = f), "start")
    expect_false(file.exists(f))

    # (c) an error in the middle of the scan: chr2's pq values are 1000.25, 1001.25, ..., so
    # with 8-row evaluation buffers the error comes after the chr1 row group was written
    assign("parquet_test_fail", function(x) {
        if (any(x >= 1010, na.rm = TRUE)) stop("parquet_test_fail saw chr2")
        x
    }, envir = globalenv())
    withr::defer(rm("parquet_test_fail", envir = globalenv()))
    withr::local_options(gbuf.size = 8)
    # capture.output keeps R's own report of the failed expression off the console
    fail_scan <- function() {
        utils::capture.output(
            gextract("parquet_test_fail(pq)", intervals = gintervals.all(), iterator = 32, file = f),
            type = "message"
        )
    }
    expect_error(fail_scan(), "saw chr2")
    expect_false(file.exists(f))
    expect_identical(list.files(dir, all.files = TRUE, no.. = TRUE), character())

    # (d) the same error does not touch an existing target
    writeLines("existing content", f)
    before <- readBin(f, "raw", file.size(f))
    expect_error(fail_scan(), "saw chr2")
    expect_identical(readBin(f, "raw", file.size(f)), before)
    expect_identical(list.files(dir, all.files = TRUE, no.. = TRUE), "x.parquet")
})

test_that("tsv output unchanged", {
    local_pq_track_db()
    f <- file.path(withr::local_tempdir(), "x.tsv")
    gextract("pq", "pq * 2", intervals = gintervals.all(), iterator = 32, colnames = c("a", "b"), file = f)
    lines <- readLines(f)
    expect_identical(lines[1:2], c("chrom\tstart\tend\ta\tb", "chr1\t0\t32\t0.25\t0.5"))
    expect_length(lines, 1 + 32 + 19)
})
