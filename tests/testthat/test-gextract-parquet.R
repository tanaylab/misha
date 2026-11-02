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
import pyarrow as pa
import pyarrow.parquet as pq

path, out_tsv, out_json, coord_type = sys.argv[1:5]

def py(v):
    if isinstance(v, bytes):
        return v.decode()
    if isinstance(v, float):
        return "NaN" if v != v else v.hex()
    return v

f = pq.ParquetFile(path)
md = f.metadata
t = f.read()

# The same table written by pyarrow itself: chrom string, start and end int32 (int64 if
# coord_type says so), the values float32, all nullable
coord = pa.int64() if coord_type == "int64" else pa.int32()
ref_schema = pa.schema([("chrom", pa.string()), ("start", coord), ("end", coord)] + [(n, pa.float32()) for n in t.column_names[3:]])
ref_path = out_json + ".ref.parquet"
pq.write_table(pa.Table.from_pydict({n: t.column(n).to_pylist() for n in t.column_names}, schema=ref_schema), ref_path)
ref = pq.read_table(ref_path)
try:
    pa.concat_tables([t, ref])
    concat_ok = True
except Exception:
    concat_ok = False

meta = {"schema": [{"name": x.name, "type": str(x.type), "nullable": x.nullable} for x in f.schema_arrow], "row_groups": [],
        "schema_equals_reference": t.schema.equals(ref.schema), "concat_ok": concat_ok}
for i in range(md.num_row_groups):
    rg = md.row_group(i)
    cols = []
    for j in range(rg.num_columns):
        c = rg.column(j)
        col = {"name": c.path_in_schema, "encodings": list(c.encodings), "compression": c.compression}
        s = c.statistics
        if s is not None and s.has_min_max:
            col["min"] = py(s.min)
            col["max"] = py(s.max)
        if s is not None and s.has_null_count:
            col["null_count"] = s.null_count
        cols.append(col)
    meta["row_groups"].append({"num_rows": rg.num_rows, "columns": cols})
# also iter_batches() and read_table(): a file with no row groups reads with read_table but not iter_batches
meta["iter_batches_rows"] = sum(b.num_rows for b in f.iter_batches())
meta["read_table_rows"] = pq.read_table(path).num_rows
with open(out_json, "w") as fh:
    json.dump(meta, fh)
    fh.write("\\n")

cols = [t.column(k).to_pylist() for k in range(t.num_columns)]
with open(out_tsv, "w") as fh:
    fh.write("\\t".join(t.column_names) + "\\n")
    for row in zip(*cols):
        fh.write("\\t".join(str(py(v)) for v in row) + "\\n")
'

# Reads a Parquet file with pyarrow: the table (value columns as doubles holding the exact
# float32 values), the arrow schema (name, type, nullable), per-row-group metadata, and how
# the file compares with pyarrow's own write of the same table (coord_type: its start and end).
read_pq <- function(path, coord_type = "int32") {
    script <- withr::local_tempfile(fileext = ".py")
    out_tsv <- withr::local_tempfile(fileext = ".tsv")
    out_json <- withr::local_tempfile(fileext = ".json")
    withr::defer(unlink(paste0(out_json, ".ref.parquet")))
    writeLines(pq_dump_script, script)
    status <- system2(pyarrow_python(), c("-I", shQuote(script), shQuote(path), shQuote(out_tsv), shQuote(out_json), coord_type))
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
        nullable = vapply(meta$schema, function(x) x$nullable, TRUE),
        row_groups = meta$row_groups,
        iter_batches_rows = meta$iter_batches_rows,
        read_table_rows = meta$read_table_rows,
        schema_equals_reference = meta$schema_equals_reference,
        concat_ok = meta$concat_ok
    )
}

rg_column <- function(rg, name) {
    Filter(function(x) identical(x$name, name), rg$columns)[[1]]
}

# Every column nullable and every column chunk GZIP, as arrow writes them: the arrow schema
# equals that of pyarrow's own write of the table, and the two tables concatenate
expect_arrow_like <- function(pq) {
    expect_true(all(pq$nullable))
    expect_true(pq$schema_equals_reference)
    expect_true(pq$concat_ok)
    codecs <- unlist(lapply(pq$row_groups, function(rg) vapply(rg$columns, function(x) x$compression, "")))
    expect_true(length(codecs) > 0 && all(codecs == "GZIP"), label = paste(unique(codecs), collapse = ", "))
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
    expect_arrow_like(pq)
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

    # created_by is "<application> version <version>", the form parquet-mr parses
    bytes <- readBin(f, "raw", file.size(f))
    footer <- rawToChar(bytes[bytes != as.raw(0)])
    expect_match(footer, paste0("misha version ", as.character(packageVersion("misha"))), fixed = TRUE, useBytes = TRUE)
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
    expect_arrow_like(pq)
    expect_false(identical(float32(inmem[["pq / 3"]]), inmem[["pq / 3"]]))
    expect_same_rows(pq$table, inmem, exprs)
})

test_that("an empty result is one empty row group, as pyarrow writes an empty table", {
    skip_without_pyarrow()
    local_pq_track_db()
    f <- file.path(withr::local_tempdir(), "x.parquet")
    # an iterator that does not meet the scope yields no rows
    scope <- gintervals("chr1", 0, 100)
    iter <- gintervals("chr2", 0, 100)

    expect_null(gextract("pq", "pq * 2", intervals = scope, iterator = iter))
    gextract("pq", "pq * 2", intervals = scope, iterator = iter, file = f)
    pq <- read_pq(f)

    expect_identical(pq$names, c("chrom", "start", "end", "pq", "pq * 2"))
    expect_identical(pq$types, c("string", "int32", "int32", "float", "float"))
    expect_arrow_like(pq)
    expect_equal(nrow(pq$table), 0)
    expect_equal(pq$iter_batches_rows, 0)
    expect_equal(pq$read_table_rows, 0)
    expect_equal(vapply(pq$row_groups, function(rg) rg$num_rows, 0), 0)
})

test_that("row groups split at the cap and at chromosome changes", {
    skip_without_pyarrow()
    local_pq_db(c(chr1 = 1100000, chr2 = 1000))
    # one writer: kids would cut chr1 into shards of their own
    withr::local_options(gmultitasking = FALSE)
    intervs <- rbind(bins_of("chr1", 1100000, 1), bins_of("chr2", 1000, 1))
    gtrack.create_dense("pos", "pos", intervs, intervs$start, 1)
    f <- file.path(withr::local_tempdir(), "x.parquet")
    exprs <- c("pos", "pos * 0")

    inmem <- gextract(exprs, intervals = gintervals.all(), iterator = 1)
    gextract(exprs, intervals = gintervals.all(), iterator = 1, file = f)
    pq <- read_pq(f)

    expect_equal(vapply(pq$row_groups, function(rg) rg$num_rows, 0), c(1048576, 51424, 1000))
    # a full row group: one data page per column chunk, its definition levels one RLE run
    expect_arrow_like(pq)
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
    pq <- read_pq(f, "int64")

    expect_identical(pq$types, c("string", "int64", "int64", "float"))
    expect_arrow_like(pq)
    expect_equal(max(pq$table$end), 2147483648)
    expect_equal(rg_column(pq$row_groups[[1]], "end")$max, 2147483648)
    expect_same_rows(pq$table, inmem, "big")
})

test_that("errors leave no file behind", {
    local_pq_track_db()
    dir <- withr::local_tempdir()
    f <- file.path(dir, "x.parquet")

    # (a) duplicate column names, (b) a value column named like a coordinate column
    expect_error(gextract("pq", "pq * 2", intervals = gintervals.all(), iterator = 32, colnames = c("v", "v"), file = f), "two columns named \"v\"")
    expect_false(file.exists(f))
    expect_error(gextract("pq", intervals = gintervals.all(), iterator = 32, colnames = "start", file = f), "value column named \"start\"")
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

# Multitasking: kids write parts that the parent merges. With gmax.processes = 4 and
# gmin.scope4process = 1000, prepare4multitasking() plans these shards:
# - an implicit iterator or a sparse one: whole chromosomes, chr1 | chr2 | chr3 chr4 | chr5 chr6
#   for six chromosomes of equal size, chr1 | chr2 | chr3 | chr4 chr5 chr6 when chr3 is twice
#   as long as the others
# - iterator = 50 over 6 x 1,300,000 bp: four ranges that cut chromosomes

# chr1..chr6 of the given sizes; dense track six (binsize 50) = start / 50 + chrom index / 8,
# so the fraction tells the chromosome. The values are small enough for gtrack.import_parquet
# to read them back exactly: it averages with float products, as gtrackimportwig does.
local_six_chrom_db <- function(sizes = rep(200000, 6), env = parent.frame()) {
    chroms <- paste0("chr", 1:6)
    local_pq_db(stats::setNames(sizes, chroms), env)
    intervs <- do.call(rbind, Map(bins_of, chroms, sizes, 50))
    gtrack.create_dense("six", "six", intervs, intervs$start / 50 + match(intervs$chrom, chroms) / 8, 50)
}

dir_entries <- function(dir) list.files(dir, all.files = TRUE, no.. = TRUE)

test_that("multitask output equals serial output", {
    skip_without_pyarrow()
    local_six_chrom_db(rep(1300000, 6))
    withr::local_options(gmax.processes = 4, gmin.scope4process = 1000)
    dir <- withr::local_tempdir()
    f <- file.path(dir, "x.parquet")
    exprs <- c("six", "six * 2", "six + 0.5")
    cols <- c("a", "b", "c")

    # NULL: shards of whole chromosomes; 50: ranges that cut chromosomes, so a chromosome's
    # rows come from two kids
    for (iter in list(NULL, 50)) {
        label <- if (is.null(iter)) "implicit iterator" else "iterator = 50"
        inmem <- gextract(exprs, intervals = gintervals.all(), iterator = iter, colnames = cols)
        gextract(exprs, intervals = gintervals.all(), iterator = iter, colnames = cols, file = f)
        pq <- read_pq(f)

        expect_equal(nrow(pq$table), 6 * 26000, label = label)
        expect_arrow_like(pq)
        expect_same_rows(pq$table, inmem, cols)
        # one chromosome per row group; the rows' order is checked above
        rg_chroms <- vapply(pq$row_groups, function(rg) rg_column(rg, "chrom")$min, "")
        expect_identical(rg_chroms, vapply(pq$row_groups, function(rg) rg_column(rg, "chrom")$max, ""), label = label)
        expect_identical(unique(rg_chroms), paste0("chr", 1:6), label = label)
        if (is.null(iter)) {
            expect_length(rg_chroms, 6)
        } else {
            expect_gt(length(rg_chroms), 6)
        }
        expect_identical(dir_entries(dir), "x.parquet", label = label)
    }

    # misha reads it back
    gtrack.import_parquet("back", f, 50, track.prefix = "back_")
    for (k in cols) {
        expect_identical(gextract(paste0("back_", k), gintervals.all(), iterator = 50)[[paste0("back_", k)]], float32(inmem[[k]]), label = k)
    }
})

test_that("a kid with no rows still merges", {
    skip_without_pyarrow()
    local_six_chrom_db(c(2, 2, 4, 2, 2, 2) * 100000)
    withr::local_options(gmax.processes = 4, gmin.scope4process = 1000)
    # no intervals on chr3, which is a shard of its own
    iter <- gintervals(paste0("chr", c(1, 2, 4, 5, 6)), 1000, 1100)
    gtrack.create_sparse("sp", "sp", iter, seq_len(nrow(iter)))
    f <- file.path(withr::local_tempdir(), "x.parquet")

    inmem <- gextract("six", "sp", intervals = gintervals.all(), iterator = "sp")
    gextract("six", "sp", intervals = gintervals.all(), iterator = "sp", file = f)
    pq <- read_pq(f)

    expect_same_rows(pq$table, inmem, c("six", "sp"))
    expect_arrow_like(pq)
    expect_identical(vapply(pq$row_groups, function(rg) rg$num_rows, 0), rep(1, 5))
})

test_that("a failing kid leaves nothing behind", {
    local_six_chrom_db()
    withr::local_options(gmax.processes = 4, gmin.scope4process = 1000)
    dir <- withr::local_tempdir()
    f <- file.path(dir, "x.parquet")
    # fails on chr5 (fraction 5 / 8), in the last of the four kids
    assign("parquet_test_fail", function(x) {
        if (any(x %% 1 == 0.625, na.rm = TRUE)) stop("parquet_test_fail saw chr5")
        x
    }, envir = globalenv())
    withr::defer(rm("parquet_test_fail", envir = globalenv()))
    fail_scan <- function() {
        utils::capture.output(
            gextract("parquet_test_fail(six)", intervals = gintervals.all(), file = f),
            type = "message"
        )
    }

    expect_error(fail_scan(), "saw chr5")
    expect_identical(.glifetime_counters(), c(ref_count = 0L, protect_count = 0L))
    expect_identical(dir_entries(dir), character())

    # an existing target is left as it was
    writeLines("existing content", f)
    expect_error(fail_scan(), "saw chr5")
    expect_identical(readLines(f), "existing content")
    expect_identical(dir_entries(dir), "x.parquet")
})

test_that("a target in a missing directory is an error naming it, not its temporary file", {
    local_six_chrom_db()
    withr::local_options(gmax.processes = 4, gmin.scope4process = 1000)
    f <- file.path(withr::local_tempdir(), "nodir", "x.parquet")
    msg <- paste0("Failed to open file ", f, " for writing: No such file or directory")
    # one process, then kids (each fails to create its part)
    for (multitask in c(FALSE, TRUE)) {
        withr::with_options(list(gmultitasking = multitask), {
            expect_error(gextract("six", intervals = gintervals.all(), file = f), msg, fixed = TRUE)
        })
    }
})

test_that("every part takes INT64 coordinates when one chromosome of the scope needs them", {
    skip_without_pyarrow()
    local_pq_db(c(chr1 = 10, chr2 = 10000))
    # the sizes misha reads; no sequence of that length is needed for this
    writeLines(c("chr1\t2147483648", "chr2\t10000"), file.path(.misha$GROOT, "chrom_sizes.txt"))
    suppressMessages(gdb.init(.misha$GROOT))
    gtrack.create_dense("big", "big", gintervals(c("chr1", "chr2"), c(2147480000, 0), c(2147483648, 10000)), c(5, 7), 1e6)
    withr::local_options(gmax.processes = 4, gmin.scope4process = 1000)
    f <- file.path(withr::local_tempdir(), "x.parquet")
    # one shard per chromosome: the chr2 kid alone would pick INT32
    scope <- gintervals(c("chr1", "chr2"), c(2147000000, 0), c(2147483648, 10000))

    inmem <- gextract("big", intervals = scope)
    gextract("big", intervals = scope, file = f)
    pq <- read_pq(f, "int64")

    expect_identical(pq$types, c("string", "int64", "int64", "float"))
    expect_arrow_like(pq)
    expect_same_rows(pq$table, inmem, "big")
})

test_that("a symlinked target is written through and the link kept", {
    skip_without_pyarrow()
    local_six_chrom_db()
    withr::local_options(gmax.processes = 4, gmin.scope4process = 1000)
    inmem <- gextract("six", intervals = gintervals.all())

    for (multitask in c(FALSE, TRUE)) {
        label <- if (multitask) "multitask" else "serial"
        parent <- withr::local_tempdir()
        link_dir <- file.path(parent, "links")
        real_dir <- file.path(parent, "real")
        dir.create(link_dir)
        dir.create(real_dir)
        for (f in file.path(real_dir, c("abs.parquet", "rel.parquet"))) {
            writeLines("old content", f)
        }
        # links/<name>.parquet -> the text of each
        links <- c(
            abs = file.path(real_dir, "abs.parquet"), # to an existing file
            rel = file.path("..", "real", "rel.parquet"), # relative, to an existing file
            dangling = "chain.parquet", # relative, to a link to a file that does not exist yet
            chain = file.path(real_dir, "new.parquet")
        )
        link_paths <- file.path(link_dir, paste0(names(links), ".parquet"))
        file.symlink(unname(links), link_paths)

        withr::with_options(list(gmultitasking = multitask), {
            for (k in c("abs", "rel", "dangling")) {
                gextract("six", intervals = gintervals.all(), file = file.path(link_dir, paste0(k, ".parquet")))
            }
        })

        expect_identical(Sys.readlink(link_paths), unname(links), label = label)
        expect_setequal(dir_entries(link_dir), basename(link_paths))
        expect_setequal(dir_entries(real_dir), c("abs.parquet", "rel.parquet", "new.parquet"))
        for (f in dir_entries(real_dir)) {
            expect_same_rows(read_pq(file.path(real_dir, f))$table, inmem, "six")
        }
    }
})

test_that("temporary files are named for this host and process", {
    local_six_chrom_db()
    withr::local_options(gmax.processes = 4, gmin.scope4process = 1000)
    out_dir <- withr::local_tempdir()
    f <- file.path(out_dir, "x.parquet")
    # what the directory holds at each evaluation of the expression, in kids too
    seen_file <- withr::local_tempfile()
    assign("parquet_test_peek", function(x) {
        cat(sprintf("%s\n", dir_entries(out_dir)), file = seen_file, sep = "", append = TRUE)
        x
    }, envir = globalenv())
    withr::defer(rm("parquet_test_peek", envir = globalenv()))
    tmp <- sprintf(".x.parquet.tmp.%s.%d", .gdb.host_tag(), Sys.getpid())
    parts <- paste0(tmp, ".", 0:3)

    withr::with_options(list(gmultitasking = FALSE), gextract("parquet_test_peek(six)", intervals = gintervals.all(), file = f))
    seen <- unique(readLines(seen_file))
    expect_identical(seen, tmp)

    unlink(c(f, seen_file))
    gextract("parquet_test_peek(six)", intervals = gintervals.all(), file = f)
    seen <- unique(readLines(seen_file))
    expect_true(all(seen %in% c(tmp, parts)), label = paste(seen, collapse = ", "))
    expect_true(any(seen %in% parts), label = paste(seen, collapse = ", "))
    expect_identical(dir_entries(out_dir), "x.parquet")
})

test_that("leftovers of killed writers on this host are removed, and nothing else", {
    local_six_chrom_db()
    withr::local_options(gmax.processes = 4, gmin.scope4process = 1000)
    host <- .gdb.host_tag()
    dead <- 4194303L # the largest pid Linux allows; running only if the pids wrapped around to it
    skip_if(dir.exists(file.path("/proc", dead)))
    # a live process of this user
    sleeper <- as.integer(system("sleep 600 >/dev/null 2>&1 & echo $!", intern = TRUE))
    withr::defer(tools::pskill(sleeper))

    for (multitask in c(FALSE, TRUE)) {
        out_dir <- withr::local_tempdir()
        stale <- sprintf(".x.parquet.tmp.%s.%d%s", host, dead, c("", ".0", ".11"))
        kept <- c(
            sprintf(".x.parquet.tmp.someotherhost.%d", dead),
            sprintf(".x.parquet.tmp.%s.%d%s", host, sleeper, c("", ".2")),
            # not exactly .<target>.tmp.<host>.<pid>[.<kid>]
            sprintf(".x.parquet.tmp.%s.0%d", host, dead),
            sprintf(".x.parquet.tmp.%s.%d%s", host, dead, c(".x", ".1.2", ".")),
            sprintf(".x.parquet.tmp.%d", dead),
            sprintf(".y.parquet.tmp.%s.%d", host, dead)
        )
        file.create(file.path(out_dir, c(stale, kept)))

        withr::with_options(list(gmultitasking = multitask), gextract("six", intervals = gintervals.all(), file = file.path(out_dir, "x.parquet")))
        expect_setequal(dir_entries(out_dir), c("x.parquet", kept))
    }
})
