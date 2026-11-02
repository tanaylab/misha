# gtrack.import_parquet: Parquet columns as dense tracks. The fixtures and the script that
# writes them (pyarrow) are in fixtures/parquet; see make_fixtures.py for what each one holds.

# Every test re-roots into a database of its own.
restore_groot_on_exit()

pq_fixture <- function(name) {
    normalizePath(test_path("fixtures", "parquet", name), mustWork = FALSE)
}

# chr1 (1000 bp) and chr2 (600 bp): 32 + 19 bins of 32 bp
local_pq_db <- function(env = parent.frame()) {
    local_db_state(env)
    tmp <- withr::local_tempdir(.local_envir = env)
    db <- file.path(tmp, "pqdb")
    create_test_db(db, data.frame(chrom = c("chr1", "chr2"), size = c(1000, 600)))
    gdb.init(db)
    db
}

bin_k <- c(0:31, 0:18)
bin_chrom_index <- rep(0:1, c(32, 19))

# The fixtures' closed form: a = k + 0.25, b = chrom_index * 1000 + k (k = start / 32)
expected_values <- function(value_col) {
    if (value_col == "a") bin_k + 0.25 else as.numeric(as.single(bin_chrom_index * 1000 + bin_k))
}

track_values <- function(track) {
    gextract(track, gintervals.all(), iterator = 32)[[track]]
}

# Hidden staging dirs left in the tracks dir (.trash.* entries are being removed already)
staging_leftovers <- function(db) {
    entries <- list.files(file.path(db, "tracks"), all.files = TRUE, no.. = TRUE)
    entries[grepl("\\.tmp\\.", entries) & !grepl("^\\.trash\\.", entries)]
}

# Nothing of tracks is left: not listed, not in GTRACKS, no track dir, no staging dir
expect_no_tracks_left <- function(db, tracks) {
    expect_length(gtrack.ls(), 0)
    expect_false(any(tracks %in% get("GTRACKS", envir = .misha)))
    expect_false(any(dir.exists(file.path(db, "tracks", paste0(tracks, ".track")))))
    expect_length(staging_leftovers(db), 0)
}

# Copies fixtures (a list of new name = fixture) into a fresh directory; returns the directory
local_fixture_dir <- function(fixtures, env = parent.frame()) {
    src <- withr::local_tempdir(.local_envir = env)
    for (i in seq_along(fixtures)) {
        file.copy(pq_fixture(fixtures[[i]]), file.path(src, names(fixtures)[i]))
    }
    src
}

test_that("imports every fixture with exact values", {
    db <- local_pq_db()
    fixtures <- c(
        "pandas_default", "plain_uncompressed", "gzip_v2", "mixed_rowgroups",
        "dict_fallback", "byte_stream_split", "unknown_fields", "shard_chr*"
    )
    for (fx in fixtures) {
        prefix <- paste0("t_", gsub("*", "", fx, fixed = TRUE), "_")
        created <- gtrack.import_parquet("test", pq_fixture(paste0(fx, ".parquet")), 32, track.prefix = prefix)
        expect_identical(created, paste0(prefix, c("a", "b")), info = fx)
        for (value_col in c("a", "b")) {
            expect_identical(track_values(paste0(prefix, value_col)), expected_values(value_col), info = paste(fx, value_col))
        }
    }
})

test_that("nulls become NaN", {
    db <- local_pq_db()
    gtrack.import_parquet("test", pq_fixture("nulls.parquet"), 32, track.prefix = "n_")
    expected <- expected_values("a")
    expected[bin_k %% 7 == 0] <- NaN
    expect_identical(track_values("n_a"), expected)
    expect_identical(track_values("n_b"), expected_values("b"))
})

test_that("a column chunk that is entirely null takes defval", {
    db <- local_pq_db()
    for (fx in c("allnull_dict", "allnull_plain")) {
        for (defval in c(NaN, -7)) {
            prefix <- paste0(fx, "_", if (is.nan(defval)) "nan" else "num", "_")
            gtrack.import_parquet("test", pq_fixture(paste0(fx, ".parquet")), 32, track.prefix = prefix, defval = defval)
            expected <- expected_values("a")
            expected[bin_chrom_index == 1] <- defval
            expect_identical(track_values(paste0(prefix, "a")), expected, label = prefix)
            expect_identical(track_values(paste0(prefix, "b")), expected_values("b"), label = prefix)
        }
    }
})

test_that("a column chunk that claims more values than its row group has is refused up front", {
    db <- local_pq_db()
    expect_error(
        gtrack.import_parquet("test", pq_fixture("huge_num_values.parquet"), 32, track.prefix = "h_"),
        "column chrom: .*number of values differs from the number of rows"
    )
    expect_no_tracks_left(db, c("h_a"))
})

test_that("a definition level the bit width cannot hold is corrupt data", {
    db <- local_pq_db()
    expect_error(
        gtrack.import_parquet("test", pq_fixture("bad_def_level.parquet"), 32, track.prefix = "d_"),
        "bad_def_level.parquet, row group 0, column a: corrupt data \\(RLE run value exceeds the bit width\\)"
    )
    expect_no_tracks_left(db, "d_a")
})

test_that("a row group of more rows than the genome has bases is refused before it is read", {
    db <- local_pq_db()
    expect_error(
        gtrack.import_parquet("test", pq_fixture("more_rows_than_bases.parquet"), 32, track.prefix = "m_"),
        "more_rows_than_bases.parquet, row group 0: 1601 rows, more than the 1600 bases of the genome"
    )
    expect_no_tracks_left(db, "m_a")
})

test_that("description is validated up front", {
    db <- local_pq_db()
    expect_error(gtrack.import_parquet(NULL, pq_fixture("single.parquet"), 32), "Usage: gtrack.import_parquet")
    expect_no_tracks_left(db, "a")
})

test_that("columns selects and orders tracks", {
    db <- local_pq_db()
    created <- gtrack.import_parquet("test", pq_fixture("pandas_default.parquet"), 32, track.prefix = "s_", columns = c("b", "a"))
    expect_identical(created, c("s_b", "s_a"))
    expect_identical(track_values("s_b"), expected_values("b"))
    expect_identical(track_values("s_a"), expected_values("a"))

    created <- gtrack.import_parquet("test", pq_fixture("pandas_default.parquet"), 32, track.prefix = "o_", columns = "b")
    expect_identical(created, "o_b")
    expect_false(gtrack.exists("o_a"))
})

test_that("columns = NULL with a string column errors and creates nothing", {
    db <- local_pq_db()
    expect_error(
        gtrack.import_parquet("test", pq_fixture("with_name.parquet"), 32, track.prefix = "w_"),
        "columns"
    )
    expect_length(gtrack.ls(), 0)
    expect_length(staging_leftovers(db), 0)

    created <- gtrack.import_parquet("test", pq_fixture("with_name.parquet"), 32, track.prefix = "w_", columns = c("a", "b"))
    expect_identical(created, c("w_a", "w_b"))
    expect_identical(track_values("w_a"), expected_values("a"))
})

test_that("aliases resolve; unknown chromosomes are reported", {
    db <- local_pq_db()
    gtrack.import_parquet("test", pq_fixture("alias.parquet"), 32, track.prefix = "al_")
    expect_identical(track_values("al_a"), expected_values("a"))
    expect_identical(track_values("al_b"), expected_values("b"))

    pq_file <- pq_fixture("unknown_chrom.parquet")
    pq_msg <- tryCatch(
        gtrack.import_parquet("test", pq_file, 32, track.prefix = "u_", columns = "a"),
        message = function(m) conditionMessage(m)
    )
    expect_match(pq_msg, "scaffold_9")

    # the message gtrack.import gives for the same rows in a TSV, but for the file name
    tsv <- withr::local_tempfile(fileext = ".tsv")
    rows <- data.frame(
        chrom = c(rep("chr1", 32), rep("scaffold_9", 3)),
        start = c(seq(0, 992, by = 32), 0, 10, 20),
        end = c(pmin(seq(32, 1024, by = 32), 1000), 10, 20, 30),
        value = c(0:31 + 0.25, 1, 1, 1)
    )
    write.table(rows, tsv, sep = "\t", quote = FALSE, row.names = FALSE)
    tsv_msg <- tryCatch(gtrack.import("u_tsv", "test", tsv, 32), message = function(m) conditionMessage(m))
    expect_identical(sub(pq_file, "FILE", pq_msg, fixed = TRUE), sub(tsv, "FILE", tsv_msg, fixed = TRUE))

    # the message was raised after the import completed
    suppressMessages(gtrack.import_parquet("test", pq_file, 32, track.prefix = "u2_"))
    expect_identical(track_values("u2_a")[1:32], expected_values("a")[1:32])
    expect_true(all(is.nan(track_values("u2_a")[33:51])))
})

test_that("unsupported files error and create nothing", {
    db <- local_pq_db()
    cases <- list(
        list(file = "zstd.parquet", binsize = 32, msg = "ZSTD"),
        list(file = "delta.parquet", binsize = 32, msg = "DELTA_BINARY_PACKED"),
        list(file = "overlap.parquet", binsize = 32, msg = "overlap"),
        list(file = "beyond_end.parquet", binsize = 32, msg = "exceeds"),
        list(file = "pandas_default.parquet", binsize = 0, msg = "Parquet import creates dense tracks; binsize must be > 0")
    )
    for (case in cases) {
        expect_error(
            gtrack.import_parquet("test", pq_fixture(case$file), case$binsize, track.prefix = "x_"),
            case$msg,
            fixed = TRUE,
            info = case$file
        )
        expect_length(gtrack.ls(), 0)
        expect_false(dir.exists(file.path(db, "tracks", "x_a.track")))
        expect_false(dir.exists(file.path(db, "tracks", "x_b.track")))
        expect_length(staging_leftovers(db), 0)
    }
    expect_error(gtrack.import_parquet("test", pq_fixture("no_such_*.parquet"), 32), "No files")
})

test_that("arguments are checked before anything is read", {
    db <- local_pq_db()
    pq_path <- pq_fixture("pandas_default.parquet")
    expect_error(gtrack.import_parquet("test", c(pq_path, pq_path), 32), "path must be a single string")
    for (binsize in list(Inf, NaN, NA, -32, 1.5, 2^32, "32", c(32, 64))) {
        expect_error(gtrack.import_parquet("test", pq_path, binsize), "binsize must be a positive integer", info = format(binsize))
    }
    expect_error(gtrack.import_parquet("test", pq_path, 32, defval = "x"), "defval must be a number")
    expect_no_tracks_left(db, c("a", "b"))
})

test_that("integer annotations are honoured; decimals, dates and timestamps are refused", {
    db <- local_pq_db()
    pq_path <- pq_fixture("annotated.parquet")
    info <- .gcall("parquet_columns", pq_path, pq_path, .misha_env())
    expect_identical(
        setNames(info$type, info$name),
        c(
            chrom = "string", start = "int64", end = "int64", u32 = "int32", u64 = "int64", i8 = "int32",
            dec = "other", ts = "other", d = "other"
        )
    )

    # columns = NULL: the decimal, timestamp and date columns are not numbers
    expect_error(gtrack.import_parquet("test", pq_path, 32, track.prefix = "an_"), "dec, ts, d .*columns")
    expect_no_tracks_left(db, paste0("an_", c("u32", "u64", "i8")))

    # unsigned values past the signed range; k = start / 32
    created <- gtrack.import_parquet("test", pq_path, 32, track.prefix = "an_", columns = c("u32", "u64", "i8"))
    expect_identical(created, c("an_u32", "an_u64", "an_i8"))
    expect_identical(track_values("an_u32"), 3e9 + bin_k * 256)
    expect_identical(track_values("an_u64"), 2^63 + bin_k * 2^40)
    expect_identical(track_values("an_i8"), -as.numeric(bin_k))

    for (value_col in c("dec", "ts", "d")) {
        expect_error(
            gtrack.import_parquet("test", pq_path, 32, track.prefix = "bad_", columns = value_col),
            c(dec = "DECIMAL", ts = "TIMESTAMP", d = "DATE")[[value_col]],
            info = value_col
        )
        expect_false(gtrack.exists(paste0("bad_", value_col)))
    }
    expect_length(staging_leftovers(db), 0)
})

test_that("all files must have the same value columns", {
    db <- local_pq_db()
    src <- local_fixture_dir(list(shard_1.parquet = "shard_chr1.parquet", shard_2.parquet = "extra_column.parquet"))
    expect_error(
        gtrack.import_parquet("test", file.path(src, "shard_*.parquet"), 32, track.prefix = "x_"),
        "shard_2.parquet has value columns a, b, c"
    )
    expect_error(
        gtrack.import_parquet("test", file.path(src, "shard_*.parquet"), 32, track.prefix = "x_", columns = c("a", "b")),
        "shard_2.parquet has value columns a, b, c"
    )
    expect_no_tracks_left(db, c("x_a", "x_b", "x_c"))
})

test_that("error messages name the file and the row group", {
    db <- local_pq_db()
    expect_error(
        gtrack.import_parquet("test", pq_fixture("overlap.parquet"), 32),
        "overlap.parquet, row group 0"
    )
    expect_error(
        gtrack.import_parquet("test", pq_fixture("zstd.parquet"), 32),
        "zstd.parquet, row group 0, column chrom"
    )
})

test_that("existing or invalid track names error before anything is written", {
    db <- local_pq_db()
    gtrack.import_parquet("test", pq_fixture("pandas_default.parquet"), 32, track.prefix = "t_x_", columns = "a")
    expect_error(
        gtrack.import_parquet("test", pq_fixture("pandas_default.parquet"), 32, track.prefix = "t_x_"),
        "already exists"
    )
    expect_false(gtrack.exists("t_x_b"))

    expect_error(
        gtrack.import_parquet("test", pq_fixture("pandas_default.parquet"), 32, track.prefix = "nodir."),
        "does not exist"
    )
    expect_identical(gtrack.ls(), "t_x_a")
    expect_length(staging_leftovers(db), 0)
})

test_that("a column that is not a valid track name is an error, never renamed", {
    db <- local_pq_db()
    # refused before anything is written: track2path() would also refuse the name, but only
    # once the tracks have been imported, with no word about 'columns'
    msg <- tryCatch(
        gtrack.import_parquet("test", pq_fixture("bad_name.parquet"), 32, track.prefix = "bn_"),
        error = function(e) conditionMessage(e)
    )
    expect_match(msg, "\"x * 2\"", fixed = TRUE)
    expect_match(msg, "columns =", fixed = TRUE)
    expect_length(gtrack.ls(), 0)
    expect_length(staging_leftovers(db), 0)

    created <- gtrack.import_parquet("test", pq_fixture("bad_name.parquet"), 32, track.prefix = "bn_", columns = c("a", "b"))
    expect_identical(created, c("bn_a", "bn_b"))
    expect_identical(track_values("bn_b"), expected_values("b"))
})

test_that("all or nothing when the second file fails", {
    db <- local_pq_db()

    # a truncated file: refused before anything is staged
    src <- local_fixture_dir(list(shard_1.parquet = "shard_chr1.parquet"))
    pq_bytes <- readBin(pq_fixture("shard_chr2.parquet"), "raw", 1e6)
    writeBin(pq_bytes[seq_len(length(pq_bytes) %/% 2)], file.path(src, "shard_2.parquet"))
    expect_error(gtrack.import_parquet("test", file.path(src, "shard_*.parquet"), 32, track.prefix = "x_"), "is not a Parquet file, or is truncated")
    expect_no_tracks_left(db, c("x_a", "x_b"))

    # a corrupt page in the second file: fails inside the import, after chr1 was written
    src <- local_fixture_dir(list(shard_1.parquet = "shard_chr1.parquet", shard_2.parquet = "corrupt_page.parquet"))
    expect_error(gtrack.import_parquet("test", file.path(src, "shard_*.parquet"), 32, track.prefix = "x_"), "shard_2.parquet, row group 0, column a: corrupt")
    expect_no_tracks_left(db, c("x_a", "x_b"))
})

test_that("all or nothing when a track fails after it was renamed into place", {
    db <- local_pq_db()
    set_attrs <- .gtrack_set_created_attrs
    local_mocked_bindings(.gtrack_set_created_attrs = function(trackstr, ...) {
        if (trackstr == "x_b") stop("mocked attribute failure")
        set_attrs(trackstr, ...)
    })
    expect_error(gtrack.import_parquet("test", pq_fixture("pandas_default.parquet"), 32, track.prefix = "x_"), "mocked attribute failure")
    expect_no_tracks_left(db, c("x_a", "x_b"))
})

test_that("a track the cleanup cannot remove is reported", {
    db <- local_pq_db()
    local_mocked_bindings(
        .gtrack_set_created_attrs = function(...) stop("mocked attribute failure"),
        .gdb.trash = function(...) FALSE
    )
    expect_warning(
        expect_error(gtrack.import_parquet("test", pq_fixture("pandas_default.parquet"), 32, track.prefix = "x_", columns = "a"), "mocked attribute failure"),
        "post-create cleanup left residue"
    )
})

test_that("binning equals gtrack.import of the same rows", {
    db <- local_pq_db()

    # rows that are not bins (the general path), with a gap, a NaN value and short rows, so
    # that some bins are only partly covered
    for (defval in c(NaN, 0.1)) {
        prefix <- if (is.nan(defval)) "un_" else "und_"
        gtrack.import_parquet("test", pq_fixture("unaligned.parquet"), 32, track.prefix = prefix, defval = defval)
        tsv_track <- paste0(prefix, "tsv")
        gtrack.import(tsv_track, "test", pq_fixture("unaligned.tsv"), 32, defval = defval)
        expect_identical(track_values(paste0(prefix, "value")), track_values(tsv_track))
    }

    # rows that are exactly the bins (the fast path), with NaN values and a default value
    rows <- data.frame(
        chrom = rep(c("chr1", "chr2"), c(32, 19)),
        start = c(seq(0, 992, by = 32), seq(0, 576, by = 32)),
        end = c(pmin(seq(32, 1024, by = 32), 1000), pmin(seq(32, 608, by = 32), 600))
    )
    defvals <- c(NaN, 0.1, 3)
    for (i in seq_along(defvals)) {
        defval <- defvals[i]
        for (fx in c("pandas_default", "nulls")) {
            prefix <- sprintf("f%d_%s_", i, fx)
            gtrack.import_parquet("test", pq_fixture(paste0(fx, ".parquet")), 32, track.prefix = prefix, columns = "a", defval = defval)
            value <- expected_values("a")
            if (fx == "nulls") {
                value[bin_k %% 7 == 0] <- NaN
            }
            tsv <- withr::local_tempfile(fileext = ".tsv")
            write.table(cbind(rows, value = value), tsv, sep = "\t", quote = FALSE, row.names = FALSE, na = "NaN")
            tsv_track <- paste0(prefix, "tsv")
            gtrack.import(tsv_track, "test", tsv, 32, defval = defval)
            expect_identical(track_values(paste0(prefix, "a")), track_values(tsv_track), info = prefix)
        }
    }
})

test_that("attributes and indexed database", {
    db <- local_pq_db()
    gtrack.import_parquet("my description", pq_fixture("pandas_default.parquet"), 32, track.prefix = "at_")
    expect_identical(gtrack.attr.get("at_a", "description"), "my description")
    expect_match(gtrack.attr.get("at_b", "created.by"), "gtrack.import_parquet", fixed = TRUE)

    gdb.convert_to_indexed(force = TRUE, verbose = FALSE)
    gtrack.import_parquet("test", pq_fixture("pandas_default.parquet"), 32, track.prefix = "ix_")
    expect_true(file.exists(file.path(db, "tracks", "ix_a.track", "track.idx")))
    expect_identical(track_values("ix_a"), expected_values("a"))
    expect_identical(track_values("ix_b"), expected_values("b"))
})

# Six chromosomes of 5000 bp, a dense track "src" binned at 50 bp, and 20 expressions of it
# with their column names: with gmax.processes = 4 the import runs in kids, 6 x 2 units of work
local_pq6_db <- function(env = parent.frame()) {
    local_db_state(env)
    tmp <- withr::local_tempdir(.local_envir = env)
    db <- file.path(tmp, "pq6db")
    create_test_db(db, data.frame(chrom = paste0("chr", 1:6), size = 5000))
    gdb.init(db)
    bins <- gintervals.all()
    bins <- data.frame(chrom = rep(bins$chrom, each = 100), start = rep(seq(0, 4950, by = 50), 6))
    bins$end <- bins$start + 50
    gtrack.create_dense("src", "src", bins, sin(seq_len(nrow(bins))) * 100, 50)
    db
}

pq6_exprs <- paste0("src * ", 1:20, " + ", 1:20)
pq6_cols <- paste0("c", 1:20)

test_that("parallel import equals serial import", {
    db <- local_pq6_db()
    pq_file <- file.path(withr::local_tempdir(), "x.parquet")
    gextract(pq6_exprs, intervals = gintervals.all(), iterator = 50, colnames = pq6_cols, file = pq_file)

    withr::with_options(list(gmultitasking = FALSE), gtrack.import_parquet("test", pq_file, 50, track.prefix = "ser_"))
    t0 <- proc.time()
    withr::with_options(list(gmultitasking = TRUE, gmax.processes = 4), gtrack.import_parquet("test", pq_file, 50, track.prefix = "par_"))
    dt <- proc.time() - t0
    expect_identical(.glifetime_counters(), c(ref_count = 0L, protect_count = 0L))

    # kids were forked: their CPU time shows up once they are reaped
    expect_gt(dt[["user.child"]] + dt[["sys.child"]], 0)
    for (col in pq6_cols) {
        ser <- gextract(paste0("ser_", col), gintervals.all(), iterator = 50, colnames = "v")
        expect_identical(gextract(paste0("par_", col), gintervals.all(), iterator = 50, colnames = "v"), ser, info = col)
    }
    # and the values are those of the expressions, not NaN in both
    expect_equal(
        gextract("par_c3", gintervals.all(), iterator = 50)$par_c3,
        gextract("src * 3 + 3", gintervals.all(), iterator = 50)[["src * 3 + 3"]],
        tolerance = 1e-6
    )
})

test_that("a kid error removes everything", {
    db <- local_pq6_db()
    src <- withr::local_tempdir()
    # written by one process, so that shard_1.parquet has one row group per chromosome
    withr::with_options(list(gmultitasking = FALSE), {
        gextract(pq6_exprs, intervals = gintervals.all(), iterator = 50, colnames = pq6_cols, file = file.path(src, "shard_1.parquet"))
        # one row on chr5 again: the rows of chr5 overlap, the other chromosomes import fine
        gextract(pq6_exprs, intervals = gintervals("chr5", 1000, 1050), iterator = 50, colnames = pq6_cols, file = file.path(src, "shard_2.parquet"))
    })

    withr::with_options(list(gmultitasking = TRUE, gmax.processes = 4), {
        expect_error(
            gtrack.import_parquet("test", file.path(src, "shard_*.parquet"), 50, track.prefix = "k_"),
            "shard_1.parquet, row group 4 and file .*shard_2.parquet, row group 0: rows overlap"
        )
    })
    expect_identical(.glifetime_counters(), c(ref_count = 0L, protect_count = 0L))
    expect_identical(gtrack.ls(), "src")
    expect_false(any(dir.exists(file.path(db, "tracks", paste0("k_", pq6_cols, ".track")))))
    expect_length(staging_leftovers(db), 0)
})

test_that("gtrack.import reads a one-column Parquet file", {
    local_pq_db()
    gtrack.import("pq_one", "test", pq_fixture("single.parquet"), binsize = 32)
    expect_equal(track_values("pq_one"), expected_values("a"))
    expect_false(is.null(gtrack.attr.get("pq_one", "created.by")) || !nzchar(gtrack.attr.get("pq_one", "created.by")))
})

test_that("gtrack.import refuses a multi-column Parquet file", {
    db <- local_pq_db()
    expect_error(
        gtrack.import("pq_two", "test", pq_fixture("pandas_default.parquet"), binsize = 32),
        "gtrack.import_parquet"
    )
    expect_no_tracks_left(db, "pq_two")
})

test_that("gtrack.import needs binsize > 0 for Parquet", {
    local_pq_db()
    expect_error(gtrack.import("pq_one", "test", pq_fixture("single.parquet"), binsize = 0), "dense tracks")
    expect_error(gtrack.import("pq_one", "test", pq_fixture("single.parquet")), "dense tracks")
})

test_that("a .parquet name is not needed", {
    local_pq_db()
    src <- local_fixture_dir(list(x.dat = "single.parquet"))
    gtrack.import("pq_dat", "test", file.path(src, "x.dat"), binsize = 32)
    expect_equal(track_values("pq_dat"), expected_values("a"))
})

test_that("messages about a gzipped Parquet file name it, not its unzipped copy", {
    db <- local_pq_db()
    src <- withr::local_tempdir()
    gzip_copy <- function(from, to) {
        con <- gzfile(to, "wb")
        writeBin(readBin(from, "raw", file.size(from)), con)
        close(con)
        to
    }

    # an error of the reader (decoding a column) and one of the import (before decoding)
    gz <- gzip_copy(pq_fixture("bad_def_level.parquet"), file.path(src, "bad_def_level.parquet.gz"))
    expect_error(suppressMessages(gtrack.import("g1", "test", gz, binsize = 32)), paste0("File ", gz, ", row group 0, column a: corrupt data"), fixed = TRUE)
    gz <- gzip_copy(pq_fixture("more_rows_than_bases.parquet"), file.path(src, "more_rows_than_bases.parquet.gz"))
    expect_error(suppressMessages(gtrack.import("g2", "test", gz, binsize = 32)), paste0("File ", gz, ", row group 0: 1601 rows"), fixed = TRUE)
    expect_no_tracks_left(db, c("g1", "g2"))

    # the unknown chromosome report: a file of another database, with a chromosome this one lacks
    other <- file.path(src, "otherdb")
    create_test_db(other, data.frame(chrom = c("chr1", "scaffold_9"), size = c(1000, 100)))
    gdb.init(other)
    pq_file <- file.path(src, "x.parquet")
    gextract("1", gintervals.all(), iterator = 32, colnames = "a", file = pq_file)
    gdb.init(db)
    gz <- gzip_copy(pq_file, file.path(src, "x.parquet.gz"))
    msgs <- character(0)
    withCallingHandlers(
        gtrack.import("g3", "test", gz, binsize = 32),
        message = function(m) {
            msgs <<- c(msgs, conditionMessage(m))
            invokeRestart("muffleMessage")
        }
    )
    expect_match(msgs, paste("chromosome name(s) in", gz, "do not exist"), fixed = TRUE, all = FALSE)
    expect_equal(track_values("g3")[1:32], rep(1, 32))
})

test_that("a Parquet input with no rows is an error naming the file", {
    db <- local_pq_db()
    empty <- pq_fixture("no_rowgroups.parquet")
    expect_error(gtrack.import_parquet("test", empty, 32, track.prefix = "e_"), paste(empty, "has no rows"), fixed = TRUE)
    expect_no_tracks_left(db, "e_a")
    expect_error(gtrack.import("e_one", "test", empty, binsize = 32), paste(empty, "has no rows"), fixed = TRUE)
    expect_no_tracks_left(db, "e_one")

    # misha's own output for an empty scope: one row group of zero rows
    own <- file.path(withr::local_tempdir(), "own.parquet")
    gtrack.create_sparse("dense_src", "src", gintervals("chr1", 0, 10), 1)
    gextract("dense_src", intervals = gintervals("chr1", 0, 100), iterator = gintervals("chr2", 0, 100), file = own)
    expect_error(gtrack.import_parquet("test", own, 32, track.prefix = "o_"), paste(own, "has no rows"), fixed = TRUE)
    expect_length(staging_leftovers(db), 0)
    expect_identical(gtrack.ls(), "dense_src")
})

test_that("an empty shard among non-empty shards imports", {
    db <- local_pq_db()
    src <- local_fixture_dir(list(shard_1.parquet = "no_rowgroups.parquet", shard_2.parquet = "single.parquet"))
    gtrack.import_parquet("test", file.path(src, "shard_*.parquet"), 32, track.prefix = "x_")
    expect_identical(track_values("x_a"), expected_values("a"))
    # all shards empty: an error naming all of them
    src <- local_fixture_dir(list(s1.parquet = "no_rowgroups.parquet", s2.parquet = "no_rowgroups.parquet"))
    expect_error(gtrack.import_parquet("test", file.path(src, "s*.parquet"), 32, track.prefix = "y_"), "s1.parquet, .*s2.parquet have no rows")
    expect_length(staging_leftovers(db), 0)
})

test_that("a directory as path suggests the glob", {
    db <- local_pq_db()
    src <- local_fixture_dir(list(a.parquet = "single.parquet"))
    expect_error(
        gtrack.import_parquet("test", src, 32, track.prefix = "d_"),
        paste0("is a directory; to import the Parquet files in it use path = \"", file.path(src, "*.parquet"), "\""),
        fixed = TRUE
    )
    expect_no_tracks_left(db, "d_a")
    expect_error(gtrack.import_parquet("test", file.path(src, "nothing"), 32), "No files match")
})
