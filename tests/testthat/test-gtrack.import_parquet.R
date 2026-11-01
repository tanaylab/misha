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
    info <- .gcall("parquet_columns", pq_path, .misha_env())
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
