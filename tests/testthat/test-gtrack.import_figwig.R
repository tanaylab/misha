# gtrack.import_figwig: the tracks of figwig files as misha tracks. The fixtures were written by the
# figwig package itself (fixtures/figwig/make_fixtures.py, which also says what each holds), and
# the expected values are those figwig's own reader reads from them (*.values.tsv).

restore_groot_on_exit()

fw_fixture <- function(name) {
    normalizePath(test_path("fixtures", "figwig", name), mustWork = TRUE)
}

FW_SIZES <- c(chr1 = 300, chr2 = 100, chr3 = 64)

local_fw_db <- function(sizes = FW_SIZES, indexed = FALSE, env = parent.frame()) {
    local_db_state(env)
    tmp <- withr::local_tempdir(.local_envir = env)
    db <- file.path(tmp, "fwdb")
    create_test_db(db, data.frame(chrom = names(sizes), size = as.numeric(sizes)))
    gdb.init(db)
    if (indexed) {
        gdb.convert_to_indexed(force = TRUE, verbose = FALSE)
    }
    db
}

fw_values <- function(name) {
    v <- utils::read.delim(fw_fixture(name), stringsAsFactors = FALSE, colClasses = c("character", "character", "numeric", "numeric", "character"))
    v$value <- as.numeric(v$value)
    v
}

# The values figwig reads for a track on a chromosome of 'size' bases: 0 where no run says otherwise
base_values <- function(values, track, chrom, size) {
    x <- numeric(size)
    r <- values[values$track == track & values$chrom == chrom, ]
    for (i in seq_len(nrow(r))) {
        x[(r$start[i] + 1):r$end[i]] <- r$value[i]
    }
    x
}

# x rounded to float, as a track file holds it
f32 <- function(x) {
    readBin(writeBin(as.numeric(x), raw(), size = 4), "numeric", size = 4, n = length(x))
}

# The intervals of the sparse track of x: maximal runs of equal float values other than 0 and NaN
expected_sparse <- function(x) {
    f <- f32(x)
    if (!length(f)) {
        return(data.frame(start = numeric(0), end = numeric(0), value = numeric(0)))
    }
    r <- rle(ifelse(is.nan(f), NA, f))
    ends <- cumsum(r$lengths)
    starts <- ends - r$lengths
    keep <- !is.na(r$values) & r$values != 0
    data.frame(start = starts[keep], end = ends[keep], value = r$values[keep])
}

sparse_of <- function(track, chrom, size) {
    d <- gextract(track, gintervals(chrom, 0, size))
    if (is.null(d)) {
        return(data.frame(start = numeric(0), end = numeric(0), value = numeric(0)))
    }
    data.frame(start = d$start, end = d$end, value = d[[track]])
}

# The bins of the dense track of x (the values of the first length(x) bases of a chromosome of
# 'size'): the mean of the non-NaN values, the other bases counting as defval unless it is NaN
expected_dense <- function(x, size, binsize, defval) {
    nb <- ceiling(size / binsize)
    vapply(seq_len(nb) - 1, function(k) {
        lo <- k * binsize
        hi <- min(lo + binsize, size)
        v <- if (lo < length(x)) x[(lo + 1):min(hi, length(x))] else numeric(0)
        v <- v[!is.nan(v)]
        if (!is.nan(defval)) {
            v <- c(v, rep(defval, (hi - lo) - length(v)))
        }
        if (length(v)) f32(sum(v) / length(v)) else NaN
    }, numeric(1))
}

dense_of <- function(track, chrom, size) {
    gextract(track, gintervals(chrom, 0, size))[[track]]
}

expect_close <- function(actual, expected, info = NULL) {
    expect_equal(length(actual), length(expected), info = info)
    same <- (is.nan(actual) & is.nan(expected)) | (!is.nan(actual) & !is.nan(expected) & abs(actual - expected) <= 2e-7 * abs(expected))
    expect_true(all(same), info = info)
}

FORMAT_TRACKS <- c(
    "ones", "random", "random2", "runs", "runs2", "neg8", "neg16", "neg32", "u16", "u32", "wide", "signal",
    "dense16", "runs_int", "empty", "dense_f32", "dense_i16", "dense_f64", "cell0", "cell1", "cell30", "cell90",
    "site0", "site60", "site119"
)

test_that("every layout of the format imports the values figwig reads, as sparse tracks", {
    db <- local_fw_db()
    values <- fw_values("zlib.values.tsv")
    for (fx in c("zlib", "stored", "tracks", "tracks_zlib", "appended")) {
        prefix <- paste0(fx, "_")
        tracks <- if (fx == "appended") c(FORMAT_TRACKS, "late") else FORMAT_TRACKS
        created <- gtrack.import_figwig(fw_fixture(paste0(fx, ".fw")), tracks = tracks, track.prefix = prefix)
        expect_identical(created, paste0(prefix, tracks))
        for (track in tracks) {
            expect_identical(gtrack.info(paste0(prefix, track))$type, "sparse")
            for (chrom in names(FW_SIZES)) {
                x <- if (track == "late") {
                    replace(numeric(FW_SIZES[[chrom]]), 66, if (chrom == "chr1") 9 else 0)
                } else {
                    base_values(values, track, chrom, FW_SIZES[[chrom]])
                }
                expect_equal(sparse_of(paste0(prefix, track), chrom, FW_SIZES[[chrom]]), expected_sparse(x), info = paste(fx, track, chrom))
            }
        }
    }
})

test_that("dense tracks are the import mean of the values figwig reads", {
    db <- local_fw_db()
    values <- fw_values("zlib.values.tsv")
    tracks <- c("random", "runs", "neg32", "u32", "wide", "signal", "dense_f32", "dense_f64", "cell30", "empty")
    for (fx in c("zlib", "tracks")) {
        for (args in list(list(7, NaN), list(64, 0.5), list(1, NaN), list(1000, -1))) {
            prefix <- sprintf("%s_%d_", fx, args[[1]])
            gtrack.import_figwig(fw_fixture(paste0(fx, ".fw")), binsize = args[[1]], defval = args[[2]], tracks = tracks, track.prefix = prefix)
            for (track in tracks) {
                expect_equal(gtrack.info(paste0(prefix, track))$bin.size, args[[1]])
                for (chrom in names(FW_SIZES)) {
                    size <- FW_SIZES[[chrom]]
                    expect_close(
                        dense_of(paste0(prefix, track), chrom, size),
                        expected_dense(base_values(values, track, chrom, size), size, args[[1]], args[[2]]),
                        info = paste(fx, args[[1]], track, chrom)
                    )
                }
            }
        }
    }
})

test_that("dense bins are those gtrack.import makes from a bedGraph of the same values", {
    db <- local_fw_db()
    values <- fw_values("zlib.values.tsv")
    for (track in c("signal", "dense_f32", "neg16")) {
        rows <- do.call(rbind, lapply(names(FW_SIZES), function(chrom) {
            x <- base_values(values, track, chrom, FW_SIZES[[chrom]])
            keep <- !is.nan(x)
            data.frame(chrom = chrom, start = which(keep) - 1, end = which(keep), value = sprintf("%.17g", x[keep]))
        }))
        bg <- withr::local_tempfile(fileext = ".bedgraph")
        utils::write.table(format(rows, scientific = FALSE, trim = TRUE), bg, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
        for (defval in c(NaN, 0)) {
            fw_track <- sprintf("fw_%s_%s", track, if (is.nan(defval)) "nan" else "zero")
            bg_track <- sprintf("bg_%s_%s", track, if (is.nan(defval)) "nan" else "zero")
            gtrack.import_figwig(fw_fixture("zlib.fw"), binsize = 9, defval = defval, tracks = track, names = fw_track)
            gtrack.import(bg_track, "", bg, binsize = 9, defval = defval)
            for (chrom in names(FW_SIZES)) {
                expect_identical(dense_of(fw_track, chrom, FW_SIZES[[chrom]]), dense_of(bg_track, chrom, FW_SIZES[[chrom]]), info = paste(track, chrom))
            }
        }
    }
})

test_that("blocks of track order, and chunks written from entries", {
    db <- local_fw_db(c(chrA = 70000, chrB = 10))
    values <- fw_values("blocks.values.tsv")
    gtrack.import_figwig(fw_fixture("blocks.fw"), track.prefix = "s_")
    gtrack.import_figwig(fw_fixture("blocks.fw"), track.prefix = "d_", binsize = 1000, defval = 0)
    ones <- c(rep(1, 66000), numeric(4000))
    for (chrom in c("chrA", "chrB")) {
        size <- c(chrA = 70000, chrB = 10)[[chrom]]
        for (track in c("t", "u", "ones")) {
            x <- if (track == "ones") (if (chrom == "chrA") ones else numeric(10)) else base_values(values, track, chrom, size)
            expect_equal(sparse_of(paste0("s_", track), chrom, size), expected_sparse(x), info = paste(track, chrom))
            expect_close(dense_of(paste0("d_", track), chrom, size), expected_dense(x, size, 1000, 0), info = paste(track, chrom))
        }
    }
    expect_equal(sparse_of("s_ones", "chrA", 70000), data.frame(start = 0, end = 66000, value = 1))

    db <- local_fw_db()
    values <- fw_values("entries.values.tsv")
    gtrack.import_figwig(fw_fixture("entries.fw"), track.prefix = "e_")
    for (track in c("steps", "random", "gaps", "scattered")) {
        for (chrom in names(FW_SIZES)) {
            x <- base_values(values, track, chrom, FW_SIZES[[chrom]])
            expect_equal(sparse_of(paste0("e_", track), chrom, FW_SIZES[[chrom]]), expected_sparse(x), info = paste(track, chrom))
        }
    }
})

test_that("an indexed database gets indexed tracks", {
    # chr3 of 128 bp: the dummy sequence of a shorter chromosome fails the conversion's check
    sizes <- c(chr1 = 300, chr2 = 100, chr3 = 128)
    db <- local_fw_db(sizes, indexed = TRUE)
    values <- fw_values("zlib.values.tsv")
    expect_warning(gtrack.import_figwig(fw_fixture("zlib.fw"), tracks = c("random", "u16"), track.prefix = "s_"), "differ in length")
    expect_warning(gtrack.import_figwig(fw_fixture("zlib.fw"), tracks = c("random", "u16"), track.prefix = "d_", binsize = 10), "differ in length")
    expect_true(file.exists(file.path(db, "tracks", "s_random.track", "track.idx")))
    expect_true(file.exists(file.path(db, "tracks", "d_u16.track", "track.idx")))
    for (chrom in names(sizes)) {
        size <- sizes[[chrom]]
        x <- base_values(values, "u16", chrom, FW_SIZES[[chrom]])
        expect_equal(sparse_of("s_u16", chrom, size), expected_sparse(x))
        expect_close(dense_of("d_u16", chrom, size), expected_dense(x, size, 10, NaN))
    }
})

test_that("names are sanitized, attributes copied, chromosomes matched", {
    db <- local_fw_db(c(chr1 = 300, chr2 = 100, chr3 = 64))
    values <- fw_values("semantics.values.tsv")
    expect_warning(
        expect_message(
            created <- gtrack.import_figwig(fw_fixture("semantics.fw"), binsize = 10, defval = 7),
            "chrUn_x"
        ),
        "differ in length"
    )
    expect_identical(created, c("atac_plus", "atac_minus", "H3K27ac_rep1", "X1x"))
    for (i in seq_along(created)) {
        fw_name <- c("atac.+", "atac.-", "H3K27ac-rep1", "1x")[i]
        # "1" is chr1 by alias; chr2 is 90 bp in the file; the file has no chr3
        expect_close(dense_of(created[i], "chr1", 300), expected_dense(base_values(values, fw_name, "1", 300), 300, 10, 7))
        expect_close(dense_of(created[i], "chr2", 100), expected_dense(base_values(values, fw_name, "chr2", 90), 100, 10, 7))
        expect_identical(dense_of(created[i], "chr3", 64), rep(7, 7))
    }

    expect_identical(gtrack.attr.get("atac_plus", "description"), "from the file")
    expect_identical(gtrack.attr.get("atac_plus", "nested"), "{\"k\": [1, 2.5], \"s\": \"\\u00e9\"}")
    expect_identical(gtrack.attr.get("atac_plus", "flag"), "true")
    expect_identical(gtrack.attr.get("atac_plus", "rep"), "1")
    expect_identical(gtrack.attr.get("atac_plus", "uni"), "caf\u00e9")
    for (absent in c("none", "empty", "type", "binsize", "name", "dtype", "group")) {
        expect_identical(gtrack.attr.get("atac_plus", absent), "")
    }
    expect_match(gtrack.attr.get("atac_plus", "created.by"), "^gtrack.import_figwig\\(")
    expect_identical(gtrack.attr.get("atac_minus", "description"), "")

    # a description given replaces the file's; attrs win over the file's
    suppressWarnings(suppressMessages(gtrack.import_figwig(fw_fixture("semantics.fw"),
        tracks = "atac.+", track.prefix = "d_",
        description = "mine", attrs = c(assay = "ATAC-seq", lab = "x")
    )))
    expect_identical(gtrack.attr.get("d_atac_plus", "description"), "mine")
    expect_identical(gtrack.attr.get("d_atac_plus", "assay"), "ATAC-seq")
    expect_identical(gtrack.attr.get("d_atac_plus", "lab"), "x")
    expect_equal(
        sparse_of("d_atac_plus", "chr2", 100),
        expected_sparse(base_values(values, "atac.+", "chr2", 90))
    )
})

test_that("tracks are selected by names, indices and attributes", {
    db <- local_fw_db()
    imp <- function(...) suppressWarnings(suppressMessages(gtrack.import_figwig(fw_fixture("semantics.fw"), ...)))
    expect_identical(imp(tracks = c("1x", "atac.-")), c("X1x", "atac_minus"))
    expect_identical(imp(tracks = c(1, 3), track.prefix = "i_"), c("i_atac_plus", "i_H3K27ac_rep1"))
    expect_identical(imp(select = list(assay = "ATAC"), track.prefix = "s_"), c("s_atac_plus", "s_atac_minus"))
    # compared as text: 2 matches 2 and "2"
    expect_identical(imp(select = list(rep = 2), track.prefix = "r_"), c("r_atac_minus", "r_H3K27ac_rep1"))
    expect_identical(imp(select = list(assay = c("ChIP", "x"), rep = "2"), track.prefix = "m_"), "m_H3K27ac_rep1")
    expect_identical(imp(tracks = c("1x", "atac.+", "atac.-"), select = list(assay = "ATAC"), track.prefix = "t_"), c("t_atac_plus", "t_atac_minus"))
    expect_identical(imp(select = list(dtype = "float32"), track.prefix = "dt_"), "dt_atac_minus")
    expect_identical(imp(tracks = c("atac.+", "1x"), names = c("a", "b")), c("a", "b"))
    expect_identical(imp(tracks = c("atac.+", "1x"), names = c("1x" = "one"), track.prefix = "n_"), c("n_atac_plus", "n_one"))

    expect_error(imp(tracks = "nope"), "no track")
    expect_error(imp(tracks = 5), "out of range")
    expect_error(imp(tracks = c(1, 1)), "more than once")
    expect_error(imp(select = list(assay = "DNase")), "matches")
    expect_error(imp(select = list(cell = "K562")), "has the attribute")
})

test_that("bad names and existing tracks create nothing", {
    db <- local_fw_db()
    imp <- function(...) suppressWarnings(suppressMessages(gtrack.import_figwig(fw_fixture("semantics.fw"), ...)))
    imp(tracks = "1x", names = "taken")
    before <- gtrack.ls()
    expect_error(imp(tracks = c("atac.+", "1x"), names = c("fresh", "taken")), "already exists")
    expect_error(imp(tracks = c("atac.+", "1x"), names = c("x", "x")), "same track name")
    expect_error(imp(tracks = "1x", names = "1y"), "Invalid track name")
    expect_error(imp(tracks = "1x", track.prefix = "nodir."), "does not exist")
    expect_error(imp(tracks = "1x", names = c("a", "b")), "element")
    expect_identical(gtrack.ls(), before)
    entries <- list.files(file.path(db, "tracks"), all.files = TRUE, no.. = TRUE)
    expect_identical(entries[grepl("\\.tmp\\.", entries) & !grepl("^\\.trash\\.", entries)], character(0))
})

test_that("files that are not figwig, or corrupt, are refused", {
    db <- local_fw_db()
    not_fw <- withr::local_tempfile(fileext = ".fw")
    writeLines(strrep("hello", 40), not_fw)
    expect_error(gtrack.import_figwig(not_fw), "not a figwig file")
    expect_error(gtrack.import_figwig(file.path(tempdir(), "none.fw")), "does not exist")

    bytes <- readBin(fw_fixture("zlib.fw"), "raw", file.size(fw_fixture("zlib.fw")))
    corrupt <- function(at) {
        f <- withr::local_tempfile(fileext = ".fw", .local_envir = parent.frame())
        b <- bytes
        b[at] <- xor(b[at], as.raw(0x40))
        writeBin(b, f)
        f
    }
    expect_error(gtrack.import_figwig(corrupt(9)), "corrupt") # in the header
    expect_error(gtrack.import_figwig(corrupt(length(bytes) - 200)), "corrupt") # in the tables
    truncated <- withr::local_tempfile(fileext = ".fw")
    writeBin(bytes[seq_len(length(bytes) - 100)], truncated)
    expect_error(gtrack.import_figwig(truncated), "corrupt")
    # in the data region: a chunk's zlib stream or stored CRC-32
    expect_error(gtrack.import_figwig(corrupt(1000), tracks = FORMAT_TRACKS), "corrupt")
    expect_length(gtrack.ls(), 0)
})

test_that("names are sanitized as PyMisha sanitizes them", {
    expect_identical(
        .figwig_sanitize(c("atac.+", "atac.-", "H3K27ac-rep1", "ENCFF001 (rep 2)", "1abc", "+", enc2utf8("café.x")), FALSE),
        c("atac_plus", "atac_minus", "H3K27ac_rep1", "ENCFF001_rep_2_", "X1abc", "plus", "caf_x")
    )
    expect_identical(.figwig_sanitize("1abc", TRUE), "1abc")
})
