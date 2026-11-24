fixture_dir <- testthat::test_path("fixtures", "motifs")
export_dir <- testthat::test_path("fixtures", "motifs", "export")

# The golden files in fixtures/motifs/export were written by these functions;
# PyMisha's tests/test_motif_export.py checks that it writes the same bytes from
# the same inputs (the reader fixtures, random_*.tsv and edge_motifs()).

acgt <- c("A", "C", "G", "T")

# Random motifs with metadata, stored as hex floats so R and Python read the
# same doubles.
random_motifs <- function() {
    mt <- utils::read.delim(file.path(export_dir, "random_matrices.tsv"), colClasses = "character")
    meta <- utils::read.delim(file.path(export_dir, "random_meta.tsv"),
        colClasses = "character", na.strings = character(0)
    )
    num <- function(x) if (x == "") NA_real_ else as.numeric(x)
    res <- lapply(meta$id, function(id) {
        rows <- mt[mt$id == id, ]
        m <- matrix(as.numeric(unlist(rows[, acgt])), ncol = 4, dimnames = list(NULL, acgt))
        row <- meta[meta$id == id, ]
        if (row$name != "") attr(m, "name") <- row$name
        for (key in c("nsites", "E", "log_odds_threshold", "log_p_value")) {
            v <- num(row[[key]])
            if (!is.na(v)) attr(m, key) <- v
        }
        m
    })
    stats::setNames(res, meta$id)
}

# Hand-made edge cases: width 1, zeros, whitespace in names, values that round
edge_motifs <- function() {
    w1 <- matrix(c(0.1, 0.2, 0.3, 0.4), ncol = 4, dimnames = list(NULL, acgt))
    attr(w1, "name") <- "  one\tposition  name "
    attr(w1, "nsites") <- 7.5
    attr(w1, "E") <- 0
    zeros <- matrix(c(
        1, 0, 0, 0,
        0, 0, 0, 1,
        0.5, 0.5, 0, 0,
        0, 0.5, 0, 0.5
    ), ncol = 4, byrow = TRUE, dimnames = list(NULL, acgt))
    thirds <- matrix(c(
        1 / 3, 1 / 3, 1 / 6, 1 / 6,
        1 / 7, 2 / 7, 3 / 7, 1 / 7,
        0.45, 0.4, 0.1, 0.05,
        0.3, 0.3, 0.2, 0.2
    ), ncol = 4, byrow = TRUE, dimnames = list(NULL, acgt))
    attr(thirds, "name") <- "thirds"
    attr(thirds, "log_odds_threshold") <- 1.23456789
    attr(thirds, "url") <- "http://example.org/thirds"
    list(w1 = w1, zeros = zeros, thirds = thirds)
}

# Every golden file and how it is written
golden_cases <- function() {
    meme <- gseq.read_meme(file.path(fixture_dir, "test_motifs.meme"))
    jaspar <- gseq.read_jaspar(file.path(fixture_dir, "test_motifs.jaspar"))
    homer <- gseq.read_homer(file.path(fixture_dir, "test_motifs.homer.motif"))
    rnd <- random_motifs()
    edge <- edge_motifs()
    sets <- list(from_meme = meme, from_jaspar = jaspar, from_homer = homer, random = rnd, edge = edge)
    cases <- list()
    for (nm in names(sets)) {
        local({
            x <- sets[[nm]]
            cases[[paste0(nm, ".meme")]] <<- function(f) gseq.write_meme(x, f)
            cases[[paste0(nm, ".jaspar")]] <<- function(f) gseq.write_jaspar(x, f)
            cases[[paste0(nm, ".motif")]] <<- function(f) gseq.write_homer(x, f)
        })
    }
    cases[["edge_simple.pfm"]] <- function(f) gseq.write_jaspar(edge, f, format = "simple")
    cases[["edge_args.meme"]] <- function(f) {
        gseq.write_meme(edge, f, background = c(T = 0.3, G = 0.2, C = 0.2, A = 0.3), strands = "+")
    }
    cases[["edge_args.jaspar"]] <- function(f) gseq.write_jaspar(edge, f, nsites = c(10, 1000, 3))
    cases[["edge_round.jaspar"]] <- function(f) gseq.write_jaspar(edge, f, nsites = c(10, 1000, 3), round_counts = TRUE)
    cases[["random_round.jaspar"]] <- function(f) gseq.write_jaspar(rnd, f, round_counts = TRUE)
    cases[["edge_args.motif"]] <- function(f) gseq.write_homer(edge, f, threshold_frac = 0.8)
    cases
}

strip_attrs <- function(m) {
    m <- unclass(m)
    attributes(m) <- list(dim = dim(m), dimnames = dimnames(m))
    m
}

expect_same_matrices <- function(a, b, tolerance = 1e-6) {
    expect_equal(names(a), names(b))
    for (i in seq_along(a)) {
        expect_equal(strip_attrs(a[[i]]), strip_attrs(b[[i]]), tolerance = tolerance)
    }
}

suppress_renorm <- function(expr) {
    withCallingHandlers(expr, warning = function(w) {
        if (grepl("re-normalizing", conditionMessage(w))) invokeRestart("muffleWarning")
    })
}

# === Round trips on the reader fixtures ===

test_that("gseq.write_meme round-trips the MEME fixture with its metadata", {
    x <- gseq.read_meme(file.path(fixture_dir, "test_motifs.meme"))
    f <- withr::local_tempfile(fileext = ".meme")
    expect_identical(gseq.write_meme(x, f), f)
    y <- gseq.read_meme(f)
    expect_equal(y, x, tolerance = 1e-12)
    expect_equal(attr(y[[1]], "background"), c(A = 0.29, C = 0.21, G = 0.21, T = 0.29))
    expect_equal(attr(y[[1]], "url"), "http://jaspar.genereg.net/matrix/MA0001.1")
    expect_true(is.na(attr(y[[3]], "E")))
})

test_that("gseq.write_jaspar round-trips the JASPAR fixture with its counts", {
    x <- gseq.read_jaspar(file.path(fixture_dir, "test_motifs.jaspar"))
    f <- withr::local_tempfile(fileext = ".jaspar")
    gseq.write_jaspar(x, f)
    expect_equal(gseq.read_jaspar(f), x, tolerance = 1e-12)
    # integer counts are written as integers, as in the original
    expect_equal(readLines(f)[2], "A  [ 10 12  4  1  0  5 ]")
})

test_that("gseq.write_jaspar format = 'simple' round-trips the simple PFM fixture", {
    x <- gseq.read_jaspar(file.path(fixture_dir, "test_motifs_simple.pfm"))
    d <- withr::local_tempdir()
    f <- file.path(d, "test_motifs_simple.pfm")
    # without nsites (simple files have none), counts are scaled to 100
    gseq.write_jaspar(x, f, format = "simple")
    expect_equal(gseq.read_jaspar(f), x, tolerance = 1e-8)
    expect_equal(readLines(f)[1], "66.666667        80 26.666667  6.666667")
    gseq.write_jaspar(x, f, nsites = 15, format = "simple")
    expect_equal(gseq.read_jaspar(f), x, tolerance = 1e-12)
    expect_equal(
        readLines(f),
        c("10 12  4  1", " 0  0  7  0", " 3  2  1  0", " 2  1  3 14")
    )
})

test_that("gseq.write_homer round-trips the HOMER fixture with its metadata", {
    x <- gseq.read_homer(file.path(fixture_dir, "test_motifs.homer.motif"))
    f <- withr::local_tempfile(fileext = ".motif")
    gseq.write_homer(x, f)
    expect_equal(gseq.read_homer(f), x, tolerance = 1e-12)
    expect_equal(
        readLines(f)[1],
        ">ATGACTCA\tAP1(bZIP)/ThP1-cJun\t8.036341\t-4130.668834"
    )
})

# === Round trips on random matrices ===

test_that("random motifs round-trip through all three formats", {
    x <- random_motifs()
    widths <- vapply(x, nrow, integer(1))
    expect_true(1L %in% widths)

    f <- withr::local_tempfile(fileext = ".meme")
    gseq.write_meme(x, f)
    y <- gseq.read_meme(f)
    expect_equal(names(y), names(x))
    for (id in names(x)) {
        m <- abs(strip_attrs(x[[id]]) - strip_attrs(y[[id]]))
        expect_lte(max(m), 5e-7 + 1e-12)
        expect_equal(attr(y[[id]], "name"), if (is.null(attr(x[[id]], "name"))) NA_character_ else gsub(" +", " ", attr(x[[id]], "name")))
        expect_equal(attr(y[[id]], "nsites"), if (is.null(attr(x[[id]], "nsites"))) NA_real_ else attr(x[[id]], "nsites"), tolerance = 1e-14)
        expect_equal(attr(y[[id]], "E"), if (is.null(attr(x[[id]], "E"))) NA_real_ else attr(x[[id]], "E"), tolerance = 1e-14)
    }

    f <- withr::local_tempfile(fileext = ".jaspar")
    gseq.write_jaspar(x, f)
    y <- gseq.read_jaspar(f)
    for (id in names(x)) {
        ns <- attr(x[[id]], "nsites")
        expect_lte(max(abs(strip_attrs(x[[id]]) - strip_attrs(y[[id]]))), 1.5e-6 / (if (is.null(ns)) 100 else ns))
        expect_equal(attr(y[[id]], "nsites"), if (is.null(ns)) 100 else ns, tolerance = 1e-6)
    }

    f <- withr::local_tempfile(fileext = ".motif")
    gseq.write_homer(x, f)
    y <- suppressWarnings(gseq.read_homer(f)) # random consensuses repeat
    expect_length(y, length(x))
    for (i in seq_along(x)) {
        expect_lte(max(abs(strip_attrs(x[[i]]) - strip_attrs(y[[i]]))), 5e-7 + 1e-12)
        thr <- attr(x[[i]], "log_odds_threshold")
        if (!is.null(thr)) expect_lte(abs(attr(y[[i]], "log_odds_threshold") - thr), 5e-7 + 1e-12)
        lp <- attr(x[[i]], "log_p_value")
        if (is.null(lp)) {
            expect_true(is.na(attr(y[[i]], "log_p_value")))
        } else {
            expect_lte(abs(attr(y[[i]], "log_p_value") - lp), 5e-7 + 1e-9)
        }
    }
})

test_that("a second write of what was read gives the same bytes", {
    # Not JASPAR: reading divides each position by its own count total, which
    # can differ from nsites in the last written digit.
    x <- random_motifs()
    for (fmt in c("meme", "homer")) {
        write <- get(paste0("gseq.write_", fmt))
        read <- get(paste0("gseq.read_", fmt))
        f1 <- withr::local_tempfile()
        f2 <- withr::local_tempfile()
        write(x, f1)
        write(suppressWarnings(read(f1)), f2)
        expect_identical(readLines(f1), readLines(f2), label = fmt)
    }
})

# === Golden files (shared with PyMisha) ===

test_that("the writers produce the golden files byte for byte", {
    cases <- golden_cases()
    for (nm in names(cases)) {
        f <- withr::local_tempfile()
        suppress_renorm(cases[[nm]](f))
        golden <- file.path(export_dir, nm)
        expect_true(file.exists(golden), label = nm)
        expect_identical(
            readBin(f, "raw", file.info(f)$size),
            readBin(golden, "raw", file.info(golden)$size),
            label = nm
        )
    }
})

# === Cross-format conversion ===

test_that("MEME motifs written as JASPAR and HOMER keep IDs, names and nsites", {
    x <- gseq.read_meme(file.path(fixture_dir, "test_motifs.meme"))
    f <- withr::local_tempfile()
    gseq.write_jaspar(x, f)
    y <- gseq.read_jaspar(f)
    expect_same_matrices(x, y, tolerance = 1e-9)
    expect_equal(unname(sapply(y, attr, "name")), c("AGL3", "RUNX1", "TFAP2A"))
    expect_equal(unname(sapply(y, attr, "nsites")), c(20, 15, 18))

    gseq.write_homer(x, f)
    y <- gseq.read_homer(f)
    expect_equal(names(y), c("CAGT", "TAGGTN", "CGATN"))
    expect_equal(unname(sapply(y, attr, "name")), c("AGL3", "RUNX1", "TFAP2A"))
})

test_that("HOMER motifs written as MEME keep consensus IDs and names", {
    x <- gseq.read_homer(file.path(fixture_dir, "test_motifs.homer.motif"))
    f <- withr::local_tempfile()
    gseq.write_meme(x, f)
    y <- gseq.read_meme(f)
    expect_same_matrices(x, y)
    expect_equal(attr(y[[1]], "name"), "AP1(bZIP)/ThP1-cJun")
    expect_true(is.na(attr(y[[1]], "nsites")))
    expect_equal(attr(y[[1]], "background"), c(A = 0.25, C = 0.25, G = 0.25, T = 0.25))
    expect_equal(attr(y[[1]], "strand"), "+ -")
})

# === Inputs ===

test_that("a plain matrix without names or metadata is written with defaults", {
    m <- matrix(c(0.7, 0.1, 0.1, 0.1, 0.25, 0.25, 0.25, 0.25), ncol = 4, byrow = TRUE)
    f <- withr::local_tempfile()
    gseq.write_meme(m, f)
    expect_equal(readLines(f)[10:11], c("MOTIF motif_1", "letter-probability matrix: alength= 4 w= 2"))
    y <- gseq.read_meme(f)
    expect_equal(names(y), "motif_1")
    expect_true(is.na(attr(y[[1]], "name")))
    expect_true(is.na(attr(y[[1]], "nsites")))
    expect_equal(strip_attrs(y[[1]]), matrix(as.numeric(m), ncol = 4, dimnames = list(NULL, acgt)))

    gseq.write_jaspar(m, f)
    expect_equal(readLines(f), c(">motif_1 motif_1", "A  [ 70 25 ]", "C  [ 10 25 ]", "G  [ 10 25 ]", "T  [ 10 25 ]"))

    gseq.write_homer(m, f)
    expect_equal(strsplit(readLines(f)[1], "\t")[[1]][1:2], c(">AN", "motif_1"))
})

test_that("unnamed list elements get motif_<i> IDs", {
    m <- matrix(0.25, nrow = 2, ncol = 4, dimnames = list(NULL, acgt))
    f <- withr::local_tempfile()
    gseq.write_meme(list(m, X = m, m), f)
    expect_equal(names(gseq.read_meme(f)), c("motif_1", "X", "motif_3"))
})

test_that("columns are matched by name and bases-by-positions matrices are transposed", {
    m <- matrix(c(0.7, 0.1, 0.1, 0.1, 0.1, 0.2, 0.3, 0.4), ncol = 4, byrow = TRUE, dimnames = list(NULL, acgt))
    f1 <- withr::local_tempfile()
    f2 <- withr::local_tempfile()
    f3 <- withr::local_tempfile()
    f4 <- withr::local_tempfile()
    gseq.write_meme(list(M = m), f1)
    gseq.write_meme(list(M = m[, c("T", "G", "C", "A")]), f2)
    gseq.write_meme(list(M = t(m)), f3)
    gseq.write_meme(list(M = as.data.frame(m)), f4)
    expect_identical(readLines(f2), readLines(f1))
    expect_identical(readLines(f3), readLines(f1))
    expect_identical(readLines(f4), readLines(f1))
})

test_that("rows that do not sum to 1 are re-normalized with a warning", {
    m <- matrix(c(1, 2, 3, 4, 0, 0, 0, 0, 0.5, 0.5, 0, 0), ncol = 4, byrow = TRUE)
    f <- withr::local_tempfile()
    expect_warning(gseq.write_meme(m, f), "re-normalizing 2 row")
    y <- gseq.read_meme(f)[[1]]
    expect_equal(unname(y[1, ]), c(0.1, 0.2, 0.3, 0.4))
    expect_equal(unname(y[2, ]), rep(0.25, 4))
    expect_equal(unname(y[3, ]), c(0.5, 0.5, 0, 0))
})

test_that("negative zeros are written as zeros", {
    m <- matrix(c(-0, 0.5, 0.5, -0), ncol = 4)
    f <- withr::local_tempfile()
    gseq.write_meme(m, f, background = c(0.5, 0.5, -0, 0))
    expect_false(any(grepl("-0", readLines(f), fixed = TRUE)))
    gseq.write_jaspar(m, f)
    expect_false(any(grepl("-0", readLines(f), fixed = TRUE)))
})

test_that("invalid matrices are errors", {
    f <- withr::local_tempfile()
    expect_error(gseq.write_meme(matrix(0.2, 2, 5), f), "4 columns")
    expect_error(gseq.write_meme(matrix(c(-1, 1, 0.5, 0.5), 1), f), "negative")
    expect_error(gseq.write_meme(matrix(c(NA, 1, 0, 0), 1), f), "missing")
    expect_error(gseq.write_meme(matrix(numeric(0), 0, 4), f), "no positions")
    expect_error(gseq.write_meme(list(), f), "non-empty")
    expect_error(gseq.write_meme(list(a = "ACGT"), f), "numeric matrix")
    expect_error(gseq.write_meme(matrix(0.25, 1, 4), c("a", "b")), "single file path")
    expect_error(
        gseq.write_meme(matrix(0.25, 1, 4, dimnames = list(NULL, c("A", "C", "G", "U"))), f),
        "4 columns"
    )
})

# === Edge cases ===

test_that("width-1 motifs and zero probabilities round-trip in every format", {
    x <- edge_motifs()
    for (fmt in c("meme", "jaspar", "homer")) {
        f <- withr::local_tempfile()
        get(paste0("gseq.write_", fmt))(x, f)
        y <- get(paste0("gseq.read_", fmt))(f)
        expect_length(y, 3)
        for (i in 1:3) {
            expect_lte(max(abs(strip_attrs(x[[i]]) - strip_attrs(y[[i]]))), 5e-7 + 1e-12)
        }
        expect_equal(nrow(y[[1]]), 1)
        expect_equal(unname(y[[2]][1, ]), c(1, 0, 0, 0))
    }
})

test_that("whitespace in names is collapsed and in IDs replaced", {
    m <- matrix(0.25, nrow = 1, ncol = 4)
    attr(m, "name") <- "  my \t motif  name "
    f <- withr::local_tempfile()
    expect_warning(gseq.write_meme(list("my id" = m), f), "whitespace")
    expect_equal(readLines(f)[10], "MOTIF my_id my motif name")
    y <- gseq.read_meme(f)
    expect_equal(names(y), "my_id")
    expect_equal(attr(y[[1]], "name"), "my motif name")

    suppressWarnings(gseq.write_jaspar(list("my id" = m), f))
    y <- gseq.read_jaspar(f)
    expect_equal(names(y), "my_id")
    expect_equal(attr(y[[1]], "name"), "my motif name")

    gseq.write_homer(list(id = m), f)
    expect_equal(attr(gseq.read_homer(f)[[1]], "name"), "my motif name")
})

test_that("an empty or missing name is not written to MEME, and is the ID in JASPAR", {
    m <- matrix(0.25, nrow = 1, ncol = 4)
    attr(m, "name") <- ""
    f <- withr::local_tempfile()
    gseq.write_meme(list(A1 = m), f)
    expect_equal(readLines(f)[10], "MOTIF A1")
    expect_true(is.na(attr(gseq.read_meme(f)[[1]], "name")))
    attr(m, "name") <- NA_character_
    gseq.write_jaspar(list(A1 = m), f)
    expect_equal(readLines(f)[1], ">A1 A1")
    expect_equal(attr(gseq.read_jaspar(f)[[1]], "name"), "A1")
})

test_that("duplicate IDs are written as they are and told apart on reading", {
    m <- matrix(0.25, nrow = 2, ncol = 4)
    x <- list(D = m, D = m, E = m)
    f <- withr::local_tempfile()
    gseq.write_meme(x, f)
    expect_equal(sum(readLines(f) == "MOTIF D"), 2)
    expect_warning(y <- gseq.read_meme(f), "Duplicate")
    expect_equal(names(y), c("D.1", "D.2", "E"))
    gseq.write_jaspar(x, f)
    expect_warning(y <- gseq.read_jaspar(f), "Duplicate")
    expect_equal(names(y), c("D.1", "D.2", "E"))
})

# === MEME options ===

test_that("gseq.write_meme writes the background and strands it is given", {
    m <- list(M = matrix(0.25, nrow = 1, ncol = 4))
    f <- withr::local_tempfile()
    gseq.write_meme(m, f, background = c(G = 2, C = 2, T = 3, A = 3), strands = "+")
    lines <- readLines(f)
    expect_equal(lines[5], "strands: +")
    expect_equal(lines[8], "A 0.300000 C 0.200000 G 0.200000 T 0.300000")
    y <- gseq.read_meme(f)
    expect_equal(attr(y[[1]], "background"), c(A = 0.3, C = 0.2, G = 0.2, T = 0.3))
    expect_equal(attr(y[[1]], "strand"), "+")

    gseq.write_meme(m, f, background = c(0.1, 0.4, 0.4, 0.1))
    expect_equal(readLines(f)[8], "A 0.100000 C 0.400000 G 0.400000 T 0.100000")
})

test_that("gseq.write_meme uses the background and strands the motifs share", {
    m <- matrix(0.25, nrow = 1, ncol = 4)
    attr(m, "background") <- c(A = 0.3, C = 0.2, G = 0.2, T = 0.3)
    attr(m, "strand") <- "+"
    f <- withr::local_tempfile()
    gseq.write_meme(list(a = m, b = m), f)
    expect_equal(readLines(f)[c(5, 8)], c("strands: +", "A 0.300000 C 0.200000 G 0.200000 T 0.300000"))

    m2 <- m
    attr(m2, "background") <- c(A = 0.25, C = 0.25, G = 0.25, T = 0.25)
    expect_warning(gseq.write_meme(list(a = m, b = m2), f), "different backgrounds")
    expect_equal(readLines(f)[8], "A 0.250000 C 0.250000 G 0.250000 T 0.250000")
})

test_that("gseq.write_meme rejects invalid backgrounds and strands", {
    m <- matrix(0.25, nrow = 1, ncol = 4)
    f <- withr::local_tempfile()
    expect_error(gseq.write_meme(m, f, background = c(0.5, 0.5)), "background")
    expect_error(gseq.write_meme(m, f, background = c(A = 0.5, C = 0.5, G = 0, U = 0)), "background")
    expect_error(gseq.write_meme(m, f, background = c(-1, 1, 0.5, 0.5)), "background")
    expect_error(gseq.write_meme(m, f, background = c(0, 0, 0, 0)), "background")
    expect_error(gseq.write_meme(m, f, strands = "-"), "strands")
})

# === JASPAR options ===

test_that("gseq.write_jaspar scales by nsites: argument, then attribute, then 100", {
    m <- matrix(c(0.5, 0.25, 0.25, 0), ncol = 4)
    m2 <- m
    attr(m2, "nsites") <- 8
    f <- withr::local_tempfile()
    gseq.write_jaspar(list(a = m, b = m2), f)
    lines <- readLines(f)
    expect_equal(lines[2], "A  [ 50 ]")
    expect_equal(lines[8], "A  [ 4 ]")
    gseq.write_jaspar(list(a = m, b = m2), f, nsites = 1)
    expect_equal(readLines(f)[c(2, 3, 5)], c("A  [  0.5 ]", "C  [ 0.25 ]", "T  [    0 ]"))
    expect_equal(unname(sapply(gseq.read_jaspar(f), attr, "nsites")), c(1, 1))
    gseq.write_jaspar(list(a = m, b = m2), f, nsites = c(2, 3))
    expect_equal(unname(sapply(gseq.read_jaspar(f), attr, "nsites")), c(2, 3))
    expect_error(gseq.write_jaspar(m, f, nsites = 0), "positive")
    expect_error(gseq.write_jaspar(m, f, nsites = c(1, 2)), "one per motif")
    expect_error(gseq.write_jaspar(m, f, format = "transfac"))
})

test_that("gseq.write_jaspar round_counts writes integer counts", {
    m <- matrix(c(1 / 6, 1 / 6, 1 / 3, 1 / 3, 0.45, 0.4, 0.1, 0.05), ncol = 4, byrow = TRUE)
    f <- withr::local_tempfile()
    gseq.write_jaspar(list(a = m), f, nsites = 3, round_counts = TRUE)
    # 0.5 rounds to even, as in Python
    expect_equal(readLines(f)[2:5], c("A  [ 0 1 ]", "C  [ 0 1 ]", "G  [ 1 0 ]", "T  [ 1 0 ]"))
    gseq.write_jaspar(list(a = m), f, nsites = 1000, round_counts = TRUE)
    y <- gseq.read_jaspar(f)[[1]]
    expect_lte(max(abs(strip_attrs(y) - matrix(m, ncol = 4, dimnames = list(NULL, acgt)))), 1e-3)
    expect_error(gseq.write_jaspar(m, f, round_counts = NA), "round_counts")
})

# === HOMER options ===

test_that("gseq.write_homer computes the consensus by HOMER's rule", {
    m <- matrix(c(
        0.97, 0.01, 0.01, 0.01,
        0.45, 0.40, 0.10, 0.05,
        0.40, 0.05, 0.30, 0.25,
        0.10, 0.30, 0.30, 0.30,
        0.25, 0.25, 0.25, 0.25,
        0.05, 0.40, 0.15, 0.40
    ), ncol = 4, byrow = TRUE)
    f <- withr::local_tempfile()
    gseq.write_homer(list(id = m), f)
    expect_equal(names(gseq.read_homer(f)), "AMRBNY")
})

test_that("a stored consensus is kept unless its width differs", {
    x <- gseq.read_homer(file.path(fixture_dir, "test_motifs.homer.motif"))
    m <- x[[1]]
    attr(m, "consensus") <- "NNNNNNNN"
    f <- withr::local_tempfile()
    gseq.write_homer(list(a = m), f)
    expect_equal(names(gseq.read_homer(f)), "NNNNNNNN")
    attr(m, "consensus") <- "NNN"
    gseq.write_homer(list(a = m), f)
    expect_equal(names(gseq.read_homer(f)), "DTTACGTN")
})

test_that("gseq.write_homer thresholds: argument, then attribute, then a fraction of the best score", {
    m <- matrix(c(1, 0, 0, 0, 0.5, 0.5, 0, 0), ncol = 4, byrow = TRUE)
    best <- log(1 / 0.25) + log(0.5 / 0.25)
    f <- withr::local_tempfile()
    gseq.write_homer(list(a = m), f)
    expect_equal(attr(gseq.read_homer(f)[[1]], "log_odds_threshold"), round(0.6 * best, 6))
    gseq.write_homer(list(a = m), f, threshold_frac = 1)
    expect_equal(attr(gseq.read_homer(f)[[1]], "log_odds_threshold"), round(best, 6))
    attr(m, "log_odds_threshold") <- 1.5
    gseq.write_homer(list(a = m), f)
    expect_equal(attr(gseq.read_homer(f)[[1]], "log_odds_threshold"), 1.5)
    gseq.write_homer(list(a = m, b = m), f, threshold = c(2, 3))
    expect_equal(unname(sapply(suppressWarnings(gseq.read_homer(f)), attr, "log_odds_threshold")), c(2, 3))
    expect_error(gseq.write_homer(m, f, threshold = NA_real_), "threshold")
    expect_error(gseq.write_homer(m, f, threshold_frac = "a"), "threshold_frac")
})

test_that("gseq.write_homer names motifs without a name by their ID", {
    m <- matrix(0.25, nrow = 1, ncol = 4)
    f <- withr::local_tempfile()
    gseq.write_homer(list(MY_ID = m), f)
    expect_equal(attr(gseq.read_homer(f)[[1]], "name"), "MY_ID")
    expect_equal(length(strsplit(readLines(f)[1], "\t")[[1]]), 3)
})
