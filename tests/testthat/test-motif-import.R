fixture_dir <- testthat::test_path("fixtures", "motifs")

# === MEME Format Tests ===

test_that("gseq.read_meme parses multiple motifs correctly", {
    f <- file.path(fixture_dir, "test_motifs.meme")
    motifs <- gseq.read_meme(f)
    expect_length(motifs, 3)
    expect_equal(names(motifs), c("MA0001.1", "MA0002.1", "MA0003.2"))
})

test_that("gseq.read_meme returns correct matrix dimensions", {
    f <- file.path(fixture_dir, "test_motifs.meme")
    motifs <- gseq.read_meme(f)
    expect_equal(nrow(motifs[["MA0001.1"]]), 4)
    expect_equal(nrow(motifs[["MA0002.1"]]), 6)
    expect_equal(nrow(motifs[["MA0003.2"]]), 5)
    expect_equal(ncol(motifs[["MA0001.1"]]), 4)
})

test_that("gseq.read_meme returns correct column names", {
    f <- file.path(fixture_dir, "test_motifs.meme")
    motifs <- gseq.read_meme(f)
    expect_equal(colnames(motifs[[1]]), c("A", "C", "G", "T"))
    expect_equal(colnames(motifs[[2]]), c("A", "C", "G", "T"))
    expect_equal(colnames(motifs[[3]]), c("A", "C", "G", "T"))
})

test_that("gseq.read_meme parses correct matrix values", {
    f <- file.path(fixture_dir, "test_motifs.meme")
    motifs <- gseq.read_meme(f)
    m1 <- motifs[["MA0001.1"]]
    expect_equal(m1[1, ], c(A = 0.250, C = 0.500, G = 0.125, T = 0.125))
    expect_equal(m1[2, ], c(A = 0.875, C = 0.000, G = 0.125, T = 0.000))
    expect_equal(m1[3, ], c(A = 0.000, C = 0.000, G = 1.000, T = 0.000))
    expect_equal(m1[4, ], c(A = 0.100, C = 0.100, G = 0.100, T = 0.700))
})

test_that("gseq.read_meme rows sum to 1.0", {
    f <- file.path(fixture_dir, "test_motifs.meme")
    motifs <- gseq.read_meme(f)
    for (m in motifs) {
        expect_equal(rowSums(m), rep(1.0, nrow(m)), tolerance = 1e-6)
    }
})

test_that("gseq.read_meme parses attributes correctly", {
    f <- file.path(fixture_dir, "test_motifs.meme")
    motifs <- gseq.read_meme(f)

    m1 <- motifs[["MA0001.1"]]
    expect_equal(attr(m1, "name"), "AGL3")
    expect_equal(attr(m1, "alength"), 4L)
    expect_equal(attr(m1, "w"), 4L)
    expect_equal(attr(m1, "nsites"), 20)
    expect_equal(attr(m1, "E"), 1.2e-005)
    expect_equal(attr(m1, "url"), "http://jaspar.genereg.net/matrix/MA0001.1")
    expect_equal(attr(m1, "strand"), "+ -")

    bg <- attr(m1, "background")
    expect_equal(bg, c(A = 0.29, C = 0.21, G = 0.21, T = 0.29))

    m3 <- motifs[["MA0003.2"]]
    expect_equal(attr(m3, "name"), "TFAP2A")
    expect_true(is.na(attr(m3, "E")))
    expect_true(is.na(attr(m3, "url")))
})

test_that("gseq.read_meme errors on nonexistent file", {
    expect_error(gseq.read_meme("/nonexistent/file.meme"), "File not found")
})

test_that("gseq.read_meme errors on empty file", {
    tmp <- tempfile(fileext = ".meme")
    writeLines(character(0), tmp)
    expect_error(gseq.read_meme(tmp), "No motifs found")
    unlink(tmp)
})

test_that("gseq.read_meme errors on log-odds matrix", {
    tmp <- tempfile(fileext = ".meme")
    writeLines(c(
        "MEME version 5",
        "",
        "MOTIF test_motif",
        "letter-probability matrix: alength= 4 w= 2",
        " 0.5  0.2  0.1 -0.3",
        " 0.1  0.1  0.7  0.1"
    ), tmp)
    expect_error(gseq.read_meme(tmp), "Log-odds matrices not supported")
    unlink(tmp)
})

test_that("gseq.read_meme warns and re-normalizes non-unit rows", {
    tmp <- tempfile(fileext = ".meme")
    writeLines(c(
        "MOTIF test_renorm",
        "letter-probability matrix: alength= 4 w= 2",
        " 0.3  0.3  0.3  0.3",
        " 0.25 0.25 0.25 0.25"
    ), tmp)
    expect_warning(motifs <- gseq.read_meme(tmp), "re-normalizing")
    expect_equal(rowSums(motifs[[1]]), c(1.0, 1.0), tolerance = 1e-6)
    unlink(tmp)
})

test_that("gseq.read_meme handles duplicate IDs with warning", {
    tmp <- tempfile(fileext = ".meme")
    writeLines(c(
        "MOTIF DUP1",
        "letter-probability matrix: alength= 4 w= 2",
        " 0.25 0.25 0.25 0.25",
        " 0.25 0.25 0.25 0.25",
        "",
        "MOTIF DUP1",
        "letter-probability matrix: alength= 4 w= 2",
        " 0.50 0.50 0.00 0.00",
        " 0.00 0.00 0.50 0.50"
    ), tmp)
    expect_warning(motifs <- gseq.read_meme(tmp), "Duplicate motif IDs")
    expect_equal(names(motifs), c("DUP1.1", "DUP1.2"))
    unlink(tmp)
})

test_that("gseq.read_meme handles single motif file", {
    tmp <- tempfile(fileext = ".meme")
    writeLines(c(
        "MOTIF SINGLE test_name",
        "letter-probability matrix: alength= 4 w= 3",
        " 0.25 0.25 0.25 0.25",
        " 0.10 0.70 0.10 0.10",
        " 0.70 0.10 0.10 0.10"
    ), tmp)
    motifs <- gseq.read_meme(tmp)
    expect_length(motifs, 1)
    expect_equal(names(motifs), "SINGLE")
    expect_equal(attr(motifs[[1]], "name"), "test_name")
    unlink(tmp)
})


# === JASPAR Header Format Tests ===

test_that("gseq.read_jaspar parses header format correctly", {
    f <- file.path(fixture_dir, "test_motifs.jaspar")
    motifs <- gseq.read_jaspar(f)
    expect_length(motifs, 2)
    expect_equal(names(motifs), c("MA0002.1", "MA0004.1"))
})

test_that("gseq.read_jaspar returns correct column names", {
    f <- file.path(fixture_dir, "test_motifs.jaspar")
    motifs <- gseq.read_jaspar(f)
    expect_equal(colnames(motifs[[1]]), c("A", "C", "G", "T"))
    expect_equal(colnames(motifs[[2]]), c("A", "C", "G", "T"))
})

test_that("gseq.read_jaspar rows sum to 1.0", {
    f <- file.path(fixture_dir, "test_motifs.jaspar")
    motifs <- gseq.read_jaspar(f)
    for (m in motifs) {
        expect_equal(rowSums(m), rep(1.0, nrow(m)), tolerance = 1e-6)
    }
})

test_that("gseq.read_jaspar converts counts to probabilities correctly", {
    f <- file.path(fixture_dir, "test_motifs.jaspar")
    motifs <- gseq.read_jaspar(f)
    m <- motifs[["MA0002.1"]]
    # Position 1: A=10, C=0, G=3, T=2, total=15
    expect_equal(unname(m[1, "A"]), 10 / 15, tolerance = 1e-6)
    expect_equal(unname(m[1, "C"]), 0 / 15, tolerance = 1e-6)
    expect_equal(unname(m[1, "G"]), 3 / 15, tolerance = 1e-6)
    expect_equal(unname(m[1, "T"]), 2 / 15, tolerance = 1e-6)
})

test_that("gseq.read_jaspar sets correct attributes", {
    f <- file.path(fixture_dir, "test_motifs.jaspar")
    motifs <- gseq.read_jaspar(f)
    m <- motifs[["MA0002.1"]]
    expect_equal(attr(m, "name"), "RUNX1")
    expect_equal(attr(m, "w"), 6L)
    expect_equal(attr(m, "nsites"), 15)
    expect_equal(attr(m, "format"), "jaspar")
})

test_that("gseq.read_jaspar errors on nonexistent file", {
    expect_error(gseq.read_jaspar("/nonexistent/file.jaspar"), "File not found")
})


# === JASPAR Simple PFM Format Tests ===

test_that("gseq.read_jaspar parses simple PFM format", {
    f <- file.path(fixture_dir, "test_motifs_simple.pfm")
    motifs <- gseq.read_jaspar(f)
    expect_length(motifs, 1)
    expect_equal(attr(motifs[[1]], "format"), "simple")
})

test_that("gseq.read_jaspar simple PFM rows sum to 1.0", {
    f <- file.path(fixture_dir, "test_motifs_simple.pfm")
    motifs <- gseq.read_jaspar(f)
    expect_equal(rowSums(motifs[[1]]), rep(1.0, nrow(motifs[[1]])), tolerance = 1e-6)
})

test_that("gseq.read_jaspar simple PFM correct values", {
    f <- file.path(fixture_dir, "test_motifs_simple.pfm")
    motifs <- gseq.read_jaspar(f)
    m <- motifs[[1]]
    # Position 1: A=10, C=0, G=3, T=2, total=15
    expect_equal(unname(m[1, "A"]), 10 / 15, tolerance = 1e-6)
    expect_equal(unname(m[1, "C"]), 0 / 15, tolerance = 1e-6)
    # Position 4: A=1, C=0, G=0, T=14, total=15
    expect_equal(unname(m[4, "T"]), 14 / 15, tolerance = 1e-6)
})

test_that("gseq.read_jaspar simple PFM uses filename as ID", {
    f <- file.path(fixture_dir, "test_motifs_simple.pfm")
    motifs <- gseq.read_jaspar(f)
    expect_equal(names(motifs), "test_motifs_simple")
})


# === HOMER Format Tests ===

test_that("gseq.read_homer parses multiple motifs correctly", {
    f <- file.path(fixture_dir, "test_motifs.homer.motif")
    motifs <- gseq.read_homer(f)
    expect_length(motifs, 2)
    expect_equal(names(motifs), c("ATGACTCA", "GATAAG"))
})

test_that("gseq.read_homer returns correct column names", {
    f <- file.path(fixture_dir, "test_motifs.homer.motif")
    motifs <- gseq.read_homer(f)
    expect_equal(colnames(motifs[[1]]), c("A", "C", "G", "T"))
    expect_equal(colnames(motifs[[2]]), c("A", "C", "G", "T"))
})

test_that("gseq.read_homer parses correct matrix values", {
    f <- file.path(fixture_dir, "test_motifs.homer.motif")
    motifs <- gseq.read_homer(f)
    m1 <- motifs[["ATGACTCA"]]
    expect_equal(unname(m1[1, "A"]), 0.353, tolerance = 1e-6)
    expect_equal(unname(m1[1, "C"]), 0.143, tolerance = 1e-6)
    expect_equal(nrow(m1), 8)
    m2 <- motifs[["GATAAG"]]
    expect_equal(nrow(m2), 6)
})

test_that("gseq.read_homer rows sum to 1.0", {
    f <- file.path(fixture_dir, "test_motifs.homer.motif")
    motifs <- gseq.read_homer(f)
    for (m in motifs) {
        expect_equal(rowSums(m), rep(1.0, nrow(m)), tolerance = 1e-6)
    }
})

test_that("gseq.read_homer parses attributes correctly", {
    f <- file.path(fixture_dir, "test_motifs.homer.motif")
    motifs <- gseq.read_homer(f)

    m1 <- motifs[["ATGACTCA"]]
    expect_equal(attr(m1, "consensus"), "ATGACTCA")
    expect_equal(attr(m1, "name"), "AP1(bZIP)/ThP1-cJun")
    expect_equal(attr(m1, "log_odds_threshold"), 8.036341, tolerance = 1e-6)
    expect_equal(attr(m1, "log_p_value"), -4130.668834, tolerance = 1e-4)
    expect_equal(attr(m1, "w"), 8L)
    expect_equal(attr(m1, "source"), "homer")
})

test_that("gseq.read_homer errors on nonexistent file", {
    expect_error(gseq.read_homer("/nonexistent/file.motif"), "File not found")
})

test_that("gseq.read_homer errors on file with no headers", {
    tmp <- tempfile(fileext = ".motif")
    writeLines(c("0.25\t0.25\t0.25\t0.25"), tmp)
    expect_error(gseq.read_homer(tmp), "No motifs found")
    unlink(tmp)
})


# === Integration: matrices work with gseq.pwm() ===

create_isolated_test_db()

test_that("MEME-parsed matrices work with gseq.pwm()", {
    f <- file.path(fixture_dir, "test_motifs.meme")
    motifs <- gseq.read_meme(f)
    m <- motifs[["MA0001.1"]]
    seqs <- c("ACGTACGT", "GGGGGGGG")
    scores <- gseq.pwm(seqs, m, mode = "max")
    expect_length(scores, 2)
    expect_true(all(is.numeric(scores)))
})

test_that("JASPAR-parsed matrices work with gseq.pwm()", {
    f <- file.path(fixture_dir, "test_motifs.jaspar")
    motifs <- gseq.read_jaspar(f)
    m <- motifs[["MA0002.1"]]
    seqs <- c("ACGTACGTACGT", "GGGGGGGGGGGG")
    scores <- gseq.pwm(seqs, m, mode = "max")
    expect_length(scores, 2)
    expect_true(all(is.numeric(scores)))
})

test_that("HOMER-parsed matrices work with gseq.pwm()", {
    f <- file.path(fixture_dir, "test_motifs.homer.motif")
    motifs <- gseq.read_homer(f)
    m <- motifs[["ATGACTCA"]]
    seqs <- c("ATGACTCATC", "CCCCCCCCCC")
    scores <- gseq.pwm(seqs, m, mode = "max")
    expect_length(scores, 2)
    expect_true(all(is.numeric(scores)))
})


# === Edge Cases ===

test_that("gseq.read_meme handles Windows line endings", {
    tmp <- tempfile(fileext = ".meme")
    # Write with \r\n
    con <- file(tmp, open = "wb")
    writeBin(charToRaw(paste0(
        "MOTIF WIN1\r\n",
        "letter-probability matrix: alength= 4 w= 2\r\n",
        " 0.25 0.25 0.25 0.25\r\n",
        " 0.10 0.70 0.10 0.10\r\n"
    )), con)
    close(con)
    motifs <- gseq.read_meme(tmp)
    expect_length(motifs, 1)
    expect_equal(nrow(motifs[[1]]), 2)
    unlink(tmp)
})

test_that("gseq.read_jaspar handles duplicate IDs with warning", {
    tmp <- tempfile(fileext = ".jaspar")
    writeLines(c(
        ">DUP1 MotifA",
        "A [ 10  5 ]",
        "C [  0  5 ]",
        "G [  5  0 ]",
        "T [  5 10 ]",
        ">DUP1 MotifB",
        "A [  5  5 ]",
        "C [  5  5 ]",
        "G [  5  5 ]",
        "T [  5  5 ]"
    ), tmp)
    expect_warning(motifs <- gseq.read_jaspar(tmp), "Duplicate motif IDs")
    expect_equal(names(motifs), c("DUP1.1", "DUP1.2"))
    unlink(tmp)
})

test_that("gseq.read_homer handles duplicate consensus sequences with warning", {
    tmp <- tempfile(fileext = ".motif")
    writeLines(c(
        ">ATCG\tMotifA\t1.0\t-100.0\t0",
        "0.25\t0.25\t0.25\t0.25",
        ">ATCG\tMotifB\t2.0\t-200.0\t0",
        "0.50\t0.50\t0.00\t0.00"
    ), tmp)
    expect_warning(motifs <- gseq.read_homer(tmp), "Duplicate motif IDs")
    expect_equal(names(motifs), c("ATCG.1", "ATCG.2"))
    unlink(tmp)
})

test_that("gseq.read_meme handles file with only header, no motifs", {
    tmp <- tempfile(fileext = ".meme")
    writeLines(c(
        "MEME version 5",
        "",
        "ALPHABET= ACGT",
        "",
        "strands: + -"
    ), tmp)
    expect_error(gseq.read_meme(tmp), "No motifs found")
    unlink(tmp)
})

test_that("gseq.read_jaspar errors on empty file", {
    tmp <- tempfile(fileext = ".jaspar")
    writeLines(character(0), tmp)
    expect_error(gseq.read_jaspar(tmp), "No motifs found")
    unlink(tmp)
})

test_that("gseq.read_homer errors on empty file", {
    tmp <- tempfile(fileext = ".motif")
    writeLines(character(0), tmp)
    expect_error(gseq.read_homer(tmp), "No motifs found")
    unlink(tmp)
})

test_that("gseq.read_meme errors on wrong column count", {
    tmp <- tempfile(fileext = ".meme")
    writeLines(c(
        "MOTIF BAD",
        "letter-probability matrix: alength= 4 w= 2",
        " 0.25 0.25 0.25 0.25",
        " 0.25 0.25 0.50"
    ), tmp)
    expect_error(gseq.read_meme(tmp), "Expected 4 columns")
    unlink(tmp)
})

test_that("gseq.read_meme errors when w does not match data rows", {
    tmp <- tempfile(fileext = ".meme")
    writeLines(c(
        "MOTIF WMISMATCH",
        "letter-probability matrix: alength= 4 w= 3",
        " 0.25 0.25 0.25 0.25",
        " 0.25 0.25 0.25 0.25"
    ), tmp)
    expect_error(gseq.read_meme(tmp), "Expected 3 rows but found 2")
    unlink(tmp)
})
