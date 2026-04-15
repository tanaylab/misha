# Local test helper: read a FASTA file into a named list of sequences
read_fasta <- function(fasta_path) {
    lines <- readLines(fasta_path)
    header_idx <- which(startsWith(lines, ">"))
    chrom_names <- sub("^>\\s*", "", sub("\\s.*", "", lines[header_idx]))
    result <- vector("list", length(chrom_names))
    names(result) <- chrom_names
    for (i in seq_along(header_idx)) {
        start <- header_idx[i] + 1L
        end <- if (i < length(header_idx)) header_idx[i + 1L] - 1L else length(lines)
        result[[i]] <- if (start > end) "" else paste0(lines[start:end], collapse = "")
    }
    result
}

test_that("ggenome.implant replaces intervals with literal sequences", {
    local_db_state()

    # Create a simple reference FASTA
    ref_fasta <- tempfile(fileext = ".fa")
    out_fasta <- tempfile(fileext = ".fa")
    withr::defer({
        unlink(ref_fasta)
        unlink(out_fasta)
        unlink(paste0(out_fasta, ".fai"))
    })

    cat(">chrA\nACGTACGTAC\n>chrB\nTTTTAAAA\n", file = ref_fasta)

    intervals <- data.frame(
        chrom = "chrA",
        start = 2,
        end = 6
    )
    donor <- "NNNN"

    result <- ggenome.implant(
        intervals, donor, out_fasta,
        genome_fasta = ref_fasta,
        create_trackdb = FALSE,
        overwrite = TRUE
    )

    expect_equal(result, out_fasta)

    # chrA: positions 0-9 = A C G T A C G T A C
    # Replace [2,6) = positions 2,3,4,5 with NNNN: A C N N N N G T A C
    seqs <- read_fasta(out_fasta)
    expect_equal(seqs[["chrA"]], "ACNNNNGTAC")
    expect_equal(seqs[["chrB"]], "TTTTAAAA")
})

test_that("ggenome.implant handles multiple perturbations on same chromosome", {
    local_db_state()

    ref_fasta <- tempfile(fileext = ".fa")
    out_fasta <- tempfile(fileext = ".fa")
    withr::defer({
        unlink(ref_fasta)
        unlink(out_fasta)
        unlink(paste0(out_fasta, ".fai"))
    })

    # 20 bases: AAAAABBBBBCCCCCDDDD
    cat(">chr1\nAAAAACCCCCGGGGGTTTTT\n", file = ref_fasta)

    intervals <- data.frame(
        chrom = c("chr1", "chr1"),
        start = c(0, 10),
        end = c(5, 15)
    )
    donors <- c("XXXXX", "YYYYY")

    ggenome.implant(
        intervals, donors, out_fasta,
        genome_fasta = ref_fasta,
        create_trackdb = FALSE,
        overwrite = TRUE
    )

    seqs <- read_fasta(out_fasta)
    # Original: AAAAACCCCCGGGGGTTTT
    # Replace [0,5) with XXXXX and [10,15) with YYYYY
    expect_equal(seqs[["chr1"]], "XXXXXCCCCCYYYYYTTTTT")
})

test_that("ggenome.implant creates .fai index", {
    local_db_state()

    ref_fasta <- tempfile(fileext = ".fa")
    out_fasta <- tempfile(fileext = ".fa")
    withr::defer({
        unlink(ref_fasta)
        unlink(out_fasta)
        unlink(paste0(out_fasta, ".fai"))
    })

    cat(">chrA\nACGTACGT\n>chrB\nTTTT\n", file = ref_fasta)

    intervals <- data.frame(chrom = "chrA", start = 0, end = 4)
    donor <- "NNNN"

    ggenome.implant(
        intervals, donor, out_fasta,
        genome_fasta = ref_fasta,
        create_trackdb = FALSE,
        overwrite = TRUE
    )

    fai_path <- paste0(out_fasta, ".fai")
    expect_true(file.exists(fai_path))

    fai <- read.delim(fai_path, header = FALSE, sep = "\t")
    expect_equal(nrow(fai), 2)
    expect_equal(fai$V1, c("chrA", "chrB"))
    # chrA has 8 bases, chrB has 4 bases
    expect_equal(fai$V2, c(8, 4))
})

test_that("ggenome.implant creates trackdb when requested", {
    local_db_state()

    ref_fasta <- tempfile(fileext = ".fa")
    out_fasta <- tempfile(fileext = ".fa")
    trackdb_dir <- tempfile()
    withr::defer({
        unlink(ref_fasta)
        unlink(out_fasta)
        unlink(paste0(out_fasta, ".fai"))
        unlink(trackdb_dir, recursive = TRUE)
    })

    cat(">chrA\nACGTACGT\n", file = ref_fasta)

    intervals <- data.frame(chrom = "chrA", start = 0, end = 4)
    donor <- "TTTT"

    ggenome.implant(
        intervals, donor, out_fasta,
        genome_fasta = ref_fasta,
        create_trackdb = TRUE,
        trackdb_path = trackdb_dir,
        overwrite = TRUE
    )

    expect_true(dir.exists(trackdb_dir))
    expect_true(file.exists(file.path(trackdb_dir, "chrom_sizes.txt")))
})

test_that("ggenome.implant errors on donor/interval length mismatch", {
    local_db_state()

    ref_fasta <- tempfile(fileext = ".fa")
    out_fasta <- tempfile(fileext = ".fa")
    withr::defer({
        unlink(ref_fasta)
        unlink(out_fasta)
    })

    cat(">chrA\nACGTACGT\n", file = ref_fasta)

    intervals <- data.frame(chrom = "chrA", start = 0, end = 4)
    donor <- "NN" # 2 bases, but interval is 4 bases wide

    expect_error(
        ggenome.implant(
            intervals, donor, out_fasta,
            genome_fasta = ref_fasta,
            create_trackdb = FALSE
        ),
        "does not match interval width"
    )
})

test_that("ggenome.implant errors on wrong number of donor sequences", {
    local_db_state()

    ref_fasta <- tempfile(fileext = ".fa")
    out_fasta <- tempfile(fileext = ".fa")
    withr::defer({
        unlink(ref_fasta)
        unlink(out_fasta)
    })

    cat(">chrA\nACGTACGT\n", file = ref_fasta)

    intervals <- data.frame(
        chrom = c("chrA", "chrA"),
        start = c(0, 4),
        end = c(2, 6)
    )
    donor <- "NN" # only 1 donor for 2 intervals

    expect_error(
        ggenome.implant(
            intervals, donor, out_fasta,
            genome_fasta = ref_fasta,
            create_trackdb = FALSE
        ),
        "must equal number of rows"
    )
})

test_that("ggenome.implant errors on missing chromosome", {
    local_db_state()

    ref_fasta <- tempfile(fileext = ".fa")
    out_fasta <- tempfile(fileext = ".fa")
    withr::defer({
        unlink(ref_fasta)
        unlink(out_fasta)
    })

    cat(">chrA\nACGTACGT\n", file = ref_fasta)

    intervals <- data.frame(chrom = "chrZ", start = 0, end = 4)
    donor <- "NNNN"

    expect_error(
        ggenome.implant(
            intervals, donor, out_fasta,
            genome_fasta = ref_fasta,
            create_trackdb = FALSE
        ),
        "not found in reference"
    )
})

test_that("ggenome.implant errors on out-of-bounds intervals", {
    local_db_state()

    ref_fasta <- tempfile(fileext = ".fa")
    out_fasta <- tempfile(fileext = ".fa")
    withr::defer({
        unlink(ref_fasta)
        unlink(out_fasta)
    })

    cat(">chrA\nACGT\n", file = ref_fasta)

    intervals <- data.frame(chrom = "chrA", start = 2, end = 10)
    donor <- "NNNNNNNN"

    expect_error(
        ggenome.implant(
            intervals, donor, out_fasta,
            genome_fasta = ref_fasta,
            create_trackdb = FALSE
        ),
        "out of bounds"
    )
})

test_that("ggenome.implant preserves chromosome order", {
    local_db_state()

    ref_fasta <- tempfile(fileext = ".fa")
    out_fasta <- tempfile(fileext = ".fa")
    withr::defer({
        unlink(ref_fasta)
        unlink(out_fasta)
        unlink(paste0(out_fasta, ".fai"))
    })

    # Intentional order: chrB before chrA
    cat(">chrB\nTTTT\n>chrA\nAAAA\n>chrC\nGGGG\n", file = ref_fasta)

    intervals <- data.frame(chrom = "chrA", start = 0, end = 2)
    donor <- "CC"

    ggenome.implant(
        intervals, donor, out_fasta,
        genome_fasta = ref_fasta,
        create_trackdb = FALSE,
        overwrite = TRUE
    )

    seqs <- read_fasta(out_fasta)
    expect_equal(names(seqs), c("chrB", "chrA", "chrC"))
    expect_equal(seqs[["chrB"]], "TTTT")
    expect_equal(seqs[["chrA"]], "CCAA")
    expect_equal(seqs[["chrC"]], "GGGG")
})

test_that("ggenome.implant respects line_width", {
    local_db_state()

    ref_fasta <- tempfile(fileext = ".fa")
    out_fasta <- tempfile(fileext = ".fa")
    withr::defer({
        unlink(ref_fasta)
        unlink(out_fasta)
        unlink(paste0(out_fasta, ".fai"))
    })

    cat(">chrA\nACGTACGTACGT\n", file = ref_fasta)

    intervals <- data.frame(chrom = "chrA", start = 0, end = 4)
    donor <- "NNNN"

    ggenome.implant(
        intervals, donor, out_fasta,
        genome_fasta = ref_fasta,
        create_trackdb = FALSE,
        line_width = 5L,
        overwrite = TRUE
    )

    lines <- readLines(out_fasta)
    expect_equal(lines, c(">chrA", "NNNNA", "CGTAC", "GT"))
})

test_that("ggenome.implant with donor from misha database", {
    local_db_state()

    # Create a reference FASTA and a donor FASTA
    ref_fasta <- tempfile(fileext = ".fa")
    donor_fasta <- tempfile(fileext = ".fa")
    donor_db <- tempfile()
    out_fasta <- tempfile(fileext = ".fa")
    withr::defer({
        unlink(ref_fasta)
        unlink(donor_fasta)
        unlink(donor_db, recursive = TRUE)
        unlink(out_fasta)
        unlink(paste0(out_fasta, ".fai"))
    })

    cat(">chrA\nAAAAAAAA\n>chrB\nCCCCCCCC\n", file = ref_fasta)
    cat(">chrA\nGGGGGGGG\n>chrB\nTTTTTTTT\n", file = donor_fasta)

    # Create a misha database from the donor FASTA
    suppressMessages(gdb.create(groot = donor_db, fasta = donor_fasta, verbose = FALSE))

    intervals <- data.frame(
        chrom = c("chrA", "chrB"),
        start = c(2, 0),
        end = c(6, 4)
    )

    ggenome.implant(
        intervals,
        donor = donor_db,
        output = out_fasta,
        genome_fasta = ref_fasta,
        create_trackdb = FALSE,
        overwrite = TRUE
    )

    seqs <- read_fasta(out_fasta)
    # chrA: AA + GGGG + AA = AAGGGGAA
    expect_equal(seqs[["chrA"]], "AAGGGGAA")
    # chrB: TTTT + CCCC = TTTTCCCC
    expect_equal(seqs[["chrB"]], "TTTTCCCC")
})

test_that("ggenome.transplant works as sugar for implant", {
    local_db_state()

    ref_fasta <- tempfile(fileext = ".fa")
    donor_fasta <- tempfile(fileext = ".fa")
    donor_db <- tempfile()
    out_fasta <- tempfile(fileext = ".fa")
    withr::defer({
        unlink(ref_fasta)
        unlink(donor_fasta)
        unlink(donor_db, recursive = TRUE)
        unlink(out_fasta)
        unlink(paste0(out_fasta, ".fai"))
    })

    cat(">chrA\nAAAAAAAA\n", file = ref_fasta)
    cat(">chrA\nGGGGGGGG\n", file = donor_fasta)

    suppressMessages(gdb.create(groot = donor_db, fasta = donor_fasta, verbose = FALSE))

    intervals <- data.frame(chrom = "chrA", start = 0, end = 4)

    ggenome.transplant(
        intervals,
        source_genome = donor_db,
        target_genome = ref_fasta,
        output = out_fasta,
        create_trackdb = FALSE,
        overwrite = TRUE
    )

    seqs <- read_fasta(out_fasta)
    expect_equal(seqs[["chrA"]], "GGGGAAAA")
})

test_that("ggenome.implant overwrite guard works", {
    local_db_state()

    ref_fasta <- tempfile(fileext = ".fa")
    out_fasta <- tempfile(fileext = ".fa")
    withr::defer({
        unlink(ref_fasta)
        unlink(out_fasta)
    })

    cat(">chrA\nACGT\n", file = ref_fasta)
    writeLines("placeholder", out_fasta)

    intervals <- data.frame(chrom = "chrA", start = 0, end = 2)
    donor <- "NN"

    expect_error(
        ggenome.implant(
            intervals, donor, out_fasta,
            genome_fasta = ref_fasta,
            create_trackdb = FALSE
        ),
        "already exists"
    )

    # With overwrite = TRUE it should work
    expect_no_error(
        ggenome.implant(
            intervals, donor, out_fasta,
            genome_fasta = ref_fasta,
            create_trackdb = FALSE,
            overwrite = TRUE
        )
    )
})
