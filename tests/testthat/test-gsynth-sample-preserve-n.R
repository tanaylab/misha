# Tests for the preserve_n option of gsynth.sample().
#
# Test fixture relies on the misha example db chr1, which contains an N
# stretch starting at 0-based position 167280 and continuing through 217279.
# We train a tiny Markov-2 model on a non-N region, then sample on intervals
# that span the ACGT -> N boundary.

# chr1 (example db) N-region facts:
#   chr1[0, 167280)        : non-N
#   chr1[167280, 217280)   : 50000 bp of 'N'
N_START <- 167280L

# Boundary-crossing sample window: 20bp ACGT then 20bp N
SAMPLE_BOUNDARY <- gintervals(1, N_START - 20L, N_START + 20L)

# Pure-N window
SAMPLE_ALL_N <- gintervals(1, N_START + 1000L, N_START + 1100L)

# Training region (well clear of the N tract)
TRAIN_INTERVALS <- gintervals(1, 0, 100000)

.train_tiny_model <- function() {
    if ("gc_vt" %in% gvtrack.ls()) gvtrack.rm("gc_vt")
    gvtrack.create("gc_vt", NULL, "kmer.frac", kmer = "G")
    gsynth.train(
        list(expr = "gc_vt", breaks = seq(0, 1, length.out = 6)),
        intervals = TRAIN_INTERVALS,
        iterator = 200,
        k = 2L
    )
}

.read_fasta_seq <- function(path) {
    lines <- readLines(path)
    paste(lines[!grepl("^>", lines)], collapse = "")
}

.read_original <- function(intervals) {
    # Read the literal bytes of the misha test db reference for these intervals.
    seq_dir <- file.path(.misha$GROOT, "seq")
    chrom_name <- as.character(intervals$chrom[1])
    seq_path <- file.path(seq_dir, paste0(chrom_name, ".seq"))
    con <- file(seq_path, "rb")
    on.exit(close(con))
    seek(con, intervals$start[1])
    readChar(con, intervals$end[1] - intervals$start[1], useBytes = TRUE)
}

test_that("gsynth.sample preserves N positions by default", {
    gdb.init_examples()
    model <- .train_tiny_model()

    out <- tempfile(fileext = ".fa")
    gsynth.sample(
        model, out,
        output_format = "fasta",
        intervals = SAMPLE_BOUNDARY,
        seed = 60427
    )

    seq <- .read_fasta_seq(out)
    expect_equal(nchar(seq), 40)

    # Original has 20 N's then 20 ACGT
    orig <- .read_original(SAMPLE_BOUNDARY)
    expect_equal(nchar(orig), 40)

    # All original-N positions must be N in the output
    n_pos <- which(strsplit(orig, "")[[1]] == "N")
    out_chars <- strsplit(seq, "")[[1]]
    expect_true(length(n_pos) > 0)
    expect_true(all(out_chars[n_pos] == "N"))

    # Non-N positions must be ACGT (uppercase from the Markov sampler)
    non_n_pos <- setdiff(seq_along(out_chars), n_pos)
    expect_true(all(out_chars[non_n_pos] %in% c("A", "C", "G", "T")))

    unlink(out)
    gvtrack.rm("gc_vt")
})

test_that("gsynth.sample with preserve_n = FALSE samples over N (legacy behavior)", {
    gdb.init_examples()
    model <- .train_tiny_model()

    out <- tempfile(fileext = ".fa")
    gsynth.sample(
        model, out,
        output_format = "fasta",
        intervals = SAMPLE_ALL_N,
        preserve_n = FALSE,
        seed = 60427
    )

    seq <- .read_fasta_seq(out)
    expect_equal(nchar(seq), 100)

    # With preserve_n disabled and a 100-bp all-N input, the output is sampled
    # uniformly (k-mer context contains N -> uniform fallback). With prob
    # ~ 4^-100 of seeing zero ACGT bases, we can safely require at least one.
    expect_true(any(strsplit(seq, "")[[1]] %in% c("A", "C", "G", "T")))

    unlink(out)
    gvtrack.rm("gc_vt")
})

test_that("gsynth.sample preserve_n = TRUE produces all-N output for an all-N region", {
    gdb.init_examples()
    model <- .train_tiny_model()

    out <- tempfile(fileext = ".fa")
    gsynth.sample(
        model, out,
        output_format = "fasta",
        intervals = SAMPLE_ALL_N,
        preserve_n = TRUE,
        seed = 60427
    )

    seq <- .read_fasta_seq(out)
    expect_equal(seq, strrep("N", 100))

    unlink(out)
    gvtrack.rm("gc_vt")
})

test_that("mask_copy interval covering an N region keeps N regardless of preserve_n", {
    gdb.init_examples()
    model <- .train_tiny_model()

    mask_copy_iv <- SAMPLE_BOUNDARY # full boundary span copies original

    for (preserve in c(TRUE, FALSE)) {
        out <- tempfile(fileext = ".fa")
        gsynth.sample(
            model, out,
            output_format = "fasta",
            intervals = SAMPLE_BOUNDARY,
            mask_copy = mask_copy_iv,
            preserve_n = preserve,
            seed = 60427
        )

        seq <- .read_fasta_seq(out)
        # mask_copy copies original bytes verbatim, including soft-masked case.
        expect_equal(seq, .read_original(SAMPLE_BOUNDARY))

        unlink(out)
    }

    gvtrack.rm("gc_vt")
})

test_that("preserve_n works for vector and misha output formats", {
    gdb.init_examples()
    model <- .train_tiny_model()

    # Vector output
    seqs <- gsynth.sample(
        model,
        output_format = "vector",
        intervals = SAMPLE_BOUNDARY,
        preserve_n = TRUE,
        seed = 60427
    )
    expect_length(seqs, 1L)
    expect_equal(nchar(seqs[1]), 40)
    orig <- .read_original(SAMPLE_BOUNDARY)
    n_pos <- which(strsplit(orig, "")[[1]] == "N")
    expect_true(all(strsplit(seqs[1], "")[[1]][n_pos] == "N"))

    # misha .seq binary output
    out_seq <- tempfile(fileext = ".seq")
    gsynth.sample(
        model, out_seq,
        output_format = "misha",
        intervals = SAMPLE_BOUNDARY,
        preserve_n = TRUE,
        seed = 60427
    )
    raw <- readBin(out_seq, "raw", n = file.info(out_seq)$size)
    chars <- rawToChar(raw)
    expect_equal(nchar(chars), 40)
    expect_true(all(strsplit(chars, "")[[1]][n_pos] == "N"))

    unlink(out_seq)
    gvtrack.rm("gc_vt")
})

test_that("preserve_n = FALSE matches legacy reproducibility on N-free intervals", {
    # Sanity guard: when there are no N's in the sampled interval, preserve_n
    # has no effect and the two settings must produce identical output for the
    # same seed. Confirms the new opt-out path is byte-equivalent to legacy.
    gdb.init_examples()
    if ("test_vt" %in% gvtrack.ls()) gvtrack.rm("test_vt")
    gvtrack.create("test_vt", "dense_track", "avg")

    intervals <- gintervals(1, 0, 5000) # chr1[0,5000) is N-free
    model <- gsynth.train(
        list(expr = "test_vt", breaks = seq(0, 1, length.out = 6)),
        intervals = intervals,
        iterator = 200,
        k = 2L
    )

    out_a <- tempfile(fileext = ".fa")
    out_b <- tempfile(fileext = ".fa")
    gsynth.sample(model, out_a,
        output_format = "fasta",
        intervals = intervals, preserve_n = TRUE, seed = 60427
    )
    gsynth.sample(model, out_b,
        output_format = "fasta",
        intervals = intervals, preserve_n = FALSE, seed = 60427
    )
    expect_identical(readLines(out_a), readLines(out_b))

    unlink(c(out_a, out_b))
    gvtrack.rm("test_vt")
})
