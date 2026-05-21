# Helpers for BAM auto-import tests. Skipped when samtools isn't on PATH.

skip_if_no_samtools <- function() {
    if (!nzchar(Sys.which("samtools"))) {
        testthat::skip("samtools not on PATH")
    }
}

# Write a synthetic SAM string and convert it to BAM via samtools view -b.
# Returns the path to the BAM file (under tempfile()).
make_test_bam <- function(sam_text) {
    skip_if_no_samtools()
    sam <- tempfile(fileext = ".sam")
    bam <- tempfile(fileext = ".bam")
    writeLines(sam_text, sam)
    res <- system2("samtools", c("view", "-b", "-o", bam, sam),
        stdout = TRUE, stderr = TRUE
    )
    if (!file.exists(bam) || file.info(bam)$size == 0) {
        stop("samtools view -b failed: ", paste(res, collapse = "\n"))
    }
    bam
}

# A SAM fixture covering chrom "1" from the test DB (which uses bare numeric
# names, not "chr1"), with two reads on different strands plus one unmapped
# read with an unknown chrom (exercises the chrom-mismatch path).
default_sam_text <- function() {
    c(
        "@HD\tVN:1.6",
        "@SQ\tSN:1\tLN:247249719",
        "r1\t0\t1\t100\t30\t10M\t*\t0\t0\tAAAAAAAAAA\t*",
        "r2\t16\t1\t200\t30\t10M\t*\t0\t0\tAAAAAAAAAA\t*",
        "r3\t4\tunknown_chrom\t1\t0\t*\t*\t0\t0\t*\t*"
    )
}
