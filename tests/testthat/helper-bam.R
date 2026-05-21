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

# Memoized BAM for the default fixture: 4 tests use it; building it once
# avoids 4 samtools forks (~100-300ms saved per test-file run).
default_bam_path <- local({
    cached <- NULL
    function() {
        skip_if_no_samtools()
        if (is.null(cached) || !file.exists(cached)) {
            cached <<- make_test_bam(default_sam_text())
        }
        cached
    }
})

# SAM fixture for chrom "chr1" (misha normalizes test-DB chrom names to
# "chr<N>" internally). Two mapped reads on opposite strands and one
# unmapped record whose RNAME field will be rewritten to "*" by samtools
# when converted to BAM, exercising the chrom-mismatch path.
default_sam_text <- function() {
    c(
        "@HD\tVN:1.6",
        "@SQ\tSN:chr1\tLN:247249719",
        "r1\t0\tchr1\t100\t30\t10M\t*\t0\t0\tAAAAAAAAAA\t*",
        "r2\t16\tchr1\t200\t30\t10M\t*\t0\t0\tAAAAAAAAAA\t*",
        "r3\t4\tunknown_chrom\t1\t0\t*\t*\t0\t0\t*\t*"
    )
}
