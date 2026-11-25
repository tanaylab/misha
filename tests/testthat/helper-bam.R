# Helpers for the BAM import tests. misha reads BAM itself, so none of them needs samtools.
#
# The fixtures in fixtures/bam were made with samtools (see fixtures/bam/make-fixtures.sh). Each
# <name>.view.sam.gz there is `samtools view -h <name>.bam`: the SAM text misha imported when it
# read BAM through samtools, so importing it gives the result a BAM import must reproduce.

bam_fixture <- function(name) {
    testthat::test_path("fixtures", "bam", name)
}

# A database the fixtures fit: chr1, chr2 and chrX of 500, 300 and 200 kb. Roots the session there
# until the calling frame ends. Its tracks/test directory takes the tracks of random_track_name("test").
local_bam_db <- function(env = parent.frame()) {
    restore_groot_on_exit(env)
    db <- tempfile("bam_db_")
    create_test_db(db, data.frame(chrom = c("chr1", "chr2", "chrX"), size = c(500000, 300000, 200000)))
    dir.create(file.path(db, "tracks", "test"))
    withr::defer(unlink(db, recursive = TRUE), envir = env)
    suppressMessages(gsetroot(db))
    invisible(db)
}

# SAM fixture for chrom "chr1": two mapped reads on opposite strands and one unmapped record on a
# chromosome the database does not have. default.bam holds the same records (samtools sets the
# unmapped record's RNAME to "*").
default_sam_text <- function() {
    c(
        "@HD\tVN:1.6",
        "@SQ\tSN:chr1\tLN:500000",
        "r1\t0\tchr1\t100\t30\t10M\t*\t0\t0\tAAAAAAAAAA\t*",
        "r2\t16\tchr1\t200\t30\t10M\t*\t0\t0\tAAAAAAAAAA\t*",
        "r3\t4\tunknown_chrom\t1\t0\t*\t*\t0\t0\t*\t*"
    )
}
