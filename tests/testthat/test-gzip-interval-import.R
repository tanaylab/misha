create_isolated_test_db()

# Regression tests for gzipped BED/GFF/VCF import.
#
# .gread_table_filtered (R/intervals-import.R) piped every file straight
# into `grep -vE <header_pat> <file>` before handing the result to
# data.table::fread(cmd = ). On a gzipped file grep sees binary content and
# prints a single literal line ("Binary file p.bed.gz matches"), so fread
# returns a 1x1 data.frame. Since nrow > 0, the code never falls through to
# the utils::read.table() fallback -- which would have worked, because R
# file connections auto-decompress gzip content regardless of extension.
# All three importers document gz support; none of them worked on gz input.

bed_lines <- function() {
    c(
        "track name=\"x\" description=\"y\"",
        "chr1\t100\t200\tn1\t0.5\t+",
        "chr1\t300\t400\tn2\t0.7\t-",
        "chr2\t500\t600\tn3\t.\t."
    )
}

gff_lines <- function() {
    c(
        "##gff-version 3",
        "chr1\ttest\texon\t101\t200\t.\t+\t.\tgene=A",
        "chr1\ttest\texon\t301\t400\t0.5\t-\t.\tgene=B"
    )
}

vcf_lines <- function() {
    c(
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "chr1\t101\trs1\tA\tG\t30\tPASS\tAC=1",
        "chr1\t201\trs2\tACG\tA\t40\tPASS\tAC=2"
    )
}

write_plain <- function(lines, ext) {
    f <- tempfile(fileext = ext)
    writeLines(lines, f)
    f
}

write_gz <- function(lines, ext) {
    f <- tempfile(fileext = paste0(ext, ".gz"))
    con <- gzfile(f, "wt")
    writeLines(lines, con)
    close(con)
    f
}

test_that("gintervals.import_bed: gzipped BED matches uncompressed BED", {
    plain <- write_plain(bed_lines(), ".bed")
    gz <- write_gz(bed_lines(), ".bed")
    on.exit(unlink(c(plain, gz)))

    out_plain <- gintervals.import_bed(plain)
    out_gz <- gintervals.import_bed(gz)
    expect_equal(nrow(out_plain), 3)
    expect_equal(out_gz, out_plain)
})

test_that("gintervals.import_gff: gzipped GFF matches uncompressed GFF", {
    plain <- write_plain(gff_lines(), ".gff")
    gz <- write_gz(gff_lines(), ".gff")
    on.exit(unlink(c(plain, gz)))

    out_plain <- gintervals.import_gff(plain)
    out_gz <- gintervals.import_gff(gz)
    expect_equal(nrow(out_plain), 2)
    expect_equal(out_gz, out_plain)
})

test_that("gintervals.import_vcf: gzipped VCF matches uncompressed VCF", {
    plain <- write_plain(vcf_lines(), ".vcf")
    gz <- write_gz(vcf_lines(), ".vcf")
    on.exit(unlink(c(plain, gz)))

    out_plain <- gintervals.import_vcf(plain)
    out_gz <- gintervals.import_vcf(gz)
    expect_equal(nrow(out_plain), 2)
    expect_equal(out_gz, out_plain)
})

# --- data.table-unavailable fallback ---------------------------------------
# .gread_table_filtered's final fallback (utils::read.table(), reached
# whenever data.table isn't installed) already auto-decompresses gz content
# via R's connection layer and needed no change. Force that branch by
# mocking the internal .gdata_table_available() indirection to FALSE, so
# this test exercises and locks in that behaviour regardless of whether
# data.table happens to be installed on the machine running the suite.
test_that("gzipped import works via the data.table-unavailable fallback", {
    testthat::skip_if_not_installed("testthat", "3.1.7")
    testthat::local_mocked_bindings(
        .gdata_table_available = function() FALSE,
        .package = "misha"
    )

    plain_bed <- write_plain(bed_lines(), ".bed")
    gz_bed <- write_gz(bed_lines(), ".bed")
    plain_gff <- write_plain(gff_lines(), ".gff")
    gz_gff <- write_gz(gff_lines(), ".gff")
    plain_vcf <- write_plain(vcf_lines(), ".vcf")
    gz_vcf <- write_gz(vcf_lines(), ".vcf")
    on.exit(unlink(c(plain_bed, gz_bed, plain_gff, gz_gff, plain_vcf, gz_vcf)))

    expect_equal(gintervals.import_bed(gz_bed), gintervals.import_bed(plain_bed))
    expect_equal(gintervals.import_gff(gz_gff), gintervals.import_gff(plain_gff))
    expect_equal(gintervals.import_vcf(gz_vcf), gintervals.import_vcf(plain_vcf))
})
