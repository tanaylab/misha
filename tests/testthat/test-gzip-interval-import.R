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
    # add = TRUE: local_mocked_bindings() above already registered its own
    # restore-the-binding on.exit() handler on this same test frame. A bare
    # on.exit() call (default add = FALSE) would silently discard that
    # handler instead of appending to it, leaking the .gdata_table_available
    # mock into whichever test_that() block runs next in this file -- which
    # is exactly what happened here until this was caught: it made the
    # "missing gunzip" test below vacuously pass by skipping the fread/grep
    # branch entirely, regardless of whether that branch was actually fixed.
    on.exit(unlink(c(plain_bed, gz_bed, plain_gff, gz_gff, plain_vcf, gz_vcf)), add = TRUE)

    expect_equal(gintervals.import_bed(gz_bed), gintervals.import_bed(plain_bed))
    expect_equal(gintervals.import_gff(gz_gff), gintervals.import_gff(plain_gff))
    expect_equal(gintervals.import_vcf(gz_vcf), gintervals.import_vcf(plain_vcf))
})

# --- missing gunzip binary ---------------------------------------------------
# When gunzip cannot be run at all, .gunzip_to_tempfile() gets status 127 from
# the shell and returns NULL, and .gread_table_filtered() then skips the
# fread/grep branch entirely and reads the still-compressed file with
# utils::read.table(), which decompresses gzip content itself, independent of
# PATH. That routing matters: handing the gzipped bytes to fread instead (as an
# earlier version of this fix did, via a Sys.which("gunzip") guard) drops the
# header-line pre-filter, and fread locks its column count from the first line,
# so a gz BED with a leading "track" header came back as a single (wrong)
# column -- i.e. it silently reintroduced the exact bug this file exists to
# catch. Simulate "gunzip missing" by restricting PATH, for the duration of
# this test only, to a directory containing nothing but a symlink to the
# real grep -- so grep keeps working and gunzip genuinely cannot be found.
test_that("gzipped BED import still works when gunzip is unavailable on PATH", {
    grep_path <- Sys.which("grep")
    skip_if(!nzchar(grep_path), "grep not found on this system")

    restricted_dir <- tempfile("path_no_gunzip_")
    dir.create(restricted_dir)
    file.symlink(grep_path, file.path(restricted_dir, "grep"))
    old_path <- Sys.getenv("PATH")
    Sys.setenv(PATH = restricted_dir)
    on.exit(
        {
            Sys.setenv(PATH = old_path)
            unlink(restricted_dir, recursive = TRUE)
        },
        add = TRUE
    )

    # Sanity-check the simulated environment before trusting the assertions
    # below: grep must still resolve, gunzip must not.
    expect_true(nzchar(Sys.which("grep")))
    expect_false(nzchar(Sys.which("gunzip")))

    plain <- write_plain(bed_lines(), ".bed")
    gz <- write_gz(bed_lines(), ".bed")
    on.exit(unlink(c(plain, gz)), add = TRUE)

    out_plain <- gintervals.import_bed(plain)
    out_gz <- gintervals.import_bed(gz)
    expect_equal(nrow(out_plain), 3)
    expect_equal(out_gz, out_plain)
})

# --- truncated / corrupt gz --------------------------------------------------
# A shell pipeline exits with the status of its *last* command, so
# "gunzip -q -c f | grep -vE pat" reported grep's status, not gunzip's: on a
# truncated archive gunzip died mid-stream, grep still exited 0 on the bytes it
# had already received, fread parsed them, and the import silently returned a
# partial interval set. Decompressing to a temp file first (as gtrack.import()
# has always done) makes the failure visible.
test_that("a truncated gzipped BED errors instead of importing part of the file", {
    skip_if_not_installed("data.table")
    skip_if(!nzchar(Sys.which("gunzip")), "gunzip not found on this system")
    skip_if(!nzchar(Sys.which("grep")), "grep not found on this system")

    # Needs to be big enough that gunzip emits a useful chunk before hitting
    # the truncation - that partial output is what used to be imported.
    n <- 3000
    lines <- c(
        "track name=\"x\"",
        sprintf(
            "chr1\t%d\t%d\tn%d\t0.5\t+",
            seq(0, by = 200, length.out = n), seq(100, by = 200, length.out = n), seq_len(n)
        )
    )
    full <- write_gz(lines, ".bed")
    truncated <- tempfile(fileext = ".bed.gz")
    on.exit(unlink(c(full, truncated)), add = TRUE)
    bytes <- readBin(full, "raw", n = file.info(full)$size)
    writeBin(bytes[seq_len(floor(length(bytes) * 0.6))], truncated)

    # The intact archive imports every interval ...
    expect_equal(nrow(gintervals.import_bed(full)), n)
    # ... and the truncated one is refused rather than silently truncated.
    expect_error(gintervals.import_bed(truncated), "decompress")
})
