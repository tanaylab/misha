create_isolated_test_db()

# Test DB canonical chromosome names are "chr1", "chr2", ... (gsetroot prepends
# "chr" to bare names from chrom_sizes.txt). So Ensembl-style bare inputs like
# "1" should be auto-aliased to "chr1".

test_that(".gchroms adds chr prefix when DB has prefixed names", {
    expect_equal(as.character(misha:::.gchroms("1")), "chr1")
    expect_equal(as.character(misha:::.gchroms(c("1", "2"))), c("chr1", "chr2"))
})

test_that(".gchroms passes through canonical names unchanged", {
    expect_equal(as.character(misha:::.gchroms("chr1")), "chr1")
    expect_equal(as.character(misha:::.gchroms(c("chr1", "chr2"))), c("chr1", "chr2"))
})

test_that(".gchroms still errors on truly unknown chromosomes with helpful message", {
    expect_error(misha:::.gchroms("frobnitz"), "Known chromosomes")
})

test_that("gintervals accepts both prefix styles", {
    a <- gintervals(1, 0, 100)
    b <- gintervals("chr1", 0, 100)
    expect_equal(a, b)
    expect_equal(as.character(a$chrom), "chr1")
})

test_that(".gchroms toggles prefix on a mixed input vector", {
    out <- as.character(misha:::.gchroms(c("chr1", "2", "3", "chr20")))
    expect_equal(out, c("chr1", "chr2", "chr3", "chr20"))
})

test_that("start>=end error mentions 0-based half-open convention", {
    expect_error(
        gintervals(1, 100, 100),
        "0-based"
    )
})

test_that("gintervals.import_bed reads basic BED3", {
    f <- tempfile(fileext = ".bed")
    on.exit(unlink(f))
    writeLines(c(
        "track name=\"x\"",
        "chr1\t100\t200",
        "1\t500\t1500" # mixed-prefix input
    ), f)
    out <- gintervals.import_bed(f)
    expect_equal(nrow(out), 2)
    expect_true(all(as.character(out$chrom) == "chr1"))
    expect_equal(sort(out$start), c(100, 500))
})

test_that("gintervals.import_bed parses BED6 strand", {
    f <- tempfile(fileext = ".bed")
    on.exit(unlink(f))
    writeLines(c(
        "chr1\t100\t200\tn1\t0.5\t+",
        "chr1\t300\t400\tn2\t0.7\t-",
        "chr2\t500\t600\tn3\t.\t."
    ), f)
    out <- gintervals.import_bed(f)
    expect_equal(nrow(out), 3)
    expect_true("strand" %in% colnames(out))
    expect_true("name" %in% colnames(out))
    expect_true("score" %in% colnames(out))
    # Sorted output: strand should be numeric ±1/0
    expect_equal(out$strand, c(1, -1, 0))
    expect_equal(out$name, c("n1", "n2", "n3"))
})

test_that("gintervals.import_gff converts 1-based closed to 0-based half-open", {
    f <- tempfile(fileext = ".gff")
    on.exit(unlink(f))
    writeLines(c(
        "##gff-version 3",
        "chr1\ttest\texon\t101\t200\t.\t+\t.\tgene=A",
        "chr1\ttest\texon\t301\t400\t0.5\t-\t.\tgene=B"
    ), f)
    out <- gintervals.import_gff(f)
    expect_equal(nrow(out), 2)
    # GFF 101..200 inclusive 1-based -> 100..200 half-open 0-based
    o <- out[order(out$start), ]
    expect_equal(o$start, c(100, 300))
    expect_equal(o$end, c(200, 400))
    expect_equal(o$strand, c(1, -1))
    expect_equal(o$type, c("exon", "exon"))
    expect_equal(o$attrs, c("gene=A", "gene=B"))
})

test_that("gintervals.import_gff feature filter works", {
    f <- tempfile(fileext = ".gff")
    on.exit(unlink(f))
    writeLines(c(
        "chr1\ttest\tgene\t1\t1000\t.\t+\t.\tg",
        "chr1\ttest\texon\t101\t200\t.\t+\t.\te"
    ), f)
    out <- gintervals.import_gff(f, feature = "exon")
    expect_equal(nrow(out), 1)
    expect_equal(out$type, "exon")
})

test_that("gintervals.import_vcf converts POS to 0-based half-open via REF length", {
    f <- tempfile(fileext = ".vcf")
    on.exit(unlink(f))
    writeLines(c(
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "chr1\t101\trs1\tA\tG\t30\tPASS\tAC=1",
        "chr1\t201\trs2\tACG\tA\t40\tPASS\tAC=2"
    ), f)
    out <- gintervals.import_vcf(f)
    expect_equal(nrow(out), 2)
    o <- out[order(out$start), ]
    # POS=101 1-based, REF=A (len 1) -> 0-based [100, 101)
    # POS=201, REF=ACG (len 3) -> 0-based [200, 203)
    expect_equal(o$start, c(100, 200))
    expect_equal(o$end, c(101, 203))
    expect_equal(o$id, c("rs1", "rs2"))
    expect_equal(o$ref, c("A", "ACG"))
    expect_equal(o$alt, c("G", "A"))
})

test_that("gintervals.import_bed name=FALSE/score=FALSE/strand=FALSE drop those columns", {
    f <- tempfile(fileext = ".bed")
    on.exit(unlink(f))
    writeLines(c("chr1\t100\t200\tn\t1.0\t+"), f)
    out <- gintervals.import_bed(f, name = FALSE, score = FALSE, strand = FALSE)
    expect_equal(colnames(out), c("chrom", "start", "end"))
})
