create_isolated_test_db()

test_that("import and extract from s_7_export.txt", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    intervs <- gscreen("test.fixedbin > 0.1", gintervals(c(1, 2)))
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.rm(tmptrack, force = TRUE)
    gtrack.import_mappedseq(tmptrack, "", "/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/input_files/s_7_export.txt", remove.dups = FALSE)
    r <- gextract(tmptrack, intervs, colnames = "test.tmptrack")
    expect_regression(r, "track.import_mappedseq.s_7_export")
})

test_that("import and extract from sample-small.sam", {
    sam <- "/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/input_files/sample-small.sam"
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gtrack.rm(tmptrack, force = TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.import_mappedseq(tmptrack, "", sam, cols.order = NULL, remove.dups = FALSE)
    r <- gextract(tmptrack, gintervals.all(), colnames = "v")

    # Expected 5' bases, computed independently: POS is 1-based; a reverse read ends at POS - 1 + its CIGAR span.
    reads <- read.table(pipe(paste("grep -v '^@'", sam, "| cut -f2-4,6")),
        sep = "\t", col.names = c("flag", "chrom", "pos", "cigar"), stringsAsFactors = FALSE
    )
    reads <- reads[bitwAnd(reads$flag, 0x904) == 0, ]
    ops <- regmatches(reads$cigar, gregexpr("[0-9]+[MDN=X]", reads$cigar))
    span <- vapply(ops, function(o) sum(as.numeric(sub("[MDN=X]", "", o))), numeric(1))
    reverse <- bitwAnd(reads$flag, 0x10) != 0
    point <- ifelse(reverse, reads$pos - 1 + span - 1, reads$pos - 1)
    chroms <- gintervals.all()
    expected <- as.data.frame(table(chrom = reads$chrom, start = point), stringsAsFactors = FALSE)
    expected <- expected[expected$Freq > 0 & expected$chrom %in% chroms$chrom, ]
    expected$start <- as.numeric(expected$start)
    expected <- expected[expected$start < chroms$end[match(expected$chrom, chroms$chrom)], ]
    expected <- expected[order(expected$chrom, expected$start), ]

    got <- data.frame(chrom = as.character(r$chrom), start = r$start, v = r$v, stringsAsFactors = FALSE)
    got <- got[order(got$chrom, got$start), ]
    expect_equal(nrow(got), nrow(expected))
    expect_equal(got$chrom, expected$chrom)
    expect_equal(got$start, expected$start)
    expect_equal(got$v, as.numeric(expected$Freq))
})

test_that("import with pileup and binsize from s_7_export.txt", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    intervs <- gscreen("test.fixedbin > 0.1", gintervals(c(1, 2)))
    gtrack.rm(tmptrack, force = TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.import_mappedseq(tmptrack, "", "/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/input_files/s_7_export.txt", remove.dups = FALSE, pileup = 180, binsize = 50)
    r <- gextract(tmptrack, intervs, colnames = "test.tmptrack")
    expect_regression(r, "track.import_pileup_binsize")
})
