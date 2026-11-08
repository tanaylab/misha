create_isolated_test_db()

test_that("gseq.extract with gscreen results on test.fixedbin", {
    intervs <- gscreen("test.fixedbin > 0.6", gintervals(c(1, 2, 3)))
    expect_regression(gseq.extract(intervs), "gseq_extract_gscreen_fixedbin")
})

test_that("gseq.extract with modified strand in gscreen results on test.fixedbin", {
    intervs <- gscreen("test.fixedbin > 0.6", gintervals(c(1, 2, 3)))
    intervs$strand <- -1
    expect_regression(gseq.extract(intervs), "gseq_extract_gscreen_fixedbin_modified_strand")
})

test_that("gseq.extract with gintervals.2d", {
    expect_error(gseq.extract(gintervals.2d(1, 10, 100, 2, 20, 300)))
})

# pad = TRUE: the reference reads the part of each interval inside its chromosome with
# gseq.extract and adds the N runs with strrep(); strand -1 rows are the reverse complement of
# the whole padded forward sequence.
revcomp <- function(x) {
    vapply(x, function(s) {
        paste(rev(strsplit(chartr("ACGTacgt", "TGCAtgca", s), "")[[1]]), collapse = "")
    }, character(1), USE.NAMES = FALSE)
}

pad_ref <- function(intervals) {
    sizes <- gintervals.all()
    lens <- sizes$end[match(as.character(intervals$chrom), as.character(sizes$chrom))]
    cs <- pmax(intervals$start, 0)
    ce <- pmin(intervals$end, lens)
    inside <- cs < ce
    core <- rep("", nrow(intervals))
    if (any(inside)) {
        core[inside] <- withr::with_options(
            list(gmultitasking = FALSE),
            gseq.extract(data.frame(chrom = intervals$chrom[inside], start = cs[inside], end = ce[inside]))
        )
    }
    left <- ifelse(inside, cs - intervals$start, intervals$end - intervals$start)
    right <- ifelse(inside, intervals$end - ce, 0)
    res <- paste0(strrep("N", left), core, strrep("N", right))
    minus <- if (is.null(intervals$strand)) rep(FALSE, nrow(intervals)) else intervals$strand %in% c(-1, "-")
    res[minus] <- revcomp(res[minus])
    res
}

test_that("gseq.extract pad = TRUE fills the part past either chromosome end with N", {
    len1 <- gintervals.all()$end[1]
    chrom1 <- as.character(gintervals.all()$chrom[1])
    chrom2 <- as.character(gintervals.all()$chrom[2])
    len2 <- gintervals.all()$end[2]
    intervs <- data.frame(
        chrom = c(chrom1, chrom1, chrom2, chrom1, chrom1, chrom2, chrom1),
        start = c(-50, len1 - 30, -20, -300, len1 + 10, len2 - 1, 1000),
        end = c(70, len1 + 40, 120, -100, len1 + 60, len2 + 1, 1100)
    )
    seqs <- gseq.extract(intervs, pad = TRUE)
    expect_equal(nchar(seqs), intervs$end - intervs$start)
    expect_identical(seqs, pad_ref(intervs))
    # left end, right end, wholly outside on either side
    expect_identical(substr(seqs[1], 1, 50), strrep("N", 50))
    expect_identical(substr(seqs[2], 31, 70), strrep("N", 40))
    expect_identical(seqs[4], strrep("N", 200))
    expect_identical(seqs[5], strrep("N", 50))
    # rows within their chromosome are read as without pad
    expect_identical(seqs[7], gseq.extract(intervs[7, ]))
    expect_identical(substr(seqs[1], 51, 120), gseq.extract(gintervals(chrom1, 0, 70)))
})

test_that("gseq.extract pad = TRUE covers a window longer than the whole chromosome", {
    sizes <- gintervals.all()
    i <- which.min(sizes$end)
    chrom <- as.character(sizes$chrom[i])
    len <- sizes$end[i]
    seqs <- gseq.extract(data.frame(chrom = chrom, start = -7, end = len + 9), pad = TRUE)
    expect_identical(seqs, paste0(strrep("N", 7), gseq.extract(gintervals(chrom, 0, len)), strrep("N", 9)))
})

test_that("gseq.extract pad = TRUE reverse complements strand -1 rows as a whole", {
    len1 <- gintervals.all()$end[1]
    chrom1 <- as.character(gintervals.all()$chrom[1])
    intervs <- data.frame(
        chrom = chrom1,
        start = c(-50, len1 - 30, -300, 1000, -5),
        end = c(70, len1 + 40, -100, 1100, 20),
        strand = c(-1, -1, -1, -1, 1)
    )
    seqs <- gseq.extract(intervs, pad = TRUE)
    expect_identical(seqs, pad_ref(intervs))
    # the N run of the chromosome start comes last on strand -1
    expect_identical(substr(seqs[1], 71, 120), strrep("N", 50))
    expect_identical(substr(seqs[2], 1, 40), strrep("N", 40))
    expect_identical(seqs[1], revcomp(gseq.extract(intervs[1, 1:3], pad = TRUE)))
    expect_identical(seqs[4], gseq.extract(intervs[4, ]))

    # strand as "+"/"-" text or factor, and integer coordinates, give the same
    intervs_chr <- intervs
    intervs_chr$strand <- ifelse(intervs$strand == -1, "-", "+")
    expect_identical(gseq.extract(intervs_chr, pad = TRUE), seqs)
    intervs_chr$strand <- factor(intervs_chr$strand)
    expect_identical(gseq.extract(intervs_chr, pad = TRUE), seqs)
    intervs_int <- intervs
    intervs_int$start <- as.integer(intervs$start)
    intervs_int$end <- as.integer(intervs$end)
    expect_identical(gseq.extract(intervs_int, pad = TRUE), seqs)
})

test_that("gseq.extract pad = TRUE keeps the input row order and extra columns do not matter", {
    set.seed(11)
    sizes <- gintervals.all()
    n <- 300
    i <- sample(nrow(sizes), n, replace = TRUE)
    centers <- floor(runif(n) * sizes$end[i])
    # a third of the windows sit on a chromosome end
    edge <- sample(c(TRUE, FALSE), n, replace = TRUE, prob = c(1, 2))
    centers[edge] <- ifelse(runif(sum(edge)) < 0.5, 0, sizes$end[i][edge])
    w <- 500
    intervs <- data.frame(
        chrom = as.character(sizes$chrom[i]), start = centers - w / 2, end = centers + w / 2,
        strand = sample(c(-1, 1), n, replace = TRUE), id = seq_len(n)
    )
    seqs <- gseq.extract(intervs, pad = TRUE)
    expect_true(all(nchar(seqs) == w))
    expect_identical(seqs, pad_ref(intervs))
})

test_that("gseq.extract pad = TRUE without a row past an end is the same as pad = FALSE", {
    intervs <- gscreen("test.fixedbin > 0.6", gintervals(c(1, 2, 3)))
    intervs$strand <- rep(c(1, -1), length.out = nrow(intervs))
    expect_identical(gseq.extract(intervs, pad = TRUE), gseq.extract(intervs))
    expect_null(gseq.extract(gintervals(1, 0, 10)[0, ], pad = TRUE))
})

test_that("gseq.extract pad = TRUE pads a named intervals set", {
    chrom1 <- as.character(gintervals.all()$chrom[1])
    intervs <- data.frame(chrom = chrom1, start = c(-5, 100), end = c(10, 120))
    withr::defer(if (gintervals.exists("test.seq_pad_set")) gintervals.rm("test.seq_pad_set", force = TRUE))
    gintervals.save("test.seq_pad_set", intervs)
    expect_identical(gseq.extract("test.seq_pad_set", pad = TRUE), pad_ref(intervs))
})

test_that("gseq.extract without pad still rejects intervals past a chromosome end", {
    len1 <- gintervals.all()$end[1]
    chrom1 <- as.character(gintervals.all()$chrom[1])
    expect_error(
        gseq.extract(data.frame(chrom = chrom1, start = -5, end = 10)),
        "start coordinate must be greater or equal than zero"
    )
    expect_error(
        gseq.extract(data.frame(chrom = chrom1, start = len1 - 5, end = len1 + 5)),
        "end coordinate exceeds chromosome boundaries"
    )
    expect_error(
        gseq.extract(data.frame(chrom = chrom1, start = len1 - 5, end = len1 + 5), pad = FALSE),
        "end coordinate exceeds chromosome boundaries"
    )
})

test_that("gseq.extract pad = TRUE still validates the intervals", {
    chrom1 <- as.character(gintervals.all()$chrom[1])
    expect_error(gseq.extract(data.frame(chrom = chrom1, start = -3, end = -5), pad = TRUE), "lesser than end")
    expect_error(gseq.extract(data.frame(chrom = chrom1, start = 10, end = 5), pad = TRUE), "lesser than end")
    expect_error(gseq.extract(data.frame(chrom = "chrNoSuch", start = -3, end = 5), pad = TRUE), "does not exist")
    expect_error(
        gseq.extract(data.frame(chrom = chrom1, start = -10, end = -5, strand = 5), pad = TRUE),
        "Invalid strand value"
    )
    expect_error(gseq.extract(gintervals.2d(1, 10, 100, 2, 20, 300), pad = TRUE))
    expect_error(gseq.extract(data.frame(chrom = chrom1, start = 0, end = 5), pad = NA), "pad must be TRUE or FALSE")
    expect_error(gseq.extract(data.frame(chrom = chrom1, start = 0, end = 5), pad = "yes"), "pad must be TRUE or FALSE")
})

test_that("gseq.extract pad = TRUE rejects non-finite coordinates as pad = FALSE does", {
    chrom1 <- as.character(gintervals.all()$chrom[1])
    # alone (no row is padded) and next to a row that is padded
    for (other in list(NULL, data.frame(chrom = chrom1, start = -5, end = 10))) {
        bad <- function(start, end) rbind(other, data.frame(chrom = chrom1, start = start, end = end))
        expect_error(gseq.extract(bad(-Inf, 10), pad = TRUE), "start coordinate must be greater or equal than zero")
        expect_error(gseq.extract(bad(10, Inf), pad = TRUE), "start coordinate must be lesser than end coordinate")
        expect_error(gseq.extract(bad(-Inf, Inf), pad = TRUE), "start coordinate must be greater or equal than zero")
        expect_error(gseq.extract(bad(NA_real_, 10), pad = TRUE), "Invalid format of interval")
        expect_error(gseq.extract(bad(-5, NaN), pad = TRUE), "Invalid format of interval")
    }
    expect_error(
        gseq.extract(data.frame(chrom = chrom1, start = 10, end = Inf), pad = FALSE),
        "start coordinate must be lesser than end coordinate"
    )
})

test_that("gseq.extract pad = TRUE keeps the genome's lowercase and N, in both database formats", {
    local_db_state()
    # one FASTA per chromosome, named after it, as the per-chromosome format needs
    fasta_dir <- tempfile()
    dir.create(fasta_dir)
    withr::defer(unlink(fasta_dir, recursive = TRUE))
    fasta <- file.path(fasta_dir, c("chrA.fa", "chrB.fa"))
    cat(">chrA\nNNacgtACGTnn\n", file = fasta[1])
    cat(">chrB\nACGTAC\n", file = fasta[2])

    for (format in c("indexed", "per-chromosome")) {
        groot <- tempfile()
        withr::defer(unlink(groot, recursive = TRUE))
        suppressMessages(gdb.create(groot = groot, fasta = fasta, format = format, verbose = FALSE))
        suppressMessages(gdb.init(groot))

        intervs <- data.frame(
            chrom = c("chrA", "chrA", "chrA", "chrB", "chrB"),
            start = c(-3, 8, -3, 6, -4),
            end = c(4, 15, 15, 9, -1),
            strand = c(1, 1, -1, 1, -1)
        )
        expect_identical(
            gseq.extract(intervs, pad = TRUE),
            c("NNNNNac", "GTnnNNN", "NNNnnACGTacgtNNNNN", "NNN", "NNN"),
            info = format
        )
    }
})

test_that("gseq.extract pad = TRUE reads an interval on an empty contig as all N", {
    local_db_state()
    fasta <- tempfile(fileext = ".fa")
    groot <- tempfile()
    withr::defer(unlink(c(fasta, groot), recursive = TRUE))
    cat(">contig1\nACTG\n>empty\n>contig2\nGGGG\n", file = fasta)
    withr::with_options(list(gmulticontig.indexed_format = TRUE), {
        suppressMessages(gdb.create(groot = groot, fasta = fasta, verbose = FALSE))
    })
    suppressMessages(gdb.init(groot))
    expect_equal(gintervals.all()$end[gintervals.all()$chrom == "empty"], 0)

    intervs <- data.frame(chrom = c("contig1", "empty", "contig2"), start = c(-1, -2, 2), end = c(2, 3, 6))
    expect_identical(gseq.extract(intervs, pad = TRUE), c("NAC", "NNNNN", "GGNN"))
})
