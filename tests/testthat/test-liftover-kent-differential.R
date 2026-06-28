# Randomized differential test: gintervals.liftover() vs UCSC liftOver.
#
# For each seed it generates several single-block chains with DISJOINT source and
# target ranges (so there is exactly one valid mapping and the two tools must agree
# bit-for-bit - no chain-selection or fragment-merging ambiguity), random + and -
# strand, then lifts a batch of random and adversarial query intervals through both
# misha and `liftOver -multiple` and asserts the canonicalized target coverage per
# query is identical.
#
# Scope: this guards the coordinate-mapping core (block boundaries, strand reversal,
# partial overlaps, and the multi-interval carried-cursor path in map_interval)
# against regressions. Validated to have teeth: a 1bp offset in the mapping produces
# hundreds of mismatches here; correct code produces zero. Multi-block (gapped)
# chains are intentionally excluded because misha fragments at gaps while liftOver
# bridges them - a documented, legitimate difference, not a bug.
#
# Skips when the liftOver binary is unavailable (e.g. CI without UCSC tools).

skip_if_not(has_liftover_binary(), "liftOver binary not found")

# merge sorted/overlapping/touching intervals -> canonical coverage (single chrom)
.canon_cov <- function(df) {
    if (is.null(df) || nrow(df) == 0) {
        return(data.frame(start = integer(), end = integer()))
    }
    df <- df[order(df$start, df$end), , drop = FALSE]
    s <- df$start[1]
    e <- df$end[1]
    out <- list()
    for (i in seq_len(nrow(df))[-1]) {
        if (df$start[i] <= e) {
            e <- max(e, df$end[i])
        } else {
            out[[length(out) + 1]] <- c(s, e)
            s <- df$start[i]
            e <- df$end[i]
        }
    }
    out[[length(out) + 1]] <- c(s, e)
    r <- as.data.frame(do.call(rbind, out))
    names(r) <- c("start", "end")
    r
}

test_that("gintervals.liftover matches UCSC liftOver on randomized single-block chains", {
    local_db_state()

    SRC <- 200000L
    TGT <- 300000L

    # Target genome only - the source chromosome is just a coordinate label carried
    # in the chain file and the query intervals.
    gdir <- tempfile("kentdiff_")
    dir.create(gdir)
    fa <- file.path(gdir, "chr1.fasta")
    cat(">chr1\n", paste(rep("A", TGT), collapse = ""), "\n", sep = "", file = fa)
    tdb <- tempfile("kentdiff_db_")
    suppressMessages(gdb.create(groot = tdb, fasta = fa))
    withr::defer(
        {
            unlink(tdb, recursive = TRUE)
            unlink(gdir, recursive = TRUE)
        },
        envir = testthat::teardown_env()
    )
    gdb.init(tdb)

    # write one single-block chain mapping source [ss,se) onto forward target
    # [tf_s, tf_e); for "-" the chain stores reverse-strand query coordinates.
    write_block <- function(con, id, score, ss, se, tf_s, strand) {
        tf_e <- tf_s + (se - ss)
        if (strand == "+") {
            qs <- tf_s
            qe <- tf_e
        } else {
            qs <- TGT - tf_e
            qe <- TGT - tf_s
        }
        write_chain_entry(con, "source1", SRC, "+", ss, se, "chr1", TGT, strand, qs, qe, id, score = score)
    }

    total_mismatch <- 0L
    for (seed in 1:25) {
        set.seed(seed)
        cf <- new_chain_file()
        nchain <- sample(2:5, 1)
        ss <- 1000L
        ts <- 1000L
        src_ranges <- list()
        for (k in seq_len(nchain)) {
            blk <- as.integer(sample(40:600, 1))
            strand <- sample(c("+", "-"), 1)
            write_block(cf, k, sample(1000:9000, 1), ss, ss + blk, ts, strand)
            src_ranges[[k]] <- c(ss, ss + blk)
            ss <- ss + blk + sample(200:3000, 1) # disjoint source
            ts <- ts + blk + sample(200:3000, 1) # disjoint target
        }

        # queries: random + adversarial (global wide span, per-chain wide, narrow at
        # each chain start) to exercise the carried-cursor path in map_interval.
        rq_s <- as.integer(sample(0:(ss + 1000), 40, replace = TRUE))
        rq_w <- as.integer(sample(1:400, 40, replace = TRUE))
        q <- data.frame(chrom = "source1", start = rq_s, end = pmin(rq_s + rq_w, SRC), stringsAsFactors = FALSE)
        last_se <- src_ranges[[nchain]][2]
        adv <- data.frame(chrom = "source1", start = 0L, end = min(SRC, last_se + 50L), stringsAsFactors = FALSE)
        for (r in src_ranges) {
            adv <- rbind(
                adv,
                data.frame(chrom = "source1", start = c(max(0L, r[1] - 50L), r[1]), end = c(min(SRC, r[2] + 50L), min(SRC, r[1] + 20L)), stringsAsFactors = FALSE)
            )
        }
        q <- rbind(q, adv)
        q <- q[q$end > q$start, , drop = FALSE]
        q$qid <- seq_len(nrow(q))

        ch <- gintervals.load_chain(cf) # defaults: src=error, tgt=auto
        mish <- gintervals.liftover(q[, c("chrom", "start", "end")], ch, include_metadata = TRUE)

        bed <- tempfile(fileext = ".bed")
        cat(sprintf("source1\t%d\t%d\tq%d\n", q$start, q$end, q$qid), file = bed)
        # Call liftOver directly (the run_kent_liftover helper mis-parses the unmapped
        # file when the BED carries a name column); we only need the mapped output.
        outf <- tempfile()
        unf <- tempfile()
        system2("liftOver", c("-multiple", "-minMatch=0.0000001", bed, cf, outf, unf), stdout = FALSE, stderr = FALSE)
        kent <- if (file.exists(outf) && file.info(outf)$size > 0) {
            read.table(outf, col.names = c("chrom", "start", "end", "name", "val"), stringsAsFactors = FALSE)
        } else {
            data.frame(chrom = character(), start = integer(), end = integer(), name = character(), val = integer())
        }
        unlink(c(outf, unf))

        for (id in q$qid) {
            m <- if (!is.null(mish) && nrow(mish)) mish[mish$intervalID == id, c("start", "end"), drop = FALSE] else data.frame()
            k <- if (nrow(kent)) kent[kent$name == paste0("q", id), c("start", "end"), drop = FALSE] else data.frame()
            if (!isTRUE(all.equal(.canon_cov(m), .canon_cov(k), check.attributes = FALSE))) {
                total_mismatch <- total_mismatch + 1L
            }
        }
        unlink(bed)
    }

    expect_equal(total_mismatch, 0L)
})
