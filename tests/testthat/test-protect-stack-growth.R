# Protection-stack growth regression tests (5.11.22).
#
# convert_intervs() used to leave seven entries on R's 50,000-slot protect stack
# per call (eleven for 2D) and never release them. Both call sites below run it in
# a loop, so the stack overflowed - an Rf_error that longjmps past ~RdbInitializer
# and leaves s_ref_count stuck for the rest of the session.
#
# Both cases are sized just past the old thresholds (measured 2026-08-24: the
# writer failed between 7,000 and 7,200 contigs, gintervals.mapply between 8,300
# and 8,330 intervals), so they fail loudly on the old code and cost ~20s here.

# These re-root into throwaway databases; put a working root back afterwards.
restore_groot_on_exit()

test_that("intervals.set.out does not grow the protect stack per contig", {
    withr::with_tempdir({
        n_contigs <- 7500
        create_test_db("manycontig", data.frame(
            chrom = paste0("chr", seq_len(n_contigs)),
            size = rep(2000, n_contigs)
        ))
        gsetroot("manycontig")
        gdb.reload()

        expect_equal(nrow(gintervals.all()), n_contigs)

        # One convert_intervs() per contig, from GIntervalsBigSet1D::save_chrom_plain_intervals
        expect_no_error(
            giterator.intervals(
                intervals = gintervals.all(), iterator = 1000,
                intervals.set.out = "manybins"
            )
        )
        expect_equal(nrow(gintervals.load("manybins")), n_contigs * 2)
    })
})

test_that("gintervals.mapply does not grow the protect stack per interval", {
    withr::with_tempdir({
        create_test_db("mapplydb", data.frame(chrom = "chr1", size = 3e6))
        gsetroot("mapplydb")
        gdb.reload()

        n <- 12000
        ivs <- gintervals(
            "chr1",
            seq(0, by = 100, length.out = n),
            seq(50, by = 100, length.out = n)
        )

        # enable.gapply.intervals converts the current intervals per scope
        # interval, so the leak was 6 slots per interval regardless of genome.
        res <- gintervals.mapply(function(x) 1, "1",
            intervals = ivs, iterator = ivs,
            enable.gapply.intervals = TRUE
        )
        expect_equal(nrow(res), n)
        expect_true(all(res$value == 1))
    })
})
