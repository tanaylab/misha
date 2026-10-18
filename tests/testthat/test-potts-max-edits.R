create_isolated_test_db()

# The best anchor score over the sequence and every sequence one substitution
# away, from the pure-R oracle in helper-potts.R. It substitutes over the whole
# read region and maximises over all anchors, so it does not rely on the
# per-anchor decomposition the kernel uses.
potts_one_edit_oracle <- function(seq, model, bidirect, strand) {
    score <- function(s) {
        v <- potts_ref_anchors(s, model, bidirect = bidirect, strand = strand)
        if (all(is.na(v))) -Inf else max(v, na.rm = TRUE)
    }
    best <- score(seq)
    chars <- strsplit(toupper(seq), "")[[1]]
    for (j in seq_along(chars)) {
        if (!(chars[j] %in% POTTS_BASES)) next
        for (b in setdiff(POTTS_BASES, chars[j])) {
            s2 <- chars
            s2[j] <- b
            best <- max(best, score(paste(s2, collapse = "")))
        }
    }
    if (is.finite(best)) best else NA_real_
}

test_that("potts.max with max_edits = 1 matches the brute-force oracle", {
    cases <- expand.grid(
        W = c(5L, 7L), mode = c("full", "sparse", "none"),
        bidirect = c(TRUE, FALSE), strand = c(1L, -1L), stringsAsFactors = FALSE
    )
    cases <- cases[!(cases$bidirect & cases$strand == -1L), ]
    for (r in seq_len(nrow(cases))) {
        cs <- cases[r, ]
        m <- potts_ref_model(W = cs$W, npair_mode = cs$mode, seed = 200L + r)
        gvtrack.create("p1", NULL, "potts.max",
            params = c(m, list(bidirect = cs$bidirect, strand = cs$strand, max_edits = 1))
        )
        for (len in c(1L, 9L)) {
            start <- 5000 + 37 * r
            iv <- gintervals(1, start, start + len)
            got <- gextract("p1", iv, iterator = iv)$p1
            seq <- gseq.extract(gintervals(1, start, start + len + cs$W - 1L))
            expect_equal(got, potts_one_edit_oracle(seq, m, cs$bidirect, cs$strand),
                tolerance = 1e-5, info = sprintf("case %d, len %d", r, len)
            )
        }
        gvtrack.rm("p1")
    }
})

test_that("max_edits = 0 is potts.max and max_edits = 1 never scores below it", {
    m <- potts_ref_model(W = 6L, npair_mode = "full", seed = 11L)
    gvtrack.create("p0", NULL, "potts.max", params = c(m, list(bidirect = TRUE)))
    gvtrack.create("p0e", NULL, "potts.max", params = c(m, list(bidirect = TRUE, max_edits = 0)))
    gvtrack.create("p1", NULL, "potts.max", params = c(m, list(bidirect = TRUE, max_edits = 1)))
    x <- gextract(c("p0", "p0e", "p1"), gintervals(1, 10000, 12000), iterator = 1)
    expect_equal(x$p0e, x$p0)
    expect_true(all(x$p1 >= x$p0 - 1e-5, na.rm = TRUE))
    expect_gt(mean(x$p1 > x$p0 + 1e-6, na.rm = TRUE), 0.5)
})

test_that("max_edits is validated", {
    m <- potts_ref_model(W = 5L, npair_mode = "full", seed = 5L)
    expect_error(gvtrack.create("bad", NULL, "potts.max", params = c(m, list(max_edits = 2))), "0 or 1")
    expect_error(gvtrack.create("bad", NULL, "potts.max", params = c(m, list(max_edits = NA))), "0 or 1")
    for (f in c("potts", "potts.max.pos")) {
        expect_error(gvtrack.create("bad", NULL, f, params = c(m, list(max_edits = 1))), "potts.max only")
    }
    expect_error(
        gvtrack.create("bad", NULL, "potts.count", params = c(m, list(score.thresh = 0, max_edits = 1))),
        "potts.max only"
    )
})
