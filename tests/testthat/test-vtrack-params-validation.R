create_isolated_test_db()

# Regression tests for two silent-drop holes in gvtrack.create():
#   1. `.vtrack_params_*` handlers read only the keys they know about from
#      `params`/`...` without validating the key set, so a misspelled key
#      (e.g. bidrect instead of bidirect) is silently dropped and the
#      default value is used instead - producing a materially different,
#      wrong result with no warning.
#   2. gvtrack.create() only inspects `...` when `func` has a handler in
#      .VTRACK_PARAM_HANDLERS; for every other func (e.g. "avg"), any named
#      arguments passed via `...` (such as a misspelled `sshift`) are
#      dropped whole.

test_pssm <- function() {
    pssm <- matrix(
        c(
            0.7, 0.1, 0.1, 0.1,
            0.1, 0.7, 0.1, 0.1,
            0.1, 0.1, 0.1, 0.7
        ),
        ncol = 4, byrow = TRUE
    )
    colnames(pssm) <- c("A", "C", "G", "T")
    pssm
}

test_that("hole 1: misspelled params-list key for pwm errors instead of silently defaulting", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- test_pssm()
    scope <- gintervals(1, 0, 300)

    # Correctly spelled bidirect = FALSE must produce a different result than
    # the (silently-defaulted) bidirect = TRUE default, proving the parameter
    # actually has an effect - otherwise this test wouldn't catch a drop.
    gvtrack.create("correct_true", NULL, "pwm", params = list(pssm = pssm, bidirect = TRUE))
    gvtrack.create("correct_false", NULL, "pwm", params = list(pssm = pssm, bidirect = FALSE))
    res <- gextract(c("correct_true", "correct_false"), scope, iterator = 100)
    expect_false(isTRUE(all.equal(res$correct_true, res$correct_false)))

    # The one-character typo must now error, not silently fall back to
    # bidirect = TRUE (which is what happened before the fix).
    expect_error(
        gvtrack.create("typo", NULL, "pwm", params = list(pssm = pssm, bidrect = FALSE)),
        "bidrect"
    )
    expect_error(
        gvtrack.create("typo", NULL, "pwm", params = list(pssm = pssm, bidrect = FALSE)),
        "pwm"
    )
})

test_that("hole 1: misspelled named-arg (...) key for pwm errors instead of silently defaulting", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- test_pssm()

    expect_error(
        gvtrack.create("typo2", NULL, "pwm", pssm = pssm, bidrect = FALSE),
        "bidrect"
    )
})

test_that("hole 2: unknown '...' argument for a func with no handler now errors", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Verified reproducer from the review: a typo of the real formal
    # 'sshift' (here 'sshfit') plus a nonsense key, for func = 'avg' which
    # has no entry in .VTRACK_PARAM_HANDLERS.
    expect_error(
        gvtrack.create("v1", "test.fixedbin", "avg", sshfit = 100, nonsense = "x"),
        "avg"
    )
    expect_false("v1" %in% gvtrack.ls())
})

test_that("hole 2: real formals (dim/sshift/eshift/filter) for a func with no handler still work", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    expect_silent(gvtrack.create("v2", "test.fixedbin", "avg", sshift = 100, eshift = 50))
    info <- gvtrack.info("v2")
    expect_equal(info$itr$sshift, 100)
    expect_equal(info$itr$eshift, 50)
})

test_that("every handler's currently-documented keys still work after validation", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- test_pssm()

    # pwm family (pwm, pwm.max, pwm.max.pos, pwm.count): pssm, bidirect,
    # prior, extend, strand, spat_factor, spat_bin, spat_min, spat_max,
    # score.thresh
    expect_silent(gvtrack.create("v_pwm", NULL, "pwm",
        pssm = pssm, bidirect = FALSE, prior = 0.01, extend = TRUE, strand = 1,
        spat_factor = c(1, 1), spat_bin = 40L, spat_min = 1, spat_max = 100
    ))
    expect_silent(gvtrack.create("v_pwm_count", NULL, "pwm.count", pssm = pssm, score.thresh = -5))

    # kmer family: kmer, extend, strand
    expect_silent(gvtrack.create("v_kmer", NULL, "kmer.count", kmer = "AT", extend = TRUE, strand = 1))
    expect_silent(gvtrack.create("v_kmer2", NULL, "kmer.frac", params = list(kmer = "AC", strand = 0)))

    # pwm.edit_distance family: + score.thresh, max_edits, max_indels,
    # score.min, score.max, direction
    expect_silent(gvtrack.create("v_edist", NULL, "pwm.edit_distance",
        pssm = pssm, score.thresh = -5, max_edits = 2, max_indels = 1,
        score.min = -10, score.max = 0, direction = "above"
    ))

    # pwm.n_mutations: no max_edits/max_indels
    expect_silent(gvtrack.create("v_nmut", NULL, "pwm.n_mutations",
        pssm = pssm, score.thresh = -5, score.min = -10, score.max = 0, direction = "below"
    ))

    # pwm.edit_distance.lse: max_edits but no max_indels
    expect_silent(gvtrack.create("v_lse", NULL, "pwm.edit_distance.lse",
        pssm = pssm, score.thresh = -5, max_edits = 2, score.min = -10, score.max = 0
    ))

    # neighbor.count: plain numeric params, no named dots
    src <- gintervals(1, 100, 110)
    expect_silent(gvtrack.create("v_near", src, "neighbor.count", 10))
})

test_that("pwm.n_mutations and pwm.edit_distance.lse reject max_indels (not in their key set)", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- test_pssm()

    expect_error(
        gvtrack.create("v_bad1", NULL, "pwm.n_mutations", pssm = pssm, score.thresh = -5, max_indels = 1),
        "max_indels"
    )
    expect_error(
        gvtrack.create("v_bad2", NULL, "pwm.edit_distance.lse", pssm = pssm, score.thresh = -5, max_indels = 1),
        "max_indels"
    )
})

test_that("neighbor.count rejects unknown named arguments via '...'", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    src <- gintervals(1, 100, 110)
    expect_error(
        gvtrack.create("v_near_bad", src, "neighbor.count", 10, max_distance = 5),
        "max_distance"
    )
})

test_that("mixing 'params' and named '...' arguments in one call errors", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- test_pssm()

    # Every handler does `dots <- params` when params is non-NULL, so the
    # named arguments never reached the unknown-key check and were silently
    # dropped: this call used to return the bidirect = TRUE numbers.
    expect_error(
        gvtrack.create("v_mixed", NULL, "pwm", params = list(pssm = pssm), bidirect = FALSE),
        "bidirect"
    )
    expect_error(
        gvtrack.create("v_mixed", NULL, "pwm", params = list(pssm = pssm), bidirect = FALSE),
        "params"
    )
    expect_false("v_mixed" %in% gvtrack.ls())

    # Same hole in the kmer handler (params wins over dots there too).
    expect_error(
        gvtrack.create("v_mixed2", NULL, "kmer.count", params = list(kmer = "AC"), strand = 1),
        "strand"
    )

    # Either style alone still works, and they still agree with each other.
    scope <- gintervals(1, 0, 300)
    expect_silent(gvtrack.create("v_p", NULL, "pwm", params = list(pssm = pssm, bidirect = FALSE)))
    expect_silent(gvtrack.create("v_d", NULL, "pwm", pssm = pssm, bidirect = FALSE))
    res <- gextract(c("v_p", "v_d"), scope, iterator = 100)
    expect_equal(res$v_p, res$v_d)
})

test_that("masked.count/masked.frac extra-parameter behavior (warning) is unchanged", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # masked.* already validated its (empty) key set via a warning before
    # this fix; that pre-existing, documented behavior is left untouched.
    expect_warning(
        gvtrack.create("v_masked", NULL, "masked.count", extra_param = 123),
        "do not accept parameters"
    )
    expect_true("v_masked" %in% gvtrack.ls())
})

# ---------------------------------------------------------------------------
# score.thresh is validated the same way across the pwm family
#
# pwm.count coerces score.thresh through .coerce_score_thresh(): exactly one
# value, numbers and character/factor spellings of numbers accepted, anything
# else rejected by name. The pwm.edit_distance family compares score.thresh
# against the same PWM log-likelihood scale - the LSE variant against a
# log-sum-exp of those scores, which can be positive where a single score
# cannot, but neither the helper nor the C++ layer constrains the range - so
# its accepted-value contract is the same one, and it shares the helper.
#
# Before that it used an is.numeric() guard of its own, which rejected the
# character and factor thresholds a config file or a read.csv column produces,
# and let an integer through to the C++ layer, which requires a double and
# rejected it as "not numeric".
# ---------------------------------------------------------------------------

# 'thresh' must be a value at which 'func' returns something that varies over
# the scope: a column saturated at 0, or at NA because the threshold is out of
# reach, compares equal to anything and could not fail.
expect_score_thresh_coerced <- function(func, thresh, iterator, scope = gintervals(1, 200, 600)) {
    values <- function(score.thresh) {
        remove_all_vtracks()
        gvtrack.create("v_thresh", NULL, func, pssm = test_pssm(), score.thresh = score.thresh)
        gextract("v_thresh", scope, iterator = iterator)$v_thresh
    }

    numeric_values <- values(thresh)
    expect_gt(length(unique(numeric_values)), 1) # not saturated at a single value
    expect_false(anyNA(numeric_values)) # not saturated at "unreachable"
    expect_true(all(numeric_values >= 0))

    # A threshold read out of a config file or a read.csv column arrives as a
    # character, or as a factor if stringsAsFactors was left on; as.numeric()
    # on a factor returns the level index (1 here), which is a different
    # threshold and would give a different column.
    expect_equal(values(as.character(thresh)), numeric_values)
    expect_equal(values(factor(as.character(thresh))), numeric_values)
    # An integer used to reach C++ as an INTSXP, which rejected it.
    expect_equal(values(as.integer(thresh)), numeric_values)

    expect_error(values(c(thresh, thresh + 1)), "score.thresh must be a single value")
    expect_error(values(NULL), "score.thresh must be a single value, and this one is empty")
    expect_error(values(TRUE), "score.thresh must be a single number")
    expect_error(values("loose"), "score.thresh must be a single number")
    expect_error(values(NA_real_), "score.thresh must be a single number")
}

test_that("pwm.edit_distance takes the same score.thresh values as pwm.count", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # -3 leaves both 0-edit and 1-edit bins in the scope.
    expect_score_thresh_coerced("pwm.edit_distance", -3, iterator = 20)
})

test_that("pwm.n_mutations takes the same score.thresh values as pwm.count", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # At the motif length, -5 leaves both 0-mutation and 3-mutation bins.
    expect_score_thresh_coerced("pwm.n_mutations", -5, iterator = 3)
})

test_that("pwm.edit_distance.lse takes the same score.thresh values as pwm.count", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    expect_score_thresh_coerced("pwm.edit_distance.lse", 0, iterator = 20)

    # The LSE objective sums over windows, so unlike a single PWM score it can
    # be positive: a positive threshold is an ordinary value here, and a
    # fractional one has to survive coercion exactly.
    remove_all_vtracks()
    gvtrack.create("lse_num", NULL, "pwm.edit_distance.lse", pssm = test_pssm(), score.thresh = 0.5)
    gvtrack.create("lse_chr", NULL, "pwm.edit_distance.lse", pssm = test_pssm(), score.thresh = "0.5")
    res <- gextract(c("lse_num", "lse_chr"), gintervals(1, 200, 600), iterator = 20)
    expect_gt(length(unique(res$lse_num)), 1)
    expect_true(all(res$lse_num > 0))
    expect_equal(res$lse_chr, res$lse_num)
})
