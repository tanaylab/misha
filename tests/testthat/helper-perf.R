# Helpers for the opt-in performance regression tests in test-perf-regression.R.
#
# Why the framework exists: we got bitten twice by perf regressions that
# weren't caught by correctness tests (MAP_POPULATE, validation-loop init_read).
# We also got bitten ONCE by an audit "finding" that turned out to be a false
# positive — the suspect path wasn't even hit on the workload it claimed to
# slow down. So the rule is: add a perf test only when both
# (a) it actually exercises the path we care about, and
# (b) you have measured a baseline on this hardware to anchor the threshold.
#
# Threshold rule: budget = baseline_ms * factor (default 3x). The baseline for
# each test is the FASTEST observed wall-clock for that operation across
# versions we care about — if a fix beats the historical state, use the new
# (lower) number; if any historical state was faster, use that. Never raise
# the baseline for convenience: the whole point is to ratchet the bar tight.
#
# Baselines are absolute wall-clock numbers measured on the lab box (mraid20 /
# NFS test_db) and not portable to other hardware — that's why the suite is
# opt-in via MISHA_PERF_TESTS=true.

skip_unless_perf <- function() {
    if (!identical(tolower(Sys.getenv("MISHA_PERF_TESTS")), "true")) {
        skip("perf test (set MISHA_PERF_TESTS=true to run)")
    }
}

# Run `expr` `warmup` times to warm caches and one-shot R-side overhead, then
# `measure` times for timing. Returns the median elapsed seconds — median is
# more robust than mean to the occasional GC pause or NFS hiccup.
time_op <- function(expr, warmup = 1L, measure = 3L) {
    expr_q <- substitute(expr)
    env <- parent.frame()
    for (i in seq_len(warmup)) eval(expr_q, env)
    times <- numeric(measure)
    for (i in seq_len(measure)) {
        times[i] <- system.time(eval(expr_q, env))[["elapsed"]]
    }
    median(times)
}

# Assert measured wall-clock time stays within budget = baseline_ms * factor.
expect_perf_baseline <- function(measured_s, label, baseline_ms, factor = 3) {
    budget_ms <- baseline_ms * factor
    measured_ms <- measured_s * 1000
    msg <- sprintf(
        "%s: measured %.0fms, budget %.0fms (= %.0fms baseline * %g)",
        label, measured_ms, budget_ms, baseline_ms, factor
    )
    if (measured_ms >= budget_ms) {
        fail(msg)
    } else {
        succeed(msg)
    }
}
