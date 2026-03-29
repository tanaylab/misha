# Edit Distance "Below" Direction

**Date:** 2026-03-29
**Branch:** `feat/edit-distance-below`

## Motivation

Currently, PWM edit distance functions only answer "how many edits to reach **above** a threshold?" This feature adds the symmetric question: "how many edits to bring the score **below** a threshold?" This is useful for finding minimal perturbations that disrupt a binding site.

## API Design

Add a `direction` parameter (`"above"` / `"below"`, default `"above"`) to:

**Virtual track functions:**
- `pwm.edit_distance`
- `pwm.edit_distance.pos`
- `pwm.max.edit_distance`
- `pwm.edit_distance.lse`
- `pwm.edit_distance.lse.pos`

**Sequence function:**
- `gseq.pwm_edits()`

**Semantics:**
- `direction="above"` (existing behavior): return 0 if score >= threshold, else minimum edits to reach >= threshold.
- `direction="below"`: return 0 if score <= threshold, else minimum edits to bring score <= threshold.

## Algorithm: Substitution-Only

### Core idea

Mirror the "above" algorithm. Instead of computing *gain* (how much score increases per edit), compute *loss* (how much score decreases per edit).

### Steps (direction="below")

1. Compute `surplus = current_score - threshold`. If `surplus <= 0`, return 0 (already below).
2. For each position `i`, compute `loss[i] = PSSM[i][current_base] - col_min[i]` — the maximum score degradation from switching to the worst base at that position.
3. Greedily pick positions with the highest loss until the accumulated loss covers the surplus.
4. Return the number of positions picked.

### Precomputation

- `col_min[i]` — per-column minimum PSSM score (analogous to `col_max[i]` used for "above").
- Loss bins — histogram of loss values, populated identically to gain bins but using loss values.
- `S_min = sum(col_min)` — minimum achievable score with all worst bases.

### Unreachable case

If `S_min > threshold`, the score cannot be brought below threshold even with all worst bases. Return `NA`.

### Key insight

Reuse the same bin/histogram machinery in `precompute_tables()`. Populate with loss values instead of gain values. The greedy counting logic (scan bins from largest to smallest) is identical.

## Algorithm: Indels

### Banded DP

The banded DP currently maximizes score. For "below" direction, add a "minimize score" variant:

- Fill DP cells with minimum achievable score instead of maximum.
- `early_abandon_banded_dp`: abandon if the minimum achievable score exceeds the threshold (cannot get below).

### Specialized 1-indel and 2-indel solvers

Each needs a "below" variant that:
- Considers deletions/insertions that *reduce* score most.
- Checks feasibility against `S_min` (extended for indel-shifted alignment lengths).

## Algorithm: LSE

### Objective

Find minimum edits to bring `F = log(sum(exp(S_p)))` below threshold, where `S_p` are per-position scores across strands/shifts.

### Approach

- Greedy: pick edits with the largest *negative* delta-Z (most score-reducing).
- Exhaustive search for `k <= 2`, greedy for `k >= 3` (same structure as "above", with flipped objective).
- `col_min` replaces `col_max` in precomputation; loss replaces gain.

## Files to Modify

### C++

| File | Changes |
|------|---------|
| `src/PWMEditDistanceScorer.h/cpp` | Add `Direction` enum (`ABOVE`, `BELOW`). In `precompute_tables()`, compute `col_min` and loss bins alongside existing gain bins. In scoring path, flip deficit/surplus logic based on direction. |
| `src/PWMLseEditDistanceScorer.h/cpp` | Direction support for LSE mode: flip objective in greedy/exhaustive search. |
| `src/GseqPwmEdits.cpp` | Pass direction parameter through to scorer; report edits that reduce score. |
| `src/TrackExpressionVars.h/cpp` | Parse and store direction parameter from R. |
| `src/SequenceVarProcessor.cpp` | Wire direction to scorer constructors. |

### R

| File | Changes |
|------|---------|
| `R/vtrack.R` | Accept and validate `direction` parameter in PWM edit distance vtrack functions. Pass to C++. |
| `R/sequence.R` | Accept and validate `direction` parameter in `gseq.pwm_edits()`. Pass to C++. |

### Tests

| File | Purpose |
|------|---------|
| `tests/testthat/test-pwm-edit-distance-below.R` | New test file for all "below" direction tests. |
| `tests/testthat/helper-pwm.R` | R reference implementation `manual_pwm_edit_distance_below()` for validation. |

## Test Strategy

1. **R reference implementation**: `manual_pwm_edit_distance_below()` that brute-force computes minimum edits to bring score below threshold (enumerate edits for small k, verify greedy for larger).
2. **All modes**: Test `MIN_EDITS`, `MIN_EDITS_POSITION`, `PWM_MAX_EDITS` with `direction="below"`.
3. **Edge cases**:
   - Already below threshold -> return 0.
   - Unreachable (cannot get below even with all worst bases) -> return `NA`.
   - `max_edits` cap reached -> return `NA`.
4. **Indels + below**: Test with `D=1` and `D=2` indel budgets.
5. **LSE + below**: Test LSE mode with direction="below".
6. **Regression**: All existing "above" tests must continue to pass unchanged (default direction is "above").
