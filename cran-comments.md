## R CMD check results

### Local

`R -e "devtools::check()"`

Result: 1 ERROR | 1 WARNING | 5 NOTEs

Notes:
* The ERROR comes from `testthat` tests that run only when `NOT_CRAN=true` (set by `devtools::check()`), including environment/performance-sensitive tests and a `prego` API mismatch.
* These tests are guarded in `tests/testthat.R` and are not run on CRAN.

### rhub (GitHub Actions backend)

Run: https://github.com/tanaylab/misha/actions/runs/22199383282

* `ubuntu-clang`: 0 errors | 0 warnings | 1 note
  - Note: non-standard top-level directory `conda-recipe`
* `gcc15`: 0 errors | 0 warnings | 1 note
  - Note: non-standard top-level directory `conda-recipe`

The previously reported devel compiler enum-conversion warnings in `GenomeTrackBinnedTransform.cpp` are no longer present on these devel targets.

* As explained in the first submission, the need for the Unix OS stems from the package's implementation of a database that utilizes shared memory files, which involves numerous Unix-specific system calls. Furthermore, the package employs parallel algorithms that depend on the Unix forking method.

