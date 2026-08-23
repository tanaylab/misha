## R CMD check results

0 errors | 0 warnings | 0 notes

## Changes

CRAN currently serves 5.6.6. See NEWS.md for the per-version detail.

* Corrected several out-of-bounds accesses and protection-stack errors in the C++ layer.
* Corrected cache invalidation when a database is rebuilt or a track renamed.
* Argument validation across the API.
* Vignettes are now evaluated at build time.

## Note on portability

The package implements a database that is based on shared memory files and therefore
includes many unix-specific system calls. In addition, many parallel algorithms used in
the package rely on the unix forking mechanism, therefore the package is not fully
portable to Windows.
