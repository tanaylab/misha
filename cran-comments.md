## R CMD check results

0 errors | 0 warnings | 0 notes

## Changes

CRAN currently serves 5.6.6. See NEWS.md for the per-version detail.

* Corrected out-of-bounds accesses, protection-stack errors and uninitialised reads in the C++ layer.
* Native routines are now registered, so argument counts are checked at call time.
* Track writes are staged, and a failed or interrupted write no longer leaves a partial track.
* Ctrl-C now interrupts long conversions and sequence scans.
* Argument validation across the API, and the vignettes are evaluated at build time.

## Note on portability

The package implements a database that is based on shared memory files and therefore
includes many unix-specific system calls. In addition, many parallel algorithms used in
the package rely on the unix forking mechanism, therefore the package is not fully
portable to Windows.
