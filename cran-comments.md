## R CMD check results

0 errors | 0 warnings | 0 notes

## Changes

* Replaced non-API C entry point `Rf_findVar` with `R_getVar`.

## Note on portability

The package implements a database that is based on shared memory files and therefore includes many unix-specific system calls. In addition, many parallel algorithms used in the package rely on the unix forking mechanism, therefore the package is not fully portable to Windows.
