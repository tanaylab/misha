## R CMD check results

0 errors | 0 warnings | 0 notes

## Changes

CRAN currently serves 5.6.6; this submission is 5.11.20. `NEWS.md` has the per-version
detail. The changes that matter to users are fixes to defects that could return wrong
results or corrupt data:

* Out-of-bounds accesses and an unclamped `partial_sort` in the multitasking code paths,
  which could return silently wrong quantiles or abort the session.
* A protection-stack leak, and C++ paths that raised R warnings mid-`.Call`, so that an
  exiting handler (`tryCatch(warning = )`, `options(warn = 2)`) longjmped past a
  destructor and left the session unusable.
* Stale cache invalidation: reading a track after its database was rebuilt, or a track
  renamed, under the same path could fail or report the wrong track type.
* Argument validation across the API, where several functions accepted arguments they
  then silently discarded.

The vignettes are now evaluated at build time; the broken examples that exposed are fixed.

## Note on portability

The package implements a database that is based on shared memory files and therefore
includes many unix-specific system calls. In addition, many parallel algorithms used in
the package rely on the unix forking mechanism, therefore the package is not fully
portable to Windows.
