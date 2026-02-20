## R CMD check results

0 errors | 0 warnings | 1 note

* Fixed devel compiler warnings in `GenomeTrackBinnedTransform.cpp` by avoiding arithmetic between different anonymous enum types.
* rhub check on `ubuntu-clang` and `gcc15` completed successfully with only one note (non-standard top-level directory `conda-recipe`): https://github.com/tanaylab/misha/actions/runs/22199383282

* As explained in the first submission, the need for the Unix OS stems from the package's implementation of a database that utilizes shared memory files, which involves numerous Unix-specific system calls. Furthermore, the package employs parallel algorithms that depend on the Unix forking method.
