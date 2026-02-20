## R CMD check results

0 errors | 0 warnings | 0 notes

* Fixed devel compiler warnings in `GenomeTrackBinnedTransform.cpp` by avoiding arithmetic between different anonymous enum types.
* rhub check on `ubuntu-clang` and `gcc15` completed successfully with status OK: https://github.com/tanaylab/misha/actions/runs/22214488815

* As explained in the first submission, the need for the Unix OS stems from the package's implementation of a database that utilizes shared memory files, which involves numerous Unix-specific system calls. Furthermore, the package employs parallel algorithms that depend on the Unix forking method.
