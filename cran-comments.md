## R CMD check results

0 errors | 0 warnings | 0 notes

* Fixed *noRemap* additional issue.
* Removed non-API calls to R.
* Fixed clang-UBSAN error.
* As explained in the first submission, the need for the Unix OS stems from the package's implementation of a database that utilizes shared memory files, which involves numerous Unix-specific system calls. Furthermore, the package employs parallel algorithms that depend on the Unix forking method.


