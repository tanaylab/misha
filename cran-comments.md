## R CMD check results

0 errors | 0 warnings | 0 notes

* Fixed memory alignment issues in clang-UBSAN (hopefully, on rhub dockers the problem is solved)
* As explained in the first submission, the need for the Unix OS stems from the package's implementation of a database that utilizes shared memory files, which involves numerous Unix-specific system calls. Furthermore, the package employs parallel algorithms that depend on the Unix forking method.


