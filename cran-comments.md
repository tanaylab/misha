## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
* The need for the Unix OS stems from the package's implementation of a database that utilizes shared memory files, which involves numerous Unix-specific system calls. Furthermore, the package employs numerous parallel algorithms that depend on the Unix forking method.

