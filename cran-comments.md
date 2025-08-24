## R CMD check results

0 errors | 0 warnings | 0 notes

* Removed C++14 requirement.
* Fixed a crucial bug in `pwm` and `kmer` virtual track functions: iterator shifts were not applied. 
* Adressed CRAN rchk notes and improved memory protection.
* Removed non-API calls to Rf_GetOption and Rf_isFrame.
* As explained in the first submission, the need for the Unix OS stems from the package's implementation of a database that utilizes shared memory files, which involves numerous Unix-specific system calls. Furthermore, the package employs parallel algorithms that depend on the Unix forking method.


