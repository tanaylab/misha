## R CMD check results

0 errors | 0 warnings | 1 note

### Resubmission 2

This version contains a few bug fixes, a new function and some updated documentation. 

To answer Victoria Wimmer's question regarding the 'UNIX only' note: 

The need for the Unix OS stems from the package's implementation of a database that utilizes shared memory files, which involves numerous Unix-specific system calls. Furthermore, the package employs parallel algorithms that depend on the Unix forking method. It would require, therefore, a large and non-trivial effort to re-implement those in Windows and mark those which cannot be implemented, which is, unfortunately, not within our current bandwidth. 
