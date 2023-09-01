## R CMD check results

0 errors | 0 warnings | 1 note

### Resubmission

* Added missing \value tags to all exported functions in the package.
* The need for the Unix OS stems from the package's implementation of a database that utilizes shared memory files, which involves numerous Unix-specific system calls. Furthermore, the package employs numerous parallel algorithms that depend on the Unix forking method.
* The `misha` package, while utilized in numerous studies over the years, does not have a dedicated publication for citation.
* A link to macbuilder test results: https://mac.R-project.org/macbuilder/results/1693508419-aba69eca0b78a073/

