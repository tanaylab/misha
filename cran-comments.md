## R CMD check results

0 errors | 0 warnings | 1 note

### Resubmission

See response to previous submission below.

* The `misha` package, while utilized in numerous studies over the years, does not have a dedicated publication for citation.
* Added missing \value tags to all exported functions in the package.
* Replaced \dontrun with \donttest in examples.
* Replaced `cat` calls with `message`.
* Added `on.exit` calls right before any call to `setwd` or `options` as instructed.
* Removed the `options(warn=-1)` calls. 
