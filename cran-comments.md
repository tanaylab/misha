## R CMD check results

0 errors | 0 warnings | 1 note

### Resubmission

See response to previous submission below.

* Added Eitan Yaffe (eitany) as an author.
* Added Weizmann Institute of Science as a cph.
* The `misha` package, while utilized in numerous studies over the years, does not have a dedicated publication for citation, but I have added citation for the 2D genome algorithms. 
* Added missing \value tags to all exported functions in the package.
* Replaced \dontrun with \donttest in examples.
* Replaced `cat` calls with `message`.
* Added `on.exit` calls right before any call to `setwd` or `options` as instructed.
* Removed the `options(warn=-1)` calls. 
* Changed the example of a function that wrote files to the current working directory to a temporary directory.
