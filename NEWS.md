# misha 4.1.1

* Use roxygen2 for documentation 
* Changed gsetroot/gdb.init to work with `devtools::load_all`
* Added a `NEWS.md` file to track changes to the package.

# misha 4.1.0

* Major bug fix (bug first appears in 4.0.8, released on June 5, 2020): virtual tracks based on `global.percentile`, `global.percentile.min` or `global.percentile.max` might occasionally return unexpected results and/or cause crashes due to faulty memory management.

# misha 4.0.9

* Fixed a minor resource leak.
* Redirect all messages and progress to stderr instead of stdout.

# misha 4.0.8

* Fixed a resource leak that might result in "protection stack overflow" error.

# misha 4.0.6

* Crash fix in gdb.create / gintervals.import.genes.

# misha 4.0.5

* Increase the maximal number of tracks allowed in a track expression to 10.000.

# misha 4.0.4

* Bug fix: "child process ended unexpectedly" errors, crashes and hang ups whilst multitasking when running out of memory

# misha 4.0.3

* Fixed installation issues on some platforms

# misha 4.0.2 

* Switched from custom random seed control (options(grnd.seed=...)) to R standard (set.seed)

# misha 4.0.1 

* Bug fix in gsample: results differ on Linux vs. OSX even when the same random seed is used on both platforms

# misha 4.0.0 

* OSX support
* Bug fix in all functions using 2D intervals set iterator: in multitasking mode some of the chromosome pairs might be skipped. In non-multitasking mode the scope might be considered as empty. This behavior is random and recurrent calls might suddenly return correct results.
* Bug fix in multitasking: occasional hang ups when memory usage of the child processes exceeds the limit gmax.mem.usage
* Bug fix in multitasking: all functions creating new files and reporting progress might create corrupted files (gtrack.create, ....)
* Bug fix in all functions using intervals.set.out parameter: small intervals set instead of big one might be created and vice versa. Alsoan error message "result size exceeded the maximum allowed" might be mistakenly generated
* Bug fix in gintervals.load applied to 1D big intervals: "Error in if (progress && percentage < 100 && progress.percentage != percentage)"
* Bug fix in all functions returning 1D or 2D intervals: in rare random cases NULL or invalid intervals set is returned
* Bug fix in gcompute_strands_autocorr: internal buffer overflow and possible memory corruption
* gdb.reload: run-time improvements 

# misha 3.7.1

* Bug fix in gintervals.liftover and gtrack.liftover: some intervals might fail to be translated
* Bug fix in gintervals.liftover: 'object 'f' not found' error if chain intervals are used in 'chain' parameter

# misha 3.7.0 

* Ubuntu support (multitasking mode still not thorougly tested)
* Bug fix: no progress report in multitasking mode
* Bug fix in gdb.create: temporary directory is created under the current GROOT instead of the new one

# misha 3.6.0 

* Bug fix: occasional defunc processes AND/OR hanging in multitasking mode

# misha 3.5.6 

* Run-time optimizations

# misha 3.5.5

* Avoid call to gdb.reload() (slow on large DB) in various functions that create or remove tracks or intervals sets
* Bug fix in gintervals.neighbors: with 2D intervals the number of the returned neighbors might be less than "maxneighbors" parameter
* Bug fix in gintervals.neighbors: with 2D intervals NULL might be returned instead of NA if na.if.notfound=T

# misha 3.5.4 

* Improved control over total maximal memory consumption in multitasking mode via gmax.mem.usage option.

# misha 3.5.3 

* Run time optimizations when using several virtual tracks based on the same array track, differing only by slice

# misha 3.5.2

* Allow usage of sparse / arrays tracks in place of intervals
* Allow usage of big intervals sets in gintervals.diff, gintervals.intersect, gintervals.mapply, gintervals.union
* Run time optimizations when 1D big intervals set is used for scope
* Run time optimizations in various gintervals.* functions when big intervals sets are used
* Bug fix: functions might get stuck or crash when array track / sparse track / 1D big intervals set iterator is used along with 1D big * intervals scope
* Bug fix: functions might get stuck or crash when 2D big intervals set iterator is used along with 2D big intervals scope
* Bug fix in gintervals.neighbors and gintervals.intersect: result might be poorly sorted when big interval sets are used
* Bug fix in gintervals.load: on an empty set or chrom returns NULL for small intervals sets and an empty data frame for a big intervals set
* Bug fix: "100%..." or "100%" is sometimes printed as the only progress report
* Bug fix: multiple progress report in some functions

# misha 3.5.1 

* Bug fix in various gintervals.* functions: invalid output (except for the first row) when 2D track is used for intervals
* Bug fix in gextract: incorrect intervalID returned when 2D track is used for intervals

# misha 3.5.0 

* Allow usage of 2D track in place of intervals
* Add progress report to gintervals.load
* Bug fix when using 2D big intervals: in some cases some or all chromosomes of the big intervals set might be skipped when big intervals set is used as a scope
* Big intervals set: before load verify that the size of a single chromosome (or chromosome pair) does not exceed gmax.data.size
* Bug fix in gcluster.run: clean up of running processes might not be completed if Ctr+C is pressed multiple times
* Bug fix in gintervals.load: returns all 2D intervals instead of a subset if one of chrom1/chrom2 is NULL and another one is not NULL
* Bug fix in gintervals.load: invalid row names if chrom / chrom1 / chrom2 parameter is used for a small intervals set

# misha 3.4.3

* Support intervals represented by tibbles

# misha 3.4.2 

* New function: gsample; returns N random samples from the specified track expression
* Improved random seed when options(grnd.seed=0): so far two calls occurring within a second used identical random generators

# misha 3.4.1

* Fixed compilation errors on some platforms
* Run time improvement in gintervals.intersect when big intervals sets are used
* Bug fix in gintervals.neighbors: returns NULL if 2D big intervals sets are used
* Bug fix: a few point tracks in a track expression might be used without specifying an iterator

# misha 3.4.0 

* Dynamically limit memory use in multitasking mode
* Bug fix: race condition and potential crash in multitasking mode when one of the child processes exits shortly after it is launched

# misha 3.3.9 

* New function: gintervals.2d.import
* New option: gbig.intervals.size - controls the threshold when big intervals sets are created. Default value: 1000000
* Bug fix in gintervals.2d.import_contacts: utterly huge tracks might have missing areas of contacts
* Reduced the default value of gmax.processes option from 64 to 16
* Bug fix: incorrect progress report in gtrack.2d.import_contacts

# misha 3.3.8 

* New function: gintervals.rbind. Runs rbind on intervals sets including big intervals sets on disk
* New "intervals.set.out" parameter added to: gtrack.array.extract, gwilcox
* Support big intervals sets in: gseq.extract, gtrack.array.extract, gtrack.modify
* Bug fix: gtrack.2d.import_contacts does not recognize chromosomes that have "chr" prefix

# misha 3.3.7 

* Bug fix for all functions using 2D iterators: in some cases full chromosome pairs can be skipped. Bug first appeared in 3.3.0

# misha 3.3.6

* Bug fix in gdb.create: "Error in .gintervals.check_new_set(intervals.set.out)..."

# misha 3.3.5 

* Added 'opt.flags' parameter to gcluster.run. Use this parameter to add restrictions to the machines that run submitted jobs: minimal RAM * requirement, explicit hostnames list, etc. See man for qsub, "-l" flag.
* Support big intervals sets in the following functions: gpartition
* New "intervals.set.out" parameter added to: gpartition, gsegment
* Bug fix: invalid error recovery in glookup - traces from intervals.set.out might be left

# misha 3.3.4 

* Interface change: gintervals.neighbors returns a data frame containing full intervals instead of their ids
* gintervals.neighbors: colnames parameter removed
* Support SAM files in gtrack.import_mappedseq
* Support big intervals sets in the following functions: gintervals.neighbors, glookup
* New "intervals.set.out" parameter added to: gintervals.neighbors, glookup
* Removed gintervals.merge function
* Runtime optimizations when big intervals sets are used in the following functions: gintervals.annotate, gintervals.diff, gintervals.* force_range, gintervals.intersect, gintervals.mapply, gintervals.save, gintervals.union

# misha 3.3.3 

* Bug fix in gquantiles, multitasking version (which is used by default): invalid quantiles might be returned if the number of iterator * intervals exceeds gmax.data.size / number_of_child_processes. number_of_child_processes equals at most to the number of different chromosomes (or chromosome pairs for 2D) used in in "intervals" parameter
* Support big intervals sets in the following functions: gintervals.mapply
* New "intervals.set.out" parameter added to: gintervals.mapply

# misha 3.3.2 

* Support big intervals sets in the following functions: gintervals.annotate, gintervals.diff, gintervals.intersect, gintervals.union, * gsegment, gwilcox
* New "intervals.set.out" parameter added to: gintervals.annotate, gintervals.diff, gintervals.intersect, gintervals.union
* Added progress report to: gintervals.save, gintervals.force_range

# misha 3.3.1 

* Support big intervals sets in the following functions: gcis_decay (only intervals parameter), gintervals.2d.band_intersect
* New "intervals.set.out" parameter added to: gintervals.2d.band_intersect
* Bug fix in gintervals.force_range: "Error in if (size > max.data.size) { : argument is of length zero"

# misha 3.3.0 

* New concept: big intervals sets
* Big intervals sets can be used for iterator parameter in all functions
* Big intervals sets can be used in intervals parameter in the following functions: gdist, gextract, gquantiles, gscreen, gsummary, gbins.* quantiles, gbins.summary, gintervals.quantiles, gintervals.summary, giterator.intervals
* Interface change: gintervals.quantiles, gintervals.summary now return also the source intervals
* New functions: gintervals.is.bigset, gintervals.chrom_sizes, gintervals.update
* New "chrom", "chrom1", "chrom2", parameters for gintervals.load
* New "intervals.set.out" parameter added to: gextract, gscreen, gintervals.force_range, gintervals.quantiles, gintervals.summary, giterator.intervals
* Restrict gintervals.quantiles and gintervals.mapply to work with only 1D and Fixed Rectangle iterators
* Changed the position of "file" parameter in gextract
* Bug fix: gscreen on vtrack with global.percentile.max returns different number of intervals in each run
* Bug fix: in 2D iterators progress report can sometimes go backwards
* Bug fix: empty intervals set is ignored if used as an iterator
* Bug fix: crash if an empty intervals set is used for scope
* Bug fix: memory leak in giterator.cartesian_grid

# misha 3.2.7

* Bug fix: clean up on error when P-values table is loaded (new RSaneSerialize and RSaneUnserialize functions)

# misha 3.2.6

* New function: gcis_decay
