# misha 4.1.0

1.  Bug fix (bug first appears in 4.0.8): virtual tracks based on "global.percentile", "global.percentile.min" or "global.percentile.max" might occasionally return unexpected results and/or cause crashes due to faulty memory management.

# misha 4.0.11

1.  Fixed compilation errors on OSX.

# misha 4.0.10

1.  "unprotect_ptr: pointer not found" error in virtual tracks based on intervals.

# misha 4.0.9

1.  Fixed a minor resource leak.
2.  Redirect all messages and progress to stderr instead of stdout.

# misha 4.0.8

1.  Fixed a resource leak that might result in "protection stack overflow" error.

# misha 4.0.6

1.  Crash fix in gdb.create / gintervals.import.genes.

# misha 4.0.5

1.  Increase the maximal number of tracks allowed in a track expression to 10.000.

# misha 4.0.4

1.  Bug fix: "child process ended unexpectedly" errors, crashes and hang ups whilst multitasking when running out of memory

# misha 4.0.3

1.  Fixed installation issues on some platforms

# misha 4.0.2

1.  Switched from custom random seed control (options(grnd.seed=...)) to R standard (set.seed)

# misha 4.0.1

1.  Bug fix in gsample: results differ on Linux vs. OSX even when the same random seed is used on both platforms

# misha 4.0.0

1.  OSX support
2.  Bug fix in all functions using 2D intervals set iterator: in multitasking mode some of the chromosome pairs might be skipped. In non-multitasking mode the scope might be considered as empty. This behavior is random and recurrent calls might suddenly return correct results.
3.  Bug fix in multitasking: occasional hang ups when memory usage of the child processes exceeds the limit gmax.mem.usage
4.  Bug fix in multitasking: all functions creating new files and reporting progress might create corrupted files (gtrack.create, ....)
5.  Bug fix in all functions using intervals.set.out parameter: small intervals set instead of big one might be created and vice versa. Also an error message "result size exceeded the maximum allowed" might be mistakenly generated
6.  Bug fix in gintervals.load applied to 1D big intervals: "Error in if (progress && percentage < 100 && progress.percentage != percentage)"
7.  Bug fix in all functions returning 1D or 2D intervals: in rare random cases NULL or invalid intervals set is returned
8.  Bug fix in gcompute_strands_autocorr: internal buffer overflow and possible memory corruption
9.  gdb.reload: run-time improvements

# misha 3.7.1

1.  Bug fix in gintervals.liftover and gtrack.liftover: some intervals might fail to be translated
2.  Bug fix in gintervals.liftover: "object 'f' not found" error if chain intervals are used in "chain" parameter

# misha 3.7.0

1.  Ubuntu support (multitasking mode still not thorougly tested)
2.  Bug fix: no progress report in multitasking mode
3.  Bug fix in gdb.create: temporary directory is created under the current GROOT instead of the new one

# misha 3.6.0

1.  Bug fix: occasional defunc processes AND/OR hanging in multitasking mode

# misha 3.5.6

1.  Run-time optimizations

# misha 3.5.5

1.  Avoid call to gdb.reload() (slow on large DB) in various functions: that create or remove tracks or intervals sets
2.  Bug fix in gintervals.neighbors: with 2D intervals the number of the returned neighbors might be less than "maxneighbors" parameter
3.  Bug fix in gintervals.neighbors: with 2D intervals NULL might be returned instead of NA if na.if.notfound=T

# misha 3.5.4

1.  Improved control over total maximal memory consumption in multitasking mode via gmax.mem.usage option.

# misha 3.5.3

1.  Run time optimizations when using several virtual tracks based on the same array track, differing only by slice

# misha 3.5.2

1.  Allow usage of sparse / arrays tracks in place of intervals
2.  Allow usage of big intervals sets in gintervals.diff, gintervals.intersect, gintervals.mapply, gintervals.union
3.  Run time optimizations when 1D big intervals set is used for scope
4.  Run time optimizations in various gintervals.* functions when big intervals sets are used
5.  Bug fix: functions might get stuck or crash when array track / sparse track / 1D big intervals set iterator is used along with 1D big intervals scope
6.  Bug fix: functions might get stuck or crash when 2D big intervals set iterator is used along with 2D big intervals scope
7.  Bug fix in gintervals.neighbors and gintervals.intersect: result might be poorly sorted when big interval sets are used
8.  Bug fix in gintervals.load: on an empty set or chrom returns NULL for small intervals sets and an empty data frame for a big intervals set
9.  Bug fix: "100%..." or "100%" is sometimes printed as the only progress report
10.  Bug fix: multiple progress report in some functions

# misha 3.5.1

1.  Bug fix in various gintervals.* functions: invalid output (except for the first row) when 2D track is used for intervals
2.  Bug fix in gextract: incorrect intervalID returned when 2D track is used for intervals

# misha 3.5.0

1.  Allow usage of 2D track in place of intervals
2.  Add progress report to gintervals.load
3.  Bug fix when using 2D big intervals: in some cases some or all chromosomes of the big intervals set might be skipped when big intervals set is used as a scope
4.  Big intervals set: before load verify that the size of a single chromosome (or chromosome pair) does not exceed gmax.data.size
5.  Bug fix in gcluster.run: clean up of running processes might not be completed if Ctr+C is pressed multiple times
6.  Bug fix in gintervals.load: returns all 2D intervals instead of a subset if one of chrom1/chrom2 is NULL and another one is not NULL
7.  Bug fix in gintervals.load: invalid row names if chrom / chrom1 / chrom2 parameter is used for a small intervals set

# misha 3.4.3

1.  Support intervals represented by tibbles

# misha 3.4.2

1.  New function: gsample; returns N random samples from the specified track expression
2.  Improved random seed when options(grnd.seed=0): so far two calls occurring within a second used identical random generators

# misha 3.4.1

1.  Fixed compilation errors on some platforms
2.  Run time improvement in gintervals.intersect when big intervals sets are used
3.  Bug fix in gintervals.neighbors: returns NULL if 2D big intervals sets are used
4.  Bug fix: a few point tracks in a track expression might be used without specifying an iterator

# misha 3.4.0

1.  Dynamically limit memory use in multitasking mode
2.  Bug fix: race condition and potential crash in multitasking mode when one of the child processes exits shortly after it is launched

# misha 3.3.18

1.  Bug fix in gintervals.force_range: error when intervals are out of range

# misha 3.3.17

1.  Run-time optimizations in track expression evaluation
2.  Run-time optimizations in gintervals.neighbors

# misha 3.3.16

1.  Run-time optimizations when working with large data frames of intervals in: gintervals.chrom_sizes, gintervals.force_range, gintervals.save
2.  Run-time optimizations when working with large data frames of intervals and using intervals.set.out parameter in all the functions that accept this parameter

# misha 3.3.15

1.  Run-time optimizations when working with big intervals sets in: gintervals.load, gintervals.diff, gintervals.force_range, gintervals.intersect, gintervals.mapply, gintervals.neighbors, gintervals.rbind, gintervals.update, gintervals.union
2.  Bug fix in gintervals.diff, gintervals.intersect, gintervals.neighbors, gintervals.rbind, gintervals.union: "object 'intervals' not found" error in some cases when big intervals sets are used
3.  Bug fix in gintervals.rbind: result does not preserve the original order if big intervals are used
4.  New undocumented function: .grbind

# misha 3.3.14

1.  gintervals.neighbors: run-time optimizations (the answer is entirely generated in C++)

# misha 3.3.13

1.  gintervals.neighbors: sort the output by original ids of intervals1, then |distance| (Manhattan distance for 2D), then ids of intervals2

# misha 3.3.12

1.  New version of gintervals.neighbors replaces both the old gintervals.neighbors and gintervals.annotate. By default gintervals.neighbors returns the closest neighbor
2.  gintervals.neighbors: support 2D intervals
3.  New parameters in gintervals.neighbors: maxneighbors, mindist1, maxdist1, mindist2, maxdist2, na.if.notfound
4.  gintervals.neighbors: columns renamed in the output
5.  Bug fix in gintervals.neighbors: in certain 1D cases some neighbors are not stated
6.  Check interrupts (Ctr+C) in gintervals.neighbors
7.  gintervals.annotate is removed

# misha 3.3.11

1.  gintervals.annotate: change the output format (instead of annotation id fully attach annotation interval)
2.  gintevals.annotate: support 2D intervals
3.  Bug fix: some jobs might return NA when run with gcluster.run

# misha 3.3.10

1.  gtrack.2d.import_contacts: support contacts files in interval-value format
2.  gtrack.2d.import_contacts: reduce the number of simultaneously opened files
3.  gtrack.2d.import: print out coordinates of duplicated object

# misha 3.3.9

1.  New function: gintervals.2d.import
2.  New option: gbig.intervals.size - controls the threshold when big intervals sets are created. Default value: 1000000
3.  Bug fix in gintervals.2d.import_contacts: utterly huge tracks might have missing areas of contacts
4.  Reduced the default value of gmax.processes option from 64 to 16
5.  Bug fix: incorrect progress report in gtrack.2d.import_contacts

# misha 3.3.8

1.  New function: gintervals.rbind. Runs rbind on intervals sets including big intervals sets on disk
2.  New "intervals.set.out" parameter added to: gtrack.array.extract, gwilcox
3.  Support big intervals sets in: gseq.extract, gtrack.array.extract, gtrack.modify
4.  Bug fix: gtrack.2d.import_contacts does not recognize chromosomes that have "chr" prefix

# misha 3.3.7

1.  Bug fix for all functions using 2D iterators: in some cases full chromosome pairs can be skipped. Bug first appeared in 3.3.0

# misha 3.3.6

1.  Bug fix in gdb.create: "Error in .gintervals.check_new_set(intervals.set.out)..."

# misha 3.3.5

1.  Added 'opt.flags' parameter to gcluster.run. Use this parameter to add restrictions to the machines that run submitted jobs: minimal RAM requirement, explicit hostnames list, etc. See man for qsub, "-l" flag.
2.  Support big intervals sets in the following functions: gpartition
3.  New "intervals.set.out" parameter added to: gpartition, gsegment
4.  Bug fix: invalid error recovery in glookup - traces from intervals.set.out might be left

# misha 3.3.4

1.  Interface change: gintervals.neighbors returns a data frame containing full intervals instead of their ids
2.  gintervals.neighbors: colnames parameter removed
3.  Support SAM files in gtrack.import_mappedseq
4.  Support big intervals sets in the following functions: gintervals.neighbors, glookup
5.  New "intervals.set.out" parameter added to: gintervals.neighbors, glookup
6.  Removed gintervals.merge function
7.  Runtime optimizations when big intervals sets are used in the following functions: gintervals.annotate, gintervals.diff, gintervals.force_range, gintervals.intersect, gintervals.mapply, gintervals.save, gintervals.union

# misha 3.3.3

1.  Bug fix in gquantiles, multitasking version (which is used by default): invalid quantiles might be returned if the number of iterator intervals exceeds gmax.data.size / number_of_child_processes. number_of_child_processes equals at most to the number of different chromosomes (or chromosome pairs for 2D) used in in "intervals" parameter
2.  Support big intervals sets in the following functions: gintervals.mapply
3.  New "intervals.set.out" parameter added to: gintervals.mapply

# misha 3.3.2

1.  Support big intervals sets in the following functions: gintervals.annotate, gintervals.diff, gintervals.intersect, gintervals.union, gsegment, gwilcox
2.  New "intervals.set.out" parameter added to: gintervals.annotate, gintervals.diff, gintervals.intersect, gintervals.union
3.  Added progress report to: gintervals.save, gintervals.force_range

# misha 3.3.1

1.  Support big intervals sets in the following functions: gcis_decay (only _intervals_ parameter), gintervals.2d.band_intersect
2.  New "intervals.set.out" parameter added to: gintervals.2d.band_intersect
3.  Bug fix in gintervals.force_range: "Error in if (size > max.data.size) { : argument is of length zero"

# misha 3.3.0

1.  New concept: big intervals sets
2.  Big intervals sets can be used for _iterator_ parameter in all functions
3.  Big intervals sets can be used in _intervals_ parameter in the following functions: gdist, gextract, gquantiles, gscreen, gsummary, gbins.quantiles, gbins.summary, gintervals.quantiles, gintervals.summary, giterator.intervals
4.  Interface change: gintervals.quantiles, gintervals.summary now return also the source intervals
5.  New functions: gintervals.is.bigset, gintervals.chrom_sizes, gintervals.update
6.  New "chrom", "chrom1", "chrom2", parameters for gintervals.load
7.  New "intervals.set.out" parameter added to: gextract, gscreen, gintervals.force_range, gintervals.quantiles, gintervals.summary, giterator.intervals
8.  Restrict gintervals.quantiles and gintervals.mapply to work with only 1D and Fixed Rectangle iterators
9.  Changed the position of "file" parameter in gextract
10.  Bug fix: gscreen on vtrack with global.percentile.max returns different number of intervals in each run
11.  Bug fix: in 2D iterators progress report can sometimes go backwards
12.  Bug fix: empty intervals set is ignored if used as an iterator
13.  Bug fix: crash if an empty intervals set is used for scope
14.  Bug fix: memory leak in giterator.cartesian_grid

# misha 3.2.8

1.  Bug fix in gintervals.save: "variable shadows the name of the intervals set" error can be generated even if auto-completion is turned off
2.  Bug fix in all functions that create a new track: "variable shadows the name of the new track" error can be generated even if auto-completion is turned off
3.  Bug fix in gdb.create, gtrack.var.set and while creating P-values table: insufficient permissions might be given to created directories

# misha 3.2.7

1.  Bug fix: no proper clean up if error is generated while P-values table is loaded

# misha 3.2.6

1.  New function: gcis_decay

# misha 3.2.5

1.  New format for tracks created by gtrack.2d.import. This format uses on average 30% less space. Old format can still be used
2.  Old computed tracks now require conversion
3.  Bug fix: crash reading computed tracks

# misha 3.2.4

1.  Changed format of 2D tracks (rectangles and computed). Use gtrack.convert to convert the old tracks
2.  Added undocumented function: .gdb.convert_tracks
3.  Support track files larger than 2 Gb on 32-bit platforms
4.  gtrack.2d.import_contacts: support creation of huge 2D tracks (practically unlimited size) whilst constant memory usage
5.  gtrack.2d.import_contacts: allow arbitrary order of contacts within contacts file or between several contacts files
6.  gtrack.2d.import_contacts: improved progress report
7.  Ensure binary consistency of 2D track files. Previously the files representing two identical 2D tracks could differ on binary level

# misha 3.2.3

1.  gtrack.2d.import_contacts: sum up duplicated contacts if 'allow.duplicates' is TRUE
2.  gtrack.2d.import_contacts: allow contacts to be passed in multiple files

# misha 3.2.2

1.  Virtual tracks are not created anymore as variables in R environment (dummy variables are created in autocompletion mode). Virtual tracks are stored inside GVTRACKS variable.
2.  Virtual tracks are not reset anymore on gsetroot or gdir.cd
3.  New function: gvtrack.info
4.  Removed functions: gvtrack.all.load, gvtrack.all.rm, gvtrack.all.save, gvtrack.import
5.  Bug fix: GERROR_EXPR variable reported but not set by functions supporting multitasking

# misha 3.2.1

1.  Bug fix in gcluster.run: on some systems a warning is generated "bash: module: line 1: syntax error: unexpected end of file\nbash: error importing function definition for `module'\n"

# misha 3.2.0

1.  Track variables are referenced by two parameters: track, varname instead of "track.varname"
2.  Renamed and adopted for new track variable convention:  
    gvar.load => gtrack.var.get  
    gvar.save => gtrack.var.set  
    gvar.ls => gtrack.var.ls  
    gvar.rm => gtrack.var.rm
3.  Remove gvar.exists function
4.  New concept: track attributes. Use .gdb.convert_attrs() to convert old trackdb to the new format.
5.  New functions: gtrack.attr.get, gtrack.attr.set, gtrack.attr.export, gtrack.attr.import, gdb.get_readonly_attrs, gdb.set_readonly_attrs
6.  created.by and created.date are no longer track variables but rather read-only track attributes
7.  gtrack.ls: allow filtering by track attributes
8.  New obligatory "description" parameter in gtrack.2d.create, gtrack.2d.import_contacts, gtrack.array.import, gtrack.convert, gtrack.create, gtrack.create_pwm_energy, gtrack.create_sparse, gtrack.import, gtrack.import_mappedseq, gtrack.import_set, gtrack.liftover, gtrack.lookup, gtrack.smooth
9.  New function: gset_input_mode. This function replaces gparam.type option and controls auto-completion of track / intervals names
10.  By default interactive mode is switched off (equivalent to gparam.type="string" in older version). Auto-completion is switched off as well.
11.  Check parameters correctness in giterator.cartesian_grid and not only when the iterator is actually used
12.  Bug fix in gcluster.run: distributed command resets GROOT and various package options
13.  Bug fix: in non-interactive input mode ("string" var mode) gvtrack.array.slice fails
14.  Bug fix: gsetroot, gdb.reload, gdb.cd leave traces in the environment if they stop on error
15.  Bug fix: gtrack.modify incorrectly sets created.by attribute
16.  Bug fix: invalid usage printed in gvtrack.all.load

# misha 3.1.11

1.  Changed the policy for multitasking job distribution

# misha 3.1.10

1.  Multitasking for glookup, gtrack.smooth, gintervals.quantiles
2.  Bug fix: memory corruption in gintervals.mapply when intervals==ALGENOME and multitasking is turned off

# misha 3.1.9

1.  Multitasking for gintervals.mapply
2.  Bug fix: invalid format of data frame returned by gintervals.mapply when intervals==ALGENOME
3.  Bug fix: error while preparing a track for percentiles queries

# misha 3.1.8

1.  Multitasking for gtrack.create

# misha 3.1.7

1.  Multitasking for gtrack.create_pwm_energy

# misha 3.1.6

1.  Multitasking for gquantiles

# misha 3.1.5

1.  Bug fix: some genomic intervals might be missing in gintervals.load_chain
2.  Bug fix: gintervals.liftover might incorrectly convert some genomic intervals to NULL
3.  Bug fix: gtrack.liftover might incorrectly set NA for some converted genomic intervals

# misha 3.1.4

1.  Bug fix: unreleased shared memory or/and named semaphore if R/misha crashes or is terminated with a signal
2.  Bug fix: in some cases 3 seconds delay in multitasked functions
3.  Bug fix: in some cases unresponsiveness on Ctrl+C in multitasked functions
4.  Bug fix: with non-default options gquantile could return incorrect value for extreme percentiles (close to 0 or to 1)

# misha 3.1.3

1.  Bug fix: deadlock in all multitasked functions (gsummary, gextract, gdist)

# misha 3.1.2

1.  Multitasking for gscreen, gsummary
2.  New R option for the package: gmax.processes, default: 64
3.  Bug fix: "child process ended unexpectedly" error in multitasked functions when evaluation of track expression fails
4.  Bug fix: Ctrl+C might stop working in R after evaluation of track expression fails

# misha 3.1.1

1.  Multitasking for gextract
2.  Bug fix: virtual tracks do not work in gdist
3.  Bug fix: potential crash and process table bloating in gdist
4.  Bug fix: potential freezing (deadlock) in gdist
5.  Bug fix: error "2D iterator is used along with 1D intervals" in gdist with 2D iterator and intervals==ALLGENOME

# misha 3.1.0

1.  Multitasking for gdist

# misha 3.0.4

1.  Bug fix: in gintervals.2d "Error in is.null(strands) : 'strands' is missing"

# misha 3.0.3

1.  Add strand parameter to gintervals

# misha 3.0.2

1.  New function: gtrack.import
2.  Allow tab-delimited files in gtrack.import_set
3.  gintervals.import_genes / gdb.init: add kgID column to annotations
4.  gintervals.import_genes / gdb.init: eliminate identical values in overlapping intervals' annotation
5.  Bug fix: in tab-delimited files if end coordinate equals the chrom size an error is reported

# misha 3.0.1

1.  Bug fix: gintervals.import_genes / gdb.create switches between utr3 and utr5

# misha 3.00

1.  New track type: array
2.  New functions: gtrack.array.import, gtrack.array.get_colnames, gtrack.array.set_colnames, gtrack.array.extract, gvtrack.array.slice
3.  Renamed:  
    gintervals.band.intersect => gintervals.2d.band_intersect  
    giterator.cartesian.grid => giterator.cartesian_grid  
    gsetroot.examples => gdb.init_examples  
    gtrack.create_2d => gtrack.2d.create  
    gtrack.import.2d_contacts => gtrack.2d.import_contacts  
    gtrack.import.mappedseq => gtrack.import_mappedseq  
    gtrack.import.wigs => gtrack.import_set
4.  gsetroot has a new alias: gdb.init
5.  gtrack.import_set: create a sparse track if binsize==0
6.  gextract: allow saving result in a tab-delimited file
7.  gintervals.force_range: eliminate intervals with non-existent chromosome
8.  gquantile / gintervals.quantile, quantile / global.percentile / global.percentile.min / global.percentile.max functions of a virtual track: use weighted average of nearest samples instead of picking up the closest sample
9.  Bug fix: sometimes using invalid value of quantile.edge.data.size option. Result: suboptimal precision at the edges for quantile calculation OR memory bloating for quantile calculations on large sets of data.
10.  Bug fix: sometimes using invalid value for gtrack.chunk.size option. Result: suboptimal performance for newly created 2D tracks + memory bloating while reading 2D tracks.
11.  Bug fix: sometimes using invalid value for gtrack.num.chunks option. Result: slow performance while reading 2D tracks OR memory bloating.
12.  Bug fix: gsetroot / gdb.reload / gdir.cd corrupts the database state if one of the names shadows a variable in R environment.

# misha 2.75

1.  giterator.cartesian.grid: replace 'expansion' parameter with 'expansion1' and 'expansion2' parameters for each axis
2.  giterator.cartesian.grid: changed the order of the parameters
3.  Bug fix: cartesian grid iterator incorrectly restricted the expansion between two neighbouring centers C1, C2 to be (C2-C1)/2

# misha 2.74

1.  New function: gdb.create
2.  gdbreload renamed to gdb.reload
3.  Added support of ftp and zipped files in gintervals.import_genes
4.  Documentation updates

# misha 2.73

1.  New function: gintervals.import_genes
2.  Documentation updates

# misha 2.72

1.  Added User Manual in PDF, Reference Manual in PDF and HTML.
2.  Updated functions documentation.

# misha 2.71

1.  Added documentation for each function from R
2.  Do not require libR.so for installation
3.  Added "maxread" parameter to gcompute_strands_autocorr()
4.  Do not unify overlapping intervals in gintervals
5.  gintervals.apply is replaced with gintervals.mappy. The function interface changes.
6.  Renamed: gvar.get() to gvar.load() and gvar.set() to gvar.save()
7.  Bug fix: gcluster.run() did not load the package
8.  Bug fix: gcluster.run() corrupted the return value
9.  Bug fix: track expression iterator might miss an interval if intervals are not canonic and the iterator type is intervals/sparse

# misha 2.70

1.  "misha" becomes an R package
2.  gversion() removed

# misha 2.60

1.  Bug fix: gdir.cd crashes

# misha 2.59

1.  Support 2D tracks (Rectangles type) in gtrack.liftover
2.  Support 2D intervals in gintervals.liftover
3.  Bug fix: gtrack.liftover does not remove temporary files

# misha 2.58

1.  Added gintervals.liftover function
2.  Bug fix: error while trying to access a 2D track.

# misha 2.57

1.  Added gtrack.liftover and gintervals.load_chain functions

# misha 2.56

1.  Cleaned up stdout field in the result of gcluster.run

# misha 2.55

1.  Bug fix: scope might be incorrectly applied while using cartesian iterator
2.  Bug fix: gtrack.import.2d_contacts crashes when fend is out of range

# misha 2.54

1.  New functions: gbins.quantiles and gbins.summary
2.  giterator.cartesian.grid: do not implicitly add zero expansion
3.  Add R parameter in gcluster.run
4.  Add support of bedGraph extention in gtrack.import.wigs
5.  Bug fix: overlapping 2D intervals are not reported correctly

# misha 2.53

1.  New function: gcluster.run

# misha 2.52

1.  Faster gsetroot using cached list of tracks and intervals.
2.  New rescan parameter for gsetroot and gdbreload.

# misha 2.51

1.  Added 2D intervals support in gintervals.apply.
2.  Bug fix: gtrack.import.mappedseq crash.

# misha 2.50

1.  Added band parameter to gdist, gextract, glookup, gpartition, gquantiles, gscreen, gsummary, gtrack.create, gtrack.lookup, gintervals.quantiles, gintervals.summary, giterator.intervals.
2.  New function: gintervals.band.intersect
3.  Removed min.band, max.band parameters from giterator.cartesian.grid

# misha 2.43

1.  Run-time optimizations for 2D queries

# misha 2.42

1.  Bug fix: crash if "preparing track for percentiles queries" is interrupted with CTRL-C

# misha 2.41

1.  New function: gdbreload
2.  gtrack.import.wigs: create tmp directory in GROOT/downloads rather than in /tmp
3.  gwget: by default use path == GROOT/downloads
4.  Bug fix in gtrack.import.wigs: do not proceed to import if one of the previous steps (ftp/unzip/convert to wig) failed or interrupted

# misha 2.40

1.  Support BigWig / BedGraph formats in gtrack.import.wigs
2.  Bug fixes in gtrack.import.wigs

# misha 2.39

1.  New functions: gwget, gtrack.import.wigs
2.  Allow to use unsorted intervals for gtrack.modify, gtrack.create_sparse
3.  Allow to use unsorted and overlapping intervals for gintervals.intersect, gintervals.union and gintervals.diff

# misha 2.38

1.  Run-time optimizations for 2D cartesian grid iterator
2.  Bug fix: progress report goes reports 100% before the actual completion of command

# misha 2.37

1.  Bug fix: 2D cartesian grid iterator skips some of the iterator intervals

# misha 2.36

1.  Bug fix: gtrack.convert does not correctly convert old (version 1) computed tracks

# misha 2.35

1.  Fix memory leaks when command exists on error or is interrupted by CTRL-C

# misha 2.34

1.  Added "sum" virtual track function
2.  Added "quantile" virtual track function
3.  Renamed "percentile*" virtual track functions to "global.percentile*"
4.  Bug fix: when iterator interval does not intersect any intervals of sparse track "stddev" virtual track function returns last value instead of NaN

# misha 2.33

1.  Added "stddev" virtual track function

# misha 2.32

1.  Bug fix: memory corruption when the track expression contains more than one virtual track of "distance" or "distance.center" type.
2.  Bug fix: incorrect statistics in gsummary / gintervals.summary when the summary is done on 1 sample.

# misha 2.31

1.  Changed the format of computed 2D tracks
2.  Bug fixes in gtrack.convert
3.  Bug fixes in GenomeTrack::get_type
4.  New function: .gtrack.create_test_computer2d
5.  gintervals.2d.force_range removed (gintervals.force_range now works for both 1D and 2D intervals)

# misha 2.30

1.  Changed the format of 2D tracks, now using StatQuadTreeCached class
2.  Added gtrack.convert function
3.  gwrite.table, gread.table were removed

# misha 2.26

1.  Iterative algorithm in gtrack.smooth is reset once in a while to prevent loss of precision in floating point calculations
2.  Bug fix: giterator.cartesian.grid does not work correctly if min/max.band are NULL.

# misha 2.25

1.  Added gintervals.force_range and gintervals.2d.force_range

# misha 2.24

1.  Change the default of min/max.band in cartesian iterator to NULL
2.  Change the default of intervals2 in cartesian iterator to NULL
3.  Add min/max.band.idx to cartesian iterator

# misha 2.23

1.  Changed virtual track function "dist" to distance

# misha 2.22

1.  Added gtraceback

# misha 2.21

1.  Added 2D computed tracks

# misha 2.20

1.  Changed the format of 2D tracks
2.  Added "sum" and "area" functions to virtual tracks
3.  Bug fix: Error message "Cannot implicitly determine iterator policy" when used with two or more virtual tracks pointing to the same physical track

# misha 2.18

1.  Bug fix: gtrackimport_mappedseq, gcompute_strands_autocorr might skip the last portion of the input file
2.  Bug fix: gtrackimport_mappedseq, gcompute_strands_autocorr might skip the first row of the input file

# misha 2.17

1.  Added fixed rectangle iterator
2.  Added support of 2D in gpartition
3.  Added support of 2D in gtrack.lookup
4.  Added support of 2D in gintervals.canonic
5.  Added support of 2D in gintervals.intersect
6.  Added gvtrack.ls
7.  Check user interrupt in gseq.extract
8.  Removed: gtrack.cache
9.  Bug fix: incorrect error message in gwilcox when used with non fixed-bin iterator
10.  Bug fix: iterator=trackname produces an error
11.  Bug fix: gtrack.lookup produces an error
12.  Bug fix: gtrack.modify produces an error
13.  Bug fix: gvtrack.import incorrectly imports a track if called within a function
14.  Bug fix: gtrack.import.2d_contacts memory leak

# misha 2.16

1.  Add band control to cartesian grid iterator
2.  Bug fix: cartesian grid iterator might produce incorrect results while using scope

# misha 2.15

1.  Added new vtrack function: dist.center
2.  Bug fix: precision loss in gextract due to float / double conversion

# misha 2.14

1.  Bug fix: gintervals.* might fail on overlapping intervals

# misha 2.13

1.  Run-time optimizations for gintervals and gintervals.2d
2.  Run-time optimizations when using GITERATOR.INTERVALS

# misha 2.12

1.  Added cartesian grid iterator

# misha 2.11

1.  Remove canonic / original property for 2D tracks

# misha 2.10

1.  Added virtual tracks

# misha 2.02

1.  Hide the chrom1/chrom2 swap in 2D intervals from the user

# misha 2.01

1.  Added support for 2D intervals in gtrack.create
2.  Added gtrack.info
3.  Removed gtrack.binsize
4.  Added gintervals.all
5.  Added gintervals.2d.all

# misha 2.00

1.  Added gintervals.2d
2.  Added gtrack.create_2d
3.  Added support for 2D intervals in gscreen, gextract, glookup, gsummary, gsummary.intervals, gquantiles, gintervals.quantiles, gdist, giterator.intervals
4.  Add unify_touching_intervals parameter to gintervals.canonic()

# misha 1.25

1.  Automatically build PV-table
2.  Remove gtrack.makepvals()

# misha 1.24

1.  Allow maximal precision of gquantiles / gintervals.quantiles near extreme probs (0 and 1)

# misha 1.23

1.  Allow gquantiles / gintervals.quantiles work on the whole genome (use random sampling)
2.  Allow using non-canonic intervals in gseq.extract()

# misha 1.22

1.  Added ".nearest" track function
2.  Bug fix: invalid error report when fixed-bin track size does not match chrom size

# misha 1.21

1.  Make .greloaddb faster (affect gsetroot, gdir.cd, ...)
2.  Do not autocomplete tracks variables
3.  Add gvar.ls()
4.  Remove gvar.load()
5.  Added pattern matching for gtrack.ls(), gintervals.ls() and gvar.ls()
6.  Added gtrack.exists() and gintervals.exists()
7.  Print the track expression in various error messages

# misha 1.20

1.  Make .greloaddb faster (affect gsetroot, gdir.cd, ...)
2.  Bug fix: in various functions: "Error, undefined gparam.type"

# misha 1.19

1.  Require all track directories to have an extention .track
2.  Forbid creation of tracks, intervals or directories inside track directories
3.  Forbid deletion of directories inside track directories
4.  Bug fix: gtrack.cache() sometimes creates tracks shorter by one sample than the chrom size

# misha 1.18

1.  Treat +/-Inf value in a track as NA
2.  Allow "dir" argument in gsetroot

# misha 1.17

1.  Add gtrack.import.mappedseq
2.  Add gcompute_strands_autocorr
3.  Support integer breaks in gdist()
4.  Suppress error report if force==TRUE and track/interv/dir do not exist in gtrack.rm / ginterv.rm / gdir.rm

# misha 1.16

1.  Maintain virtual working directory in gsetroot / gdir.cd. Do not change the shell working directory of the user.
2.  Add gdir.cwd

# misha 1.15

1.  Add gdir.create, gdir.rm, gdir.cd
2.  Remove gtrackset.*, gintervset.*

# misha 1.14

1.  Stop creating variables with tracksets/intervsets names for TAB-completion
2.  Before addition / removal of track/intervals variables check whether they already exist
3.  Bug: gtrack.cache overrides "created.by" and "created.date" variables

# misha 1.13

1.  Rename intervals files by adding them ".interv" extention

# misha 1.12

1.  Remove global/user attribute for tracksets/intervsets

# misha 1.11

1.  Add "iterator" parameter to the relevant functions
2.  Remove giterator.policy() function
3.  Remove gapply() function
4.  Bug in gtrack.create_pwm_energy(): Error in sprintf("gtrack.create_pwm_energy(%s, %g, %s, %s, %g, track=%s)", : invalid format '%g'; use format %s for character objects
5.  Add "breaks" attribute to gdist result

# misha 1.10

1.  Support sparse tracks and track expression iterators over intervals
2.  Added gversion()
3.  Check in giterator.intervals that the memory is not blown up
4.  Do not unify touching intervals in giterator.policy()
5.  error in gtrack.makepvals(): "Expression does not produce a numeric result"
6.  Allow GAPPLY.INTERVALS and GITERATOR.INTERVALS variables

# misha 1.02

1.  Stop supporting point intervals (start coordinate == end coordinate). Automatically convert old point intervals to (chrom, start, start+1).
2.  Bug fix: gtrack.modify() fails with ".created.by.intervs" error.

# misha 1.01

1.  Replace global variables with options: GMAX_DATA_SIZE => gmax.data.size, GBINSIZE => gstep, GBUFSIZE => gbuf.size
2.  Add an option "gparam.type" controlling the type of the arguments ("var" or "string")
3.  gintervals.annotate produces incorrect results when annotation intervals contain overlapping intervals

# misha 1.00

1.  Fix an error when the expression is too long or something like that (Rami for details) Fixed by Amos 5/7/10 patched .Rprofile
2.  Fix gcreate crash
3.  If I have track names with a minus sign in them - they are built correctly - but I can't do conditions on them. e.g. : rv <- gscreen (GSE14097.ERG_AR_re_ChIP-7 > 20 && TE.WT__Input < 7). Resolve but blocking - signs in track names (and maybe other sensitive characters). Fixed: track names are not allowed to contain characters other than alphanumeric and _. The fix applies for gcreatetrack, gcachemultires, gsmooth.
4.  added support for writing multi resolution track binning (gcachemultires) - for use by tgbr.
5.  Misha: grmtrack does not remove the dataset directory when the last track is deleted
6.  Misha: add force delete command line option for grmtrack for track deletion without confirmation
7.  Misha: add grmdataset function to delete the whole dataset recursively
8.  Misha: eliminate use of the slow gsetroot() for track list refresh after track creation or removal The fix applies to grmtrack, gcreatetrack, gcachemultires, gsmooth.
9.  Rami: On failure to create a track - remove it. The fix applies for gcreatetrack, gcachemultires, gsmooth.
10.  Misha: allow a dot in the name of the track variable.
11.  should store the command that created a track in a track varible (For future reference), we should also save the date (in the future we may want to do something with dependencies). Each newly created track is assigned two track variables: created.by and created.date.
12.  Rami: Easier dump to text of info (current need to use sep="\t" rownames = F, etc) Added gwrite.table() function.
13.  Summary statistics on tracks (a simplification of the distribution feature). We want to compute the total, average, stdev, min and max of a track, number of bins, number of NAN bins. All that without going throuhg the distribution computation (maybe using R function that seat over the gdist function) New function added: gsummary.
14.  Change coordinate->bin function (should be int(coordinate/bin_size)
15.  Iterator protect/unprotect (eval_next_int, eval_next_bool etc)
16.  Misha: gsetroot() should undefine previous track variables
17.  Misha: mix of float and double in gdist() causes incorrect bin assignment for the values at the border.
18.  Misha: add smart progress report.
19.  Misha: end coordinate of the last segment of the chromosome was mistakenly extended by binsize.
20.  Misha: limit the range that intervals cover for gextract and gquantiles to prevent memory allocation failures.
21.  Misha: enhance gintervals() function to accept chroms as strings without "chr" prefix or as integers.
22.  eliminate gpath support.
23.  Allow various services (gscan, gdist) to be limited to a subset of the genome (defined by a give interval set). This will probably require changing the classes that iterate over chromosome such that all will use a common GenomeTrack interface.
24.  Misha: allow track expressions for gextract and gquantiles.
25.  Rami: Find out why you need to do 'as.vector' on output of gquantiles to get regular numbers. Misha: as designed. gquantiles returns a table, not a vector. The table is of NxM size where N is the number of intervals, and M is the number of quantiles.
26.  Return intervals as a dataframe from C code rather than a list that should be later converted to a dataframe. Misha: affects gscreen, gwilcox, gintervals.union, gsetroot. Improves run time for large sets of intervals.
27.  Misha: intervals returned by gwilcox should cover the whole area covered by the small window rather than just the center of the window.
28.  Why is the usage of: PeakIntervals_TE_AR_AR <- gwilcox (TE.TE_AR_DHT__AR, 100000, 1000, maxpval = 0.001) - giving me intervals with peak values of 0? Misha: added _what2find_ option to gwilcox. This option controls whether peaks/lows or both should be searched by gwilcox.
29.  Implemenet operations on intervals. Union: generate a new interval set, which include an interval for each closure of the two given set. Intersection: returns an interval set with all the non empty intersections of intervals in the two sets. Difference: returns the intervals of set A such that intersection with interval from B are removed. Misha: added gintervals.union(), gintervals.intersect(), gintervals.diff() function
30.  Call R eval in bulks to boost run times. Misha: added ".bufsize=1" parameter to all functions that accept track expressions. Increase "bufsize" to boost the performance.
31.  Misha: allow command interruption by "ctrl-C" in various function such as gwilcox. (Check R event loop?)
32.  Misha: if "ctrl-C" is pressed whilst gcreate/gsmooth/gcachemult track directory is not removed.
33.  Add gintervals.annotate function
34.  Add auto-completion with TAB for annotations.
35.  Add gls.annots() function
36.  Add support for annotation function "annotation.dist" in the track expressions
37.  Misha: determine .bufsize parameter automatically based on the dimensions of the first evaluation + the interval size Misha: for now interval size is not checked
38.  Misha: fix error messages format Misha: fixed with an ugly hack
39.  Add gapply
40.  Allow TrackScanner to iterate over a few track expressions
41.  Bug: message "Error: GBINSIZE variable is undefined or not numeric" is printed when invalid trackname is used in track expression.
42.  Bug: invalid results whenever the same track is used in more than one track expression in gdist()
43.  Bug: incorrect pvals when the limits are falling on the maximum
44.  Bug: incorrect annot.dist when you also use intervals=something
45.  Bug: gscreen returns one bin less for each interval
46.  create all track dirs with group open permissions
47.  Provide quick gsetroot (for scripts - e.g. not set up of completion variables? focus on subset of dirs? lazy init of gsetroot?) Misha: the performance of the original gsetroot was optimized
48.  Take care of environment handling in various eval(parse()) constructs
49.  Allow string parameters for track_expr - if the parameter is a string, just avoid parse
50.  Reverse the convention of positive / negative distance between coord / interval and annotation
51.  Change DB directory structure: tracks/tracksets/tracknames, annots/annotsets/annotnames Misha: renamed functions: gls => gls.tracks, gcreatetrack => gcreate.track, grmtrack => grm.track, grmdataset => grm.trackset. Newly added functions: gls.annots, gls.annotsets, gls.tracksets, gcreate.annotset, grm.annotset.
52.  Add trackset / annotset manipulation functions + user / global flag support
53.  Bug: if ~ symbol is used in gsetroot, track files cannot be accessed
54.  Allow custom column names for gapply and gextract
55.  Regression tests
56.  make it easier to import intervals/annotations: add ginterval.import that will take non cannonical intervals and will sort them and unify them to become disjoint. The function will return the cannonical ginterval object plus a factor mapping ids in the original data frame to the new cannoincal gintervals. We can then use tapply to import meta-data.
57.  Merge annotations and intervals concepts. Misha: renamed functions: gls.annots => gls.intervs, gls.annotset => gls.intervsets, grm.annots => grm.intervs, grm.annotsets => grm.intervsets. Newly added functions: gintervals.load, gintervals.save.
58.  Effi: segmentation fault in gdist
59.  Buf fix: calling gprepare_pvals() twice on the same track fails
60.  need a simple way to merge data onto intervals, without reordering and invalidating the intervals (i.e. doing ginterval.import after merge is annoying. This may be solvable using standard r merge options, but we need to wrap it up nicely) Misha: Added gintervals.merge function
61.  More protective use of intervals: coerce fields as factors if needed and generally be aware of potential user changes to intervals.
62.  allow colelction of pairs of interval within a maximal distance Misha: Added gintervals.neighbors function
63.  Bug fix: gcachemultires does not work
64.  Decide about the naming policy for the functions
65.  Restrict gdist to accept minval=2
66.  Add dist.XXX function
67.  check how many times gapply calls func1 in gapply(func1(x1), track, gintervals(1, 0, 100))
68.  gintervals.apply should reverse intervals of the minus strand
69.  bug: gintervals.neighbors(gintervals(1, 0, 2000, 1), gintervals(1, 500, 1500, 1), 0, 0) produces an error
70.  added colnames parameter to gintervals.neighbors
71.  file descriptors are not closed properly in onexit() causing various *.create* functions to leave junk if interrupted
72.  chromo sequence and pwms to the trackdb directory: allow import of pwms to the library (using several sets of pwm, each defined by a file?)
73.  allow PWM energy computation (using the chromo sequence and a pwm)
74.  Create tracks from 4c sites and prof files
75.  BUG: The same interval.dist could not appear more than once in an expression - for example , gdist(global.tss.dist+global.tss.dist,-1000,1000,100) produces "object 'global.tss.dist' not found"
76.  Provide the interval ID in the function provided to intervals.gapply Misha: GINTERVID variable is maintained while gapply
77.  add progress report for non-track expression functions Misha: added progress report to gintervals.annotate
78.  Allow to smooth NaN values in gtrack.smooth (add smooth_nans parameter)
79.  it might be usefull to have a direct way get sets of intervals corresponding to division of the genome according to the values of some track. it would be best if gdist could return another extra value with intervals corresponding to each of the track combinations Misha: gpartition function was added
80.  Misha: avoid memory blow up when large vectors are returned by gquantile, etc.
81.  gseq.extract() should return reverse-complementary sequence for strand -1
82.  Bug fix: gpartition does not treat NaNs correctly
83.  Change gdist and gpartition to accept breaks rather than minval/maxval/numbins
84.  Bug: with equal bin size gdist crashes
85.  Bug: gquantiles crashes with error: unprotect_ptr: pointer not found
86.  Bug: gdist called from a function does not work Misha: deal with unevaluated (=promised) values
87.  In case of track expression result mismatch save the result of the last evaluation in GERROR_EXPR variable
88.  Bug: gtrack.make_pvals() does not work
89.  Bug: gsummary() max value is incorrect for negative values
90.  Added gvar.get() function
91.  Added gvar.exists() function
92.  gvar.rm() "warning: variable was not found" when variable exists but not loaded
93.  Added gtrack.binsize()
94.  Allow gintervals(chrom) which is equal to gintervals(chrom, 0, -1)
95.  Effi: if x contains a line of NA's gintervals.import(x) causes segmentation fault
96.  Misha: optimize gextract to present chroms as factors rather than strings
97.  Misha: truncate long column names
98.  Added glookup() and gtrack.lookup() functions
99.  Renamed "origin" attribute of gintervals.import to "mapping"
100.  Allow regression tests to be invoked by the matching string
101.  Removed gvar.loadall() function
102.  Added gtrack.modify() function
103.  Misha: delete gvar.loadall() - as it might blow up the memory in large databases
104.  Misha: gvar.* functions do not work when track variable X is loaded and passed unquoted
105.  Misha: read-only functions refuse to work with a track that does not have write permissions
106.  Allow overlaps all functions except gintervals.intersect, gintervals.diff, gintervals.union, dist.XXX
107.  Renamed gquantiles to gintervals.quantiles
108.  Renamed gintervals.import to gintervals.canonic
109.  Added gintervals.summary()
110.  Added gquantiles()
111.  Bug: gvar.load does not work when the variable is given unquoted
112.  Warning "is.na() applied to non-(list or vector) of type 'NULL'" when gscreen(is.na(track)) is called
113.  Versioning and installation added
114.  Calling gintervals.canonic(i) when i is a dataframe without any rows causes R to crash Misha: other functions (gintervals.union, ...) were fixed too