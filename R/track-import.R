# Track import functions and helpers

get_bigWigToWig_bin <- function() {
    dir <- tempdir()
    utils::untar(system.file("bigWigToWig.tar.gz", package = "misha"), exdir = dir)
    return(file.path(dir, "bigWigToWig"))
}


# Internal helpers for gtrack.import
.gtrack_is_bed_path <- function(path) {
    grepl("\\.bed(\\.(gz|zip))?$", path, ignore.case = TRUE, perl = TRUE)
}

.gtrack_read_bed <- function(file) {
    # Try to use data.table::fread if available, with fallback to read.table
    # Pre-filter header lines with grep (if available) before reading
    bed <- NULL

    # Try data.table::fread (fastest option)
    if (requireNamespace("data.table", quietly = TRUE)) {
        tryCatch(
            {
                # Check if grep is available for pre-filtering
                can_grep <- nzchar(Sys.which("grep"))
                if (can_grep) {
                    # Use grep to filter out header lines BEFORE reading
                    cmd <- sprintf("grep -vE '^(track|browser|#|$)' %s", shQuote(file))
                    bed <- data.table::fread(
                        cmd = cmd,
                        header = FALSE, sep = "\t", quote = "",
                        fill = TRUE, stringsAsFactors = FALSE,
                        data.table = FALSE, showProgress = FALSE
                    )
                } else {
                    # No grep: read then filter in R
                    bed <- data.table::fread(
                        file,
                        header = FALSE, sep = "\t", quote = "",
                        fill = TRUE, stringsAsFactors = FALSE,
                        data.table = FALSE, showProgress = FALSE
                    )
                    # Filter header lines in R if we got data
                    if (!is.null(bed) && nrow(bed) > 0) {
                        v1 <- trimws(bed[[1]])
                        keep <- !(startsWith(v1, "track") | startsWith(v1, "browser") |
                            startsWith(v1, "#") | v1 == "")
                        bed <- bed[keep, , drop = FALSE]
                    }
                }

                # Validate: ensure we got reasonable data
                if (is.null(bed) || nrow(bed) == 0 || ncol(bed) < 3) {
                    bed <- NULL
                }
            },
            error = function(e) {
                bed <<- NULL
            }
        )
    }

    # Fall back to utils::read.table if fread failed
    if (is.null(bed)) {
        bed <- utils::read.table(
            file,
            header = FALSE, sep = "", quote = "", comment.char = "",
            fill = TRUE, stringsAsFactors = FALSE, colClasses = "character"
        )

        # Filter header lines
        if (nrow(bed) > 0) {
            v1 <- trimws(bed[[1]])
            keep <- !(startsWith(v1, "track") | startsWith(v1, "browser") |
                startsWith(v1, "#") | v1 == "")
            bed <- bed[keep, , drop = FALSE]
        }
    }

    # Validate the final result
    if (is.null(bed) || nrow(bed) == 0) {
        stop(sprintf("BED file %s appears to be empty or contains no data intervals", file), call. = FALSE)
    }

    if (ncol(bed) < 3) {
        stop(sprintf("BED file %s appears to be malformed (less than 3 columns)", file), call. = FALSE)
    }

    chrom <- as.character(bed[[1]])
    start <- suppressWarnings(as.numeric(bed[[2]]))
    end <- suppressWarnings(as.numeric(bed[[3]]))

    if (any(is.na(start)) || any(is.na(end))) {
        stop(sprintf("Non-numeric coordinates detected in BED file %s", file), call. = FALSE)
    }

    values <- rep(1, nrow(bed))
    if (ncol(bed) >= 5) {
        score <- suppressWarnings(as.numeric(bed[[5]]))
        if (!all(is.na(score))) {
            values <- ifelse(is.na(score), 1, score)
        }
    }

    intervals_df <- data.frame(
        chrom = chrom,
        start = start,
        end = end,
        stringsAsFactors = FALSE
    )

    list(intervals = intervals_df, values = values)
}

.gtrack_set_created_attrs <- function(trackstr, description, created_by, attrs) {
    .gtrack.attr.set(trackstr, "created.by", created_by, TRUE)
    .gtrack.attr.set(trackstr, "created.date", date(), TRUE)
    .gtrack.attr.set(trackstr, "created.user", Sys.getenv("USER"), TRUE)
    .gtrack.attr.set(trackstr, "description", description, TRUE)

    if (!is.null(attrs)) {
        if (is.null(names(attrs)) || any(names(attrs) == "")) {
            stop("attrs must be a named vector or list", call. = FALSE)
        }
        for (attr_name in names(attrs)) {
            .gtrack.attr.set(trackstr, attr_name, attrs[[attr_name]], FALSE)
        }
    }
}


#' Creates a track from WIG / BigWig / BedGraph / BED / tab-delimited file
#'
#' Creates a track from WIG / BigWig / BedGraph / BED / tab-delimited file
#'
#' This function creates a track from WIG / BigWig / BedGraph / tab-delimited
#' file. Zipped files are supported (file name must have '.gz' or '.zip' suffix).
#'
#' Tab-delimited files must start with a header line with the following column
#' names (tab-separated): 'chrom', 'start', 'end', and exactly one value column
#' name (e.g. 'value'). Each subsequent line provides a single interval:
#' - chrom: chromosome name (e.g. 'chr1')
#' - start: 0-based start coordinate (inclusive)
#' - end: 0-based end coordinate (exclusive)
#' - value: numeric value (floating point allowed); exactly one value column is supported
#'
#' Columns must be separated by tabs. Coordinates must refer to chromosomes
#' existing in the current genome. Missing values can be specified as 'NaN'.
#'
#' BED files (.bed/.bed.gz/.bed.zip) are also supported. If the BED 'score'
#' column (5th column) exists and is numeric, it is used as the interval value;
#' otherwise a constant value of 1 is used. For BED inputs, 'binsize' controls
#' the output type: if 'binsize' is 0 the track is 'Sparse'; otherwise the track
#' is 'Dense' with bin-averaged values based on overlaps with BED intervals (and
#' 'defval' for regions not covered).
#'
#' If 'binsize' is 0 the resulted track is created in 'Sparse' format.
#' Otherwise the 'Dense' format is chosen with a bin size equal to 'binsize'.
#' The values that were not defined in input file file are substituted by
#' 'defval' value.
#'
#' 'description' is added as a track attribute.
#'
#' When multiple databases are connected via \code{\link{gsetroot}}, the track
#' is created in the current working directory (.misha$GWD), which defaults to the
#' last connected database. Use \code{\link{gdir.cd}} with an absolute path to
#' change where new tracks are created.
#'
#' @param track track name
#' @param description a character string description
#' @param file file path
#' @param binsize bin size of the newly created 'Dense' track or '0' for a
#' 'Sparse' track
#' @param defval default track value
#' @param attrs a named vector or list of attributes to be set on the track after import
#' @return None.
#' @seealso \code{\link{gtrack.import_set}}, \code{\link{gtrack.rm}},
#' \code{\link{gtrack.info}}, \code{\link{gdir.create}}, \code{\link{gextract}}
#' @keywords ~wig ~bigwig ~bedgraph ~track
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' \donttest{
#' gdb.init_examples()
#'
#' # Create a simple WIG file for demonstration
#' temp_file <- tempfile(fileext = ".wig")
#' writeLines(c(
#'     "track type=wiggle_0 name=\"example track\"",
#'     "fixedStep chrom=chr1 start=1 step=1",
#'     "1.5",
#'     "2.0",
#'     "1.8",
#'     "3.2"
#' ), temp_file)
#'
#' # Basic import
#' gtrack.import("example_track", "Example track from WIG file",
#'     temp_file,
#'     binsize = 1
#' )
#' gtrack.info("example_track")
#' gtrack.rm("example_track", force = TRUE)
#'
#' # Import with custom attributes
#' attrs <- c("author" = "researcher", "version" = "1.0", "experiment" = "test")
#' gtrack.import("example_track_with_attrs", "Example track with attributes",
#'     temp_file,
#'     binsize = 1, attrs = attrs
#' )
#'
#' # Check that attributes were set
#' gtrack.attr.get("example_track_with_attrs", "author")
#' gtrack.attr.get("example_track_with_attrs", "version")
#' gtrack.attr.get("example_track_with_attrs", "experiment")
#'
#' # Clean up
#' gtrack.rm("example_track_with_attrs", force = TRUE)
#' }
#'
#' @export gtrack.import
gtrack.import <- function(track = NULL, description = NULL, file = NULL, binsize = NULL, defval = NaN, attrs = NULL) {
    if (is.null(substitute(track)) || is.null(description) || is.null(file)) {
        stop("Usage: gtrack.import(track, description, file, binsize, defval = NaN, attrs = NULL)", call. = FALSE)
    }

    .gcheckroot()

    trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())

    # Get creation context (handles prefix resolution)
    ctx <- .gconfirmtrackcreate(trackstr)
    trackdir <- .track_dir(trackstr)
    direxisted <- file.exists(trackdir)
    retv <- 0
    success <- FALSE

    tmp.dirname <- ""
    file.original <- file

    # Execute creation in the target database context
    .gwith_db_context(ctx$db_path, function() {
        tryCatch(
            {
                report.progress <- FALSE
                is_bed_input <- .gtrack_is_bed_path(file)

                if (length(grep("^.+\\.gz$", file, perl = TRUE)) || length(grep("^.+\\.zip$", file, perl = TRUE))) {
                    message("Unzipping...\n")
                    report.progress <- TRUE
                    tmp.dirname <<- tempfile()
                    if (!dir.create(tmp.dirname, recursive = TRUE, mode = "0777")) {
                        stop(sprintf("Failed to create a directory %s", tmp.dirname), call. = FALSE)
                    }

                    file.unzipped <- paste(tmp.dirname, "/", basename(gsub("\\.(gz|zip)$", "", file, perl = TRUE)), sep = "")
                    retv <- system(paste("/bin/sh -c \"gunzip -q -c", file, ">", file.unzipped, "\""), intern = FALSE)
                    if (retv != 0) {
                        stop(sprintf("Failed to unzip file %s", file), call. = FALSE)
                    }
                    file <- file.unzipped
                }

                # BED files can be imported as sparse (binsize==0) or dense (binsize>0)
                if (is_bed_input || length(grep("^.+\\.bed$", file, perl = TRUE))) {
                    message("Importing BED file...\n")
                    report.progress <- TRUE
                    bed_parsed <- .gtrack_read_bed(file)
                    if (!is.null(binsize) && !is.na(binsize) && binsize > 0) {
                        intervalData <- data.frame(
                            chrom = bed_parsed$intervals$chrom,
                            start = bed_parsed$intervals$start,
                            end = bed_parsed$intervals$end,
                            value = as.numeric(bed_parsed$values)
                        )
                        .gcall("gtrack_create_dense", ctx$base_name, intervalData, binsize, defval, .misha_env())
                        .gdb.add_track(ctx$qualified_name)
                        .gtrack_set_created_attrs(ctx$qualified_name, description, sprintf("gtrack.import(%s, description, \"%s\", %d, %g, attrs)", trackstr, file.original, binsize, defval), attrs)
                    } else {
                        .gcall("gtrack_create_sparse", ctx$base_name, bed_parsed$intervals, bed_parsed$values, .misha_env())
                        .gdb.add_track(ctx$qualified_name)
                        .gtrack_set_created_attrs(ctx$qualified_name, description, sprintf("gtrack.import(%s, description, \"%s\", %d, %g, attrs)", trackstr, file.original, 0, defval), attrs)
                    }

                    # If database is indexed, automatically convert the track to indexed format
                    if (.gdb.is_indexed()) {
                        gtrack.convert_to_indexed(ctx$qualified_name)
                    }

                    success <<- TRUE
                } else if (length(grep("^.+\\.bw$", file, perl = TRUE)) || length(grep("^.+\\.bigWig$", file, perl = TRUE)) ||
                    # looks like all bigWig files start with "fc26" in their first two bytes
                    system(sprintf("od -x -N 2 \"%s\"", file), intern = TRUE)[1] == "0000000 fc26") {
                    message("Converting from BigWig to WIG...\n")
                    report.progress <- TRUE
                    if (tmp.dirname == "") {
                        tmp.dirname <<- tempfile()
                        if (!dir.create(tmp.dirname, recursive = TRUE, mode = "0777")) {
                            stop(sprintf("Failed to create a directory %s", tmp.dirname), call. = FALSE)
                        }
                    }

                    file.noext <- basename(gsub("^(.+)\\.(.+)$", "\\1", file, perl = TRUE))
                    file.converted <- paste(tmp.dirname, "/", file.noext, ".wig", sep = "")
                    retv <- system(paste(get_bigWigToWig_bin(), file, file.converted), intern = FALSE)
                    if (retv != 0) {
                        stop("BigWigToWig conversion failed", call. = FALSE)
                    }
                    file <- file.converted
                }

                if (success) {
                    # BED path handled above; skip the generic wig importer
                } else {
                    if (report.progress) {
                        message("Converting to track...\n")
                    }

                    .gcall("gtrackimportwig", ctx$base_name, file, binsize, defval, .misha_env())
                    .gdb.add_track(ctx$qualified_name)
                    .gtrack_set_created_attrs(ctx$qualified_name, description, sprintf("gtrack.import(%s, description, \"%s\", %d, %g, attrs)", trackstr, file.original, binsize, defval), attrs)

                    # If database is indexed, automatically convert the track to indexed format
                    if (.gdb.is_indexed()) {
                        gtrack.convert_to_indexed(ctx$qualified_name)
                    }

                    success <<- TRUE
                }
            },
            finally = {
                if (tmp.dirname != "") {
                    unlink(tmp.dirname, recursive = TRUE)
                }

                if (!success && !direxisted) {
                    unlink(trackdir, recursive = TRUE)
                    .gdb.rm_track(ctx$qualified_name)
                }
            }
        )
    })
    retv <- 0 # suppress return value
}


#' Creates a track from a file of mapped sequences
#'
#' Creates a track from a file of mapped sequences.
#'
#' This function creates a track from a file of mapped sequences. The file can
#' be in SAM format or in a general TAB delimited text format where each line
#' describes a single read.
#'
#' For a SAM file 'cols.order' must be set to 'NULL'.
#'
#' For a general TAB delimited text format the following columns must be
#' presented in the file: sequence, chromosome, coordinate and strand. The
#' position of these columns should be specified in 'cols.order' argument. The
#' default value of 'cols.order' is an array of (9, 11, 13, 14) meaning that
#' sequence is expected to be found at column number 9, chromosome - at column
#' 11, coordinate - at column 13 and strand - at column 14. The column indices
#' are 1-based, i.e. the first column is referenced by 1. Chromosome needs a
#' prefix 'chr' e.g. 'chr1'. Valid strand values are '+' or 'F' for forward
#' strand and '-' or 'R' for the reverse strand.
#'
#' Each read at given coordinate can be "expanded" to cover an interval rather
#' than a single point. The length of the interval is controlled by 'pileup'
#' argument. The direction of expansion depends on the strand value. If
#' 'pileup' is '0', no expansion is performed and the read is converted to a
#' single point. The track is created in sparse format. If 'pileup' is greater
#' than zero, the output track is in dense format. 'binsize' controls the bin
#' size of the dense track.
#'
#' If 'remove.dups' is 'TRUE' the duplicated coordinates are counted only once.
#'
#' 'description' is added as a track attribute.
#'
#' 'gtrack.import_mappedseq' returns the statistics of the conversion process.
#'
#' @param track track name
#' @param description a character string description
#' @param file name of mapped sequences file
#' @param pileup interval expansion
#' @param binsize bin size of a dense track
#' @param cols.order order of sequence, chromosome, coordinate and strand
#' columns in mapped sequences file or NULL if SAM file is used
#' @param remove.dups if 'TRUE' the duplicated coordinates are counted only
#' once.
#' @return A list of conversion process statistics.
#' @seealso \code{\link{gtrack.rm}}, \code{\link{gtrack.info}},
#' \code{\link{gdir.create}}
#' @keywords ~mapped ~sequence ~track
#' @export gtrack.import_mappedseq
gtrack.import_mappedseq <- function(track = NULL, description = NULL, file = NULL, pileup = 0, binsize = -1, cols.order = c(9, 11, 13, 14), remove.dups = TRUE) {
    if (is.null(substitute(track)) || is.null(description) || is.null(file)) {
        stop("Usage: gtrack.import_mappedseq(track, description, file, pileup = 0, binsize = -1, cols.order = c(9, 11, 13, 14), remove.dups = TRUE)", call. = FALSE)
    }
    .gcheckroot()

    trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())

    # Get creation context (handles prefix resolution)
    ctx <- .gconfirmtrackcreate(trackstr)
    trackdir <- .track_dir(trackstr)
    direxisted <- file.exists(trackdir)
    retv <- 0
    success <- FALSE

    # Execute creation in the target database context
    .gwith_db_context(ctx$db_path, function() {
        tryCatch(
            {
                retv <<- .gcall("gtrackimport_mappedseq", ctx$base_name, file, pileup, binsize, cols.order, remove.dups, .misha_env())
                .gdb.add_track(ctx$qualified_name)
                .gtrack.attr.set(
                    ctx$qualified_name, "created.by",
                    sprintf("gtrack.import_mappedseq(%s, description, \"%s\", pileup=%d, binsize=%d, remove.dups=%s)", trackstr, file, pileup, binsize, remove.dups), TRUE
                )
                .gtrack.attr.set(ctx$qualified_name, "created.date", date(), TRUE)
                .gtrack.attr.set(ctx$qualified_name, "created.user", Sys.getenv("USER"), TRUE)
                .gtrack.attr.set(ctx$qualified_name, "description", description, TRUE)
                success <<- TRUE
            },
            finally = {
                if (!success && !direxisted) {
                    unlink(trackdir, recursive = TRUE)
                    .gdb.rm_track(ctx$qualified_name)
                }
            }
        )
    })
    retv
}


#' Creates one or more tracks from multiple WIG / BigWig / BedGraph /
#' tab-delimited files on disk or FTP
#'
#' Creates one or more tracks from WIG / BigWig / BedGraph / tab-delimited
#' files on disk or FTP.
#'
#' This function is similar to 'gtrack.import' however unlike the latter it can
#' create multiple tracks. Additionally the files can be fetched from an FTP
#' server.
#'
#' The files are expected to be in WIG / BigWig / BedGraph / tab-delimited
#' formats. One can learn about the format of the tab-delimited file by running
#' 'gextract' function with a 'file' parameter set to the name of the file.
#' Zipped files are supported (file name must have '.gz' or '.zip' suffix).
#'
#' Files are specified by 'path' argument. 'path' can be also a URL of an FTP
#' server in the form of 'ftp://[address]/[files]'. If 'path' is a URL, the
#' files are first downloaded from FTP server to a temporary directory and then
#' imported to tracks. The temporary directory is created at 'GROOT/downloads'.
#'
#' Regardless whether 'path' is file path or to a URL, it can contain
#' wildcards. Hence multiple files can be imported (and downloaded) at once.
#'
#' If 'binsize' is 0 the resulted tracks are created in 'Sparse' format.
#' Otherwise the 'Dense' format is chosen with a bin size equal to 'binsize'.
#' The values that were not defined in input file file are substituted by
#' 'defval' value.
#'
#' The name of a each created track is of '[track.prefix][filename]' form,
#' where 'filename' is the name of the WIG file. For example, if 'track.prefix'
#' equals to "wigs."" and an input file name is 'mydata', a track named
#' 'wigs.mydata' is created. If 'track.prefix' is 'NULL' no prefix is appended
#' to the name of the created track.
#'
#' Existing tracks are not overwritten and no new directories are automatically
#' created.
#'
#' 'description' is added to the created tracks as an attribute.
#'
#' 'gtrack.import_set' does not stop if an error occurs while importing a file.
#' It rather continues importing the rest of the files.
#'
#' 'gtrack.import_set' returns the names of the files that were successfully
#' imported and those that failed.
#'
#' @param path file path or URL (may contain wildcards)
#' @param description a character string description
#' @param binsize bin size of the newly created 'Dense' track or '0' for a
#' 'Sparse' track
#' @param track.prefix prefix for a track name
#' @param defval default track value
#' @return Names of files that were successfully imported and those that
#' failed.
#' @seealso \code{\link{gtrack.import}}, \code{\link{gwget}},
#' \code{\link{gtrack.rm}}, \code{\link{gtrack.info}},
#' \code{\link{gdir.create}}, \code{\link{gextract}}
#' @keywords ~wig ~bigwig ~bedgraph ~track
#' @export gtrack.import_set
gtrack.import_set <- function(description = NULL, path = NULL, binsize = NULL, track.prefix = NULL, defval = NaN) {
    .gcheckroot()

    if (is.null(description) || is.null(path) || is.null(binsize)) {
        stop("Usage: gtrack.import_set(description, path, binsize, track.prefix = NULL, defval = NaN)", call. = FALSE)
    }

    if (is.null(substitute(track.prefix))) {
        track.prefix <- ""
    } else {
        track.prefix <- do.call(.gexpr2str, list(substitute(track.prefix)), envir = parent.frame())
    }

    files <- c()
    tmp.dirname <- ""

    tryCatch(
        {
            tmp.dirname <- tempfile(pattern = "", tmpdir = paste(get("GROOT", envir = .misha), "/downloads", sep = ""))
            if (!dir.create(tmp.dirname, recursive = TRUE, mode = "0777")) {
                stop(sprintf("Failed to create a directory %s", tmp.dirname), call. = FALSE)
            }
            protocol <- "ftp://"
            if (substr(path, 1, nchar(protocol)) == protocol) {
                # ftp
                files <- gwget(path, tmp.dirname)

                if (!length(files)) {
                    stop("No files downloaded. Exiting.", call. = FALSE)
                }
            } else {
                # local path
                files <- system(paste("/bin/sh -c \"ls -d -A", path, "\""), intern = TRUE)
            }

            files <- files[!file.info(files)$isdir]
            if (!length(files)) {
                stop("No files to import. Exiting.", call. = FALSE)
            }

            files.imported <- c()

            for (file in files) {
                tryCatch(
                    {
                        message(sprintf("Importing file %s", file))
                        file.noext <- basename(gsub("^([^.]+)(\\..*)*$", "\\1", file, perl = TRUE))
                        trackstr <- paste(track.prefix, file.noext, sep = "")

                        .gcall_noninteractive(gtrack.import, trackstr, description, file, binsize, defval)
                        files.imported <- c(files.imported, file)
                        success <- TRUE
                    },
                    error = function(e) {
                        msg <- as.character(e)
                        if (msg == "Error: Command interrupted!\n") {
                            stop("Command interrupted!", call. = FALSE)
                        } else {
                            message(sprintf("%s", msg))
                        }
                    }
                )
            }

            files <- basename(files)
            if (length(files.imported)) {
                files.imported <- basename(files.imported)
            }
            files.failed <- setdiff(files, files.imported)
            res <- new.env()
            if (length(files.failed)) {
                res$files.failed <- files.failed
            }
            if (length(files.imported)) {
                res$files.imported <- files.imported
            }
            as.list(res)
        },
        finally = {
            unlink(tmp.dirname, recursive = TRUE)
        }
    )
}
