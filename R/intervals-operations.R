# Interval operations (intersect, union, diff, canonic, etc.)

#' Intersects two-dimensional intervals with a band
#'
#' Intersects two-dimensional intervals with a band.
#'
#' This function intersects each two-dimensional interval from 'intervals' with
#' 'band'. If the intersection is not empty, the interval is shrunk to the
#' minimal rectangle that contains the band and added to the return value.
#'
#' If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
#' set. Use this parameter if the result size exceeds the limits of the
#' physical memory.
#'
#' @param intervals two-dimensional intervals
#' @param band track expression band. If 'NULL' no band is used.
#' @param intervals.set.out intervals set name where the function result is
#' optionally outputted
#' @return If 'intervals.set.out' is 'NULL' a data frame representing the
#' intervals.
#' @seealso \code{\link{gintervals.2d}}, \code{\link{gintervals.intersect}}
#' @keywords ~band ~intersect
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gintervals.2d.band_intersect(gintervals.2d(1), c(10000, 20000))
#'
#' @export gintervals.2d.band_intersect
gintervals.2d.band_intersect <- function(intervals = NULL, band = NULL, intervals.set.out = NULL) {
    if (is.null(intervals)) {
        stop("Usage: gintervals.2d.band_intersect(intervals, band = NULL, intervals.set.out = NULL)", call. = FALSE)
    }

    intervals <- rescue_ALLGENOME(intervals, as.character(substitute(intervals)))

    intervals.set.out <- do.call(.gexpr2str, list(substitute(intervals.set.out)), envir = parent.frame())

    if (!is.null(intervals.set.out)) {
        fullpath <- .gintervals.check_new_set(intervals.set.out)
    }

    # intervals can be NULL if gextract is piped with gscreen and the latter returns NULL
    success <- FALSE
    res <- NULL
    tryCatch(
        {
            if (!is.null(intervals)) {
                res <- .gcall("ginterv_intersectband", intervals, band, intervals.set.out, .misha_env())
                if (!is.null(intervals.set.out) && .gintervals.is_bigset(intervals.set.out, FALSE) && !.gintervals.needs_bigset(intervals.set.out)) {
                    .gintervals.big2small(intervals.set.out)
                }
            }
            success <- TRUE
        },
        finally = {
            if (!success && !is.null(intervals.set.out)) {
                unlink(fullpath, recursive = TRUE)
            }
        }
    )

    # refresh the list of GINTERVS, etc.
    if (!is.null(intervals.set.out)) {
        .gdb.add_intervals.set(intervals.set.out)
        retv <- 0 # suppress return value
    } else {
        res
    }
}


#' Converts intervals to canonic form
#'
#' Converts intervals to canonic form.
#'
#' This function converts 'intervals' into a "canonic" form: properly sorted
#' with no overlaps. The result can be used later in the functions that require
#' the intervals to be in canonic form. Use 'unify_touching_intervals' to
#' control whether the intervals that touch each other (i.e. the end coordinate
#' of one equals to the start coordinate of the other) are unified.
#' 'unify_touching_intervals' is ignored if two-dimensional intervals are used.
#'
#' Since 'gintervals.canonic' unifies overlapping or touching intervals, the
#' number of the returned intervals might be less than the number of the
#' original intervals. To allow the user to find the origin of the new interval
#' 'mapping' attribute is attached to the result. It maps between the original
#' intervals and the resulted intervals. Use 'attr(retv_of_gintervals.canonic,
#' "mapping")' to retrieve the map.
#'
#' @param intervals intervals to be converted
#' @param unify_touching_intervals if 'TRUE', touching one-dimensional
#' intervals are unified
#' @return A data frame representing the canonic intervals and an attribute
#' 'mapping' that maps the original intervals to the resulted ones.
#' @seealso \code{\link{gintervals}}, \code{\link{gintervals.2d}}
#' @keywords ~intervals ~canonic
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#'
#' ## Create intervals manually by using 'data.frame'.
#' ## Note that we add an additional column 'data'.
#' ## Return value:
#' ##   chrom start   end data
#' ## 1  chr1 11000 12000   10
#' ## 2  chr1   100   200   20
#' ## 3  chr1 10000 13000   30
#' ## 4  chr1 10500 10600   40
#' intervs <- data.frame(
#'     chrom = "chr1",
#'     start = c(11000, 100, 10000, 10500),
#'     end = c(12000, 200, 13000, 10600),
#'     data = c(10, 20, 30, 40)
#' )
#'
#' ## Convert the intervals into the canonic form.
#' ## The function discards any columns besides chrom, start and end.
#' ## Return value:
#' ##  chrom start   end
#' ## 1  chr1   100   200
#' ## 2  chr1 10000 13000
#' res <- gintervals.canonic(intervs)
#'
#' ## By inspecting mapping attribute we can see how the new
#' ## intervals were created: "2 1 2 2" means that the first
#' ## interval in the result was created from the second interval in
#' ## the original set (we look for the indices in mapping where "1"
#' ## appears). Likewise the second interval in the result was
#' ## created from 3 intervals in the original set. Their indices are
#' ## 1, 3 and 4 (once again we look for the indices in mapping where
#' ## "2" appears).
#' ## Return value:
#' ## 2 1 2 2
#' attr(res, "mapping")
#'
#' ## Finally (and that is the most useful part of 'mapping'
#' ## attribute): we add a new column 'data' to our result which is
#' ## the mean value of the original data column. The trick is done
#' ## using 'tapply' on par with 'mapping' attribute. For example,
#' ## 20.00000 equals is a result of 'mean(intervs[2,]$data' while
#' ## 26.66667 is a result of 'mean(intervs[c(1,3,4),]$data)'.
#' ## 'res' after the following call:
#' ##  chrom start   end     data
#' ## 1  chr1   100   200 20.00000
#' ## 2  chr1 10000 13000 26.66667
#' res$data <- tapply(intervs$data, attr(res, "mapping"), mean)
#'
#' @export gintervals.canonic
gintervals.canonic <- function(intervals = NULL, unify_touching_intervals = TRUE) {
    if (is.null(intervals)) {
        stop("Usage: gintervals.canonic(intervals, unify_touching_intervals = TRUE)", call. = FALSE)
    }

    intervals <- rescue_ALLGENOME(intervals, as.character(substitute(intervals)))

    res <- .gcall("gintervcanonic", intervals, unify_touching_intervals, .misha_env())
    res
}

#' Mark overlapping intervals with a group ID
#'
#' @param intervals intervals set
#' @param group_col name of the column to store the overlap group IDs (default: "overlap_group")
#' @return The intervals set with an additional column containing group IDs
#' from gintervals.canonic mapping. All overlapping intervals will have the same group ID.
#'
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' # Create sample overlapping intervals
#' intervs <- data.frame(
#'     chrom = "chr1",
#'     start = c(11000, 100, 10000, 10500),
#'     end = c(12000, 200, 13000, 10600),
#'     data = c(10, 20, 30, 40)
#' )
#'
#' # Mark overlapping intervals
#' intervs_marked <- gintervals.mark_overlaps(intervs)
#'
#' # Use custom column name
#' intervs_marked <- gintervals.mark_overlaps(intervs, group_col = "my_groups")
#' @inheritParams gintervals.canonic
#' @export
gintervals.mark_overlaps <- function(intervals, group_col = "overlap_group", unify_touching_intervals = TRUE) {
    canon <- gintervals.canonic(intervals, unify_touching_intervals = unify_touching_intervals)
    mapping <- attr(canon, "mapping")
    intervals[, group_col] <- mapping

    return(intervals)
}


#' Calculates difference of two intervals sets
#'
#' Returns difference of two sets of intervals.
#'
#' This function returns a genomic space that is covered by 'intervals1' but
#' not covered by 'intervals2'.
#'
#' If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
#' set. Use this parameter if the result size exceeds the limits of the
#' physical memory.
#'
#' @param intervals1,intervals2 set of one-dimensional intervals
#' @param intervals.set.out intervals set name where the function result is
#' optionally outputted
#' @return If 'intervals.set.out' is 'NULL' a data frame representing the
#' intervals.
#' @seealso \code{\link{gintervals}}, \code{\link{gintervals.intersect}},
#' \code{\link{gintervals.union}}
#' @keywords ~diff
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#'
#' intervs1 <- gscreen("dense_track > 0.15")
#' intervs2 <- gscreen("dense_track < 0.2")
#'
#' ## 'res3' equals to 'res4'
#' res3 <- gintervals.diff(intervs1, intervs2)
#' res4 <- gscreen("dense_track >= 0.2")
#'
#' @export gintervals.diff
gintervals.diff <- function(intervals1 = NULL, intervals2 = NULL, intervals.set.out = NULL) {
    if (is.null(intervals1) || is.null(intervals2)) {
        stop("Usage: gintervals.diff(intervals1, intervals2, intervals.set.out = NULL)", call. = FALSE)
    }

    intervals.set.out <- do.call(.gexpr2str, list(substitute(intervals.set.out)), envir = parent.frame())

    if (.gintervals.is_bigset(intervals1) || .gintervals.is_bigset(intervals2) || !is.null(intervals.set.out)) {
        res <- NULL

        FUN <- function(intervals, intervals.set.out, envir) {
            intervals1 <- intervals[[1]]
            intervals2 <- intervals[[2]]
            chrom_res <- .gcall("gintervdiff", intervals1, intervals2, .misha_env())
            if (!is.null(chrom_res) && nrow(chrom_res) > 0) {
                if (is.null(intervals.set.out)) {
                    assign("res", c(get("res", envir = envir), list(chrom_res)), envir = envir)
                    .gverify_max_data_size(sum(unlist(lapply(get("res", envir), nrow))), arguments = "intervals.set.out")
                }
            }
            chrom_res
        }

        .gintervals.apply(gintervals.chrom_sizes(intervals1), list(intervals1, intervals2), intervals.set.out, FUN, intervals.set.out, environment())

        if (!is.null(res)) {
            res <- do.call(.grbind, res)
        } # much faster than calling rbind incrementally in FUN

        if (is.null(intervals.set.out)) {
            if (!is.null(res) && nrow(res)) {
                res
            } else {
                NULL
            }
        } else {
            retv <- 0
        } # suppress return value
    } else {
        res <- .gcall("gintervdiff", intervals1, intervals2, .misha_env())
        res
    }
}


#' Tests for a named intervals set existence
#'
#' Tests for a named intervals set existence.
#'
#' This function returns 'TRUE' if a named intervals set exists in Genomic
#' Database.
#'
#' @param intervals.set name of an intervals set
#' @return 'TRUE' if a named intervals set exists. Otherwise 'FALSE'.
#' @seealso \code{\link{gintervals.ls}}, \code{\link{gintervals.load}},
#' \code{\link{gintervals.rm}}, \code{\link{gintervals.save}},
#' \code{\link{gintervals}}, \code{\link{gintervals.2d}}
#' @keywords ~intervals
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gintervals.exists("annotations")
#'
#' @export gintervals.exists
gintervals.exists <- function(intervals.set = NULL) {
    if (is.null(substitute(intervals.set))) {
        stop("Usage: gintervals.exists(intervals.set)", call. = FALSE)
    }
    .gcheckroot()

    intervals.set <- do.call(.gexpr2str, list(substitute(intervals.set)), envir = parent.frame())
    !is.na(match(intervals.set, get("GINTERVS", envir = .misha)))
}


#' Returns the path on disk of an interval set
#'
#' Returns the path on disk of an interval set.
#'
#' This function returns the actual file system path where an interval set is stored.
#' The function works with a single interval set name or a vector of names.
#'
#' @param intervals.set name of an interval set or a vector of interval set names
#' @return A character vector containing the full paths to the interval sets on disk.
#' @seealso \code{\link{gintervals.exists}}, \code{\link{gintervals.ls}},
#' \code{\link{gtrack.path}}
#' @keywords ~intervals ~path
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gintervals.path("annotations")
#' gintervals.path(c("annotations", "coding"))
#'
#' @export gintervals.path
gintervals.path <- function(intervals.set = NULL) {
    if (is.null(substitute(intervals.set))) {
        stop("Usage: gintervals.path(intervals.set)", call. = FALSE)
    }
    .gcheckroot()

    intervals.set <- do.call(.gexpr2str, list(substitute(intervals.set)), envir = parent.frame())

    # Handle vectorized input
    if (length(intervals.set) == 0) {
        return(character(0))
    }

    # Construct paths for each interval set
    paths <- vapply(intervals.set, function(interv) {
        path <- gsub("\\.", "/", interv)
        sprintf("%s.interv", paste(get("GWD", envir = .misha), path, sep = "/"))
    }, character(1), USE.NAMES = FALSE)

    paths
}


#' Limits intervals to chromosomal range
#'
#' Limits intervals to chromosomal range.
#'
#' This function enforces the intervals to be within the chromosomal range [0,
#' chrom length) by altering the intervals' boundaries. Intervals that lay
#' entirely outside of the chromosomal range are eliminated. The new intervals
#' are returned.
#'
#' If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
#' set. Use this parameter if the result size exceeds the limits of the
#' physical memory.
#'
#' @param intervals intervals
#' @param intervals.set.out intervals set name where the function result is
#' optionally outputted
#' @return If 'intervals.set.out' is 'NULL' a data frame representing the
#' intervals.
#' @seealso \code{\link{gintervals}}, \code{\link{gintervals.2d}},
#' \code{\link{gintervals.canonic}}
#' @keywords ~intervals
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' intervs <- data.frame(
#'     chrom = "chr1",
#'     start = c(11000, -100, 10000, 10500),
#'     end = c(12000, 200, 13000000, 10600)
#' )
#' gintervals.force_range(intervs)
#'
#' @export gintervals.force_range
gintervals.force_range <- function(intervals = NULL, intervals.set.out = NULL) {
    if (is.null(substitute(intervals))) {
        stop("Usage: gintervals.force_range(intervals, intervals.set.out = NULL)", call. = FALSE)
    }
    .gcheckroot()

    intervals <- rescue_ALLGENOME(intervals, as.character(substitute(intervals)))

    intervals.set.out <- do.call(.gexpr2str, list(substitute(intervals.set.out)), envir = parent.frame())

    res <- NULL
    FUN <- function(intervals, intervals.set.out, envir) {
        intervals <- intervals[[1]]
        if (.gintervals.is1d(intervals)) {
            intervals$start <- pmax(intervals$start, 0)
            intervals$end <- pmin(intervals$end, gintervals.all()[match(intervals$chrom, gintervals.all()$chrom), ]$end)
            intervals <- intervals[!is.na(intervals$end) & intervals$start < intervals$end, ]
        } else {
            intervals$start1 <- pmax(intervals$start1, 0)
            intervals$end1 <- pmin(intervals$end1, gintervals.all()[match(intervals$chrom1, gintervals.all()$chrom), ]$end)
            intervals$start2 <- pmax(intervals$start2, 0)
            intervals$end2 <- pmin(intervals$end2, gintervals.all()[match(intervals$chrom2, gintervals.all()$chrom), ]$end)
            intervals <- intervals[!is.na(intervals$end1) & !is.na(intervals$end2) & intervals$start1 < intervals$end1 & intervals$start2 < intervals$end2, ]
        }
        if (is.null(intervals.set.out)) {
            assign("res", c(get("res", envir = envir), list(intervals)), envir = envir)
            .gverify_max_data_size(sum(unlist(lapply(get("res", envir), nrow))), arguments = "intervals.set.out")
        }
        intervals
    }

    .gintervals.apply(gintervals.chrom_sizes(intervals), intervals, intervals.set.out, FUN, intervals.set.out, environment())

    if (!is.null(res)) {
        res <- do.call(.grbind, res)
    } # much faster than calling rbind incrementally in FUN

    if (is.null(intervals.set.out)) {
        if (!is.null(res) && nrow(res)) {
            res
        } else {
            NULL
        }
    } else {
        retv <- 0
    } # suppress return value
}

#' Normalize intervals to a fixed size
#'
#' This function normalizes intervals by computing their centers and then expanding
#' them to a fixed size, while ensuring they don't cross chromosome boundaries.
#'
#' @param intervals intervals set
#' @param size target size for normalized intervals (must be positive integer)
#' @param intervals.set.out intervals set name where the function result is saved.
#' If NULL, the result is returned to the user.
#' @return Normalized intervals set with fixed size, or NULL if result is saved to intervals.set.out
#' @seealso \code{\link{gintervals.force_range}}
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' intervs <- gintervals(1, c(1000, 5000), c(2000, 6000))
#' gintervals.normalize(intervs, 500)
#'
#' @export gintervals.normalize
gintervals.normalize <- function(intervals = NULL, size = NULL, intervals.set.out = NULL) {
    if (is.null(substitute(intervals)) || is.null(size)) {
        stop("Usage: gintervals.normalize(intervals, size, intervals.set.out = NULL)", call. = FALSE)
    }
    .gcheckroot()

    if (!is.numeric(size) || length(size) != 1 || size <= 0) {
        stop("Size must be a positive number", call. = FALSE)
    }

    intervals <- rescue_ALLGENOME(intervals, as.character(substitute(intervals)))
    intervals.set.out <- do.call(.gexpr2str, list(substitute(intervals.set.out)), envir = parent.frame())

    res <- NULL
    FUN <- function(intervals, intervals.set.out, envir) {
        intervals <- intervals[[1]]
        if (.gintervals.is2d(intervals)) {
            stop("gintervals.normalize does not support 2D intervals", call. = FALSE)
        }

        # Get original column order and additional columns
        original_cols <- names(intervals)
        basic_cols <- c("chrom", "start", "end", "strand")
        extra_cols <- setdiff(original_cols, basic_cols)

        normalized <- .gcall("gintervals_normalize", intervals, as.integer(size), .misha_env())

        # Preserve additional columns from original intervals in the correct order
        if (!is.null(normalized) && nrow(normalized) > 0 && length(extra_cols) > 0) {
            # Ensure we have the right number of rows
            if (nrow(normalized) == nrow(intervals)) {
                for (col in extra_cols) {
                    normalized[[col]] <- intervals[[col]]
                }
            }
        }

        # Also ensure strand column exists if it was in original
        if (!is.null(normalized) && nrow(normalized) > 0 && "strand" %in% names(intervals) && !"strand" %in% names(normalized)) {
            normalized$strand <- intervals$strand
        }

        # Reorder columns to match original order
        if (!is.null(normalized) && nrow(normalized) > 0) {
            # Get the columns that exist in normalized result
            available_cols <- intersect(original_cols, names(normalized))
            # Reorder to match original order
            normalized <- normalized[, available_cols, drop = FALSE]
        }

        if (is.null(intervals.set.out)) {
            assign("res", c(get("res", envir = envir), list(normalized)), envir = envir)
            .gverify_max_data_size(sum(unlist(lapply(get("res", envir), nrow))), arguments = "intervals.set.out")
        }
        normalized
    }

    .gintervals.apply(gintervals.chrom_sizes(intervals), intervals, intervals.set.out, FUN, intervals.set.out, environment())

    if (!is.null(res) && length(res) > 0) {
        # Filter out NULL results before rbind
        res <- res[!sapply(res, is.null)]
        if (length(res) > 0) {
            res <- do.call(.grbind, res)
        } else {
            res <- NULL
        }
    } # much faster than calling rbind incrementally in FUN

    if (is.null(intervals.set.out)) {
        if (!is.null(res) && nrow(res)) {
            res
        } else {
            NULL
        }
    } else {
        retv <- 0
    } # suppress return value
}

#' Imports genes and annotations from files
#'
#' Imports genes and annotations from files.
#'
#' This function reads a definition of genes from 'genes.file' and returns four
#' sets of intervals: TSS, exons, 3utr and 5utr. In addition to the regular
#' intervals columns 'strand' column is added. It contains '1' values for '+'
#' strands and '-1' values for '-' strands.
#'
#' If annotation file 'annots.file' is given then annotations are attached too
#' to the intervals. The names of the annotations as they would appear in the
#' return value must be specified in 'annots.names' argument.
#'
#' Both 'genes.file' and 'annots.file' can be either a file path or URL in a
#' form of 'ftp://[address]/[file]'. Files that these arguments point to can be
#' zipped or unzipped.
#'
#' Examples of 'genes.file' and 'annots.file' can be found here:
#'
#' \code{ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz}
#' \code{ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/kgXref.txt.gz}
#'
#' If a few intervals overlap (for example: two TSS regions) they are all
#' unified to an interval that covers the whole overlapping region. 'strand'
#' value is set to '0' if two or more of the overlapping intervals have
#' different strands. The annotations of the overlapping intervals are
#' concatenated to a single character string separated by semicolons. Identical
#' values of overlapping intervals' annotation are eliminated.
#'
#' @param genes.file name or URL of file that contains genes
#' @param annots.file name of URL file that contains annotations. If 'NULL' no
#' annotations are imported
#' @param annots.names annotations names
#' @return A list of four intervals sets named 'tss', 'exons', 'utr3' and
#' 'utr5'. 'strand' column and annotations are attached to the intevals.
#' @seealso \code{\link{gintervals}}, \code{\link{gdb.create}}
#' @keywords ~intervals ~import ~genes
#' @export gintervals.import_genes
gintervals.import_genes <- function(genes.file = NULL, annots.file = NULL, annots.names = NULL) {
    if (is.null(genes.file)) {
        stop("Usage: gintervals.import_genes(genes.file, annots.file = NULL, annots.names = NULL)", call. = FALSE)
    }

    tmp.dirname <- tempfile(pattern = "", tmpdir = paste(get("GROOT", envir = .misha), "/downloads", sep = ""))
    if (!dir.create(tmp.dirname, recursive = TRUE, mode = "0777")) {
        stop(sprintf("Failed to create a directory %s", tmp.dirname), call. = FALSE)
    }

    files <- list(genes.file, annots.file)
    file.types <- c("genes.file", "annots.file")

    tryCatch(
        {
            for (i in 1:length(files)) {
                if (is.null(files[[i]])) {
                    next
                }

                protocol <- "ftp://"
                if (substr(files[[i]], 1, nchar(protocol)) == protocol) {
                    # ftp
                    f <- gwget(files[[i]], tmp.dirname)
                    if (length(f) != 1) {
                        stop(sprintf("More than one file matches %s argument", file.types[i]), call. = FALSE)
                    }
                    files[[i]] <- f
                }

                if (length(grep("^.+\\.gz$", files[[i]], perl = TRUE))) {
                    f.unzipped <- basename(gsub("^(.+)\\.gz$", "\\1", files[[i]], perl = TRUE))
                    f.unzipped <- paste(tmp.dirname, "/", f.unzipped, sep = "")
                    cmd <- paste("/bin/sh -c \"gunzip -q -c", files[[i]], ">", f.unzipped, "\"")
                    if (system(cmd)) {
                        stop(sprintf("Command failed: %s", cmd), call. = FALSE)
                    }
                    files[[i]] <- f.unzipped
                }
            }

            res <- .gcall("gintervals_import_genes", files[[1]], files[[2]], annots.names, .misha_env())
            res
        },
        finally = {
            unlink(tmp.dirname, recursive = TRUE)
        }
    )
}


#' Calculates an intersection of two sets of intervals
#'
#' Calculates an intersection of two sets of intervals.
#'
#' This function returns intervals that represent a genomic space which is
#' achieved by intersection of 'intervals1' and 'intervals2'.
#'
#' If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
#' set. Use this parameter if the result size exceeds the limits of the
#' physical memory.
#'
#' @param intervals1,intervals2 set of intervals
#' @param intervals.set.out intervals set name where the function result is
#' optionally outputted
#' @return If 'intervals.set.out' is 'NULL' a data frame representing the
#' intersection of intervals.
#' @seealso \code{\link{gintervals.2d.band_intersect}},
#' \code{\link{gintervals.diff}}, \code{\link{gintervals.union}},
#' \code{\link{gintervals}}, \code{\link{gintervals.2d}}
#' @keywords ~intersect
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#'
#' intervs1 <- gscreen("dense_track > 0.15")
#' intervs2 <- gscreen("dense_track < 0.2")
#'
#' ## 'intervs3' and 'intervs4' are identical
#' intervs3 <- gintervals.intersect(intervs1, intervs2)
#' intervs4 <- gscreen("dense_track > 0.15 & dense_track < 0.2")
#'
#' @export gintervals.intersect
gintervals.intersect <- function(intervals1 = NULL, intervals2 = NULL, intervals.set.out = NULL) {
    if (is.null(intervals1) || is.null(intervals2)) {
        stop("Usage: gintervals.intersect(intervals1, intervals2, intervals.set.out = NULL)", call. = FALSE)
    }

    intervals.set.out <- do.call(.gexpr2str, list(substitute(intervals.set.out)), envir = parent.frame())

    if (.gintervals.is_bigset(intervals1) || .gintervals.is_bigset(intervals2) || !is.null(intervals.set.out)) {
        res <- NULL

        FUN <- function(intervals, intervals.set.out, envir) {
            intervals1 <- intervals[[1]]
            intervals2 <- intervals[[2]]
            chrom_res <- .gcall("gintervintersect", intervals1, intervals2, .misha_env())
            if (!is.null(chrom_res) && nrow(chrom_res) > 0) {
                if (is.null(intervals.set.out)) {
                    assign("res", c(get("res", envir = envir), list(chrom_res)), envir = envir)
                    .gverify_max_data_size(sum(unlist(lapply(get("res", envir), nrow))), arguments = "intervals.set.out")
                }
            }
            chrom_res
        }

        chroms1 <- gintervals.chrom_sizes(intervals1)
        chroms1$size <- NULL
        chroms2 <- gintervals.chrom_sizes(intervals2)
        chroms2$size <- NULL
        .gintervals.apply(merge(chroms1, chroms2), list(intervals1, intervals2), intervals.set.out, FUN, intervals.set.out, environment())

        if (!is.null(res)) {
            res <- do.call(.grbind, res)
        } # much faster than calling rbind incrementally in FUN

        if (is.null(intervals.set.out)) {
            if (!is.null(res) && nrow(res)) {
                res
            } else {
                NULL
            }
        } else {
            retv <- 0
        } # suppress return value
    } else {
        res <- .gcall("gintervintersect", intervals1, intervals2, .misha_env())
        res
    }
}
#' Returns number of intervals per chromosome
#'
#' Returns number of intervals per chromosome (or chromosome pair).
#'
#' This function returns number of intervals per chromosome (for 1D intervals)
#' or chromosome pair (for 2D intervals).
#'
#' @param intervals intervals set
#' @return Data frame representing number of intervals per chromosome (for 1D
#' intervals) or chromosome pair (for 2D intervals).
#' @seealso \code{\link{gintervals.load}}, \code{\link{gintervals.save}},
#' \code{\link{gintervals.exists}}, \code{\link{gintervals.ls}},
#' \code{\link{gintervals}}, \code{\link{gintervals.2d}}
#' @keywords ~intervals
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gintervals.chrom_sizes("annotations")
#'
#' @export gintervals.chrom_sizes
gintervals.chrom_sizes <- function(intervals = NULL) {
    if (is.null(intervals)) {
        stop("Usage: gintervals.chrom_sizes(intervals)", call. = FALSE)
    }
    .gcheckroot()

    intervals <- rescue_ALLGENOME(intervals, as.character(substitute(intervals)))

    if (.gintervals.is_bigset(intervals)) {
        if (.gintervals.big.is1d(intervals)) {
            stats <- .gintervals.big.meta(intervals)$stats
            res <- stats[, match(c("chrom", "size"), colnames(stats))]
        } else {
            stats <- .gintervals.big.meta(intervals)$stats
            res <- stats[, match(c("chrom1", "chrom2", "size"), colnames(stats))]
        }
    } else {
        res <- .gcall("gintervals_chrom_sizes", .gintervals.load_ext(intervals), .misha_env())
    }

    if (nrow(res) > 1) {
        rownames(res) <- 1:nrow(res)
    }
    res
}


#' Tests for big intervals set
#'
#' Tests for big intervals set.
#'
#' This function tests whether 'intervals.set' is a big intervals set.
#' Intervals set is big if it is stored in big intervals set format and given
#' the current limits it cannot be fully loaded into memory.
#'
#' Memory limit is controlled by 'gmax.data.size' option (see:
#' 'getOption("gmax.data.size")').
#'
#' @param intervals.set name of an intervals set
#' @return 'TRUE' if intervals set is big, otherwise 'FALSE'.
#' @seealso \code{\link{gintervals.load}}, \code{\link{gintervals.save}},
#' \code{\link{gintervals.exists}}, \code{\link{gintervals.ls}}
#' @keywords ~intervals
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gintervals.is.bigset("annotations")
#'
#' @export gintervals.is.bigset
gintervals.is.bigset <- function(intervals.set = NULL) {
    if (is.null(intervals.set)) {
        stop("Usage: gintervals.is.bigset(intervals.set)", call. = FALSE)
    }
    .gcheckroot()

    .gintervals.is_bigset(intervals.set) && !.gintervals.loadable(intervals.set)
}


#' Calculates a union of two sets of intervals
#'
#' Calculates a union of two sets of intervals.
#'
#' This function returns intervals that represent a genomic space covered by
#' either 'intervals1' or 'intervals2'.
#'
#' If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
#' set. Use this parameter if the result size exceeds the limits of the
#' physical memory.
#'
#' @param intervals1,intervals2 set of one-dimensional intervals
#' @param intervals.set.out intervals set name where the function result is
#' optionally outputted
#' @return If 'intervals.set.out' is 'NULL' a data frame representing the union
#' of intervals.
#' @seealso \code{\link{gintervals.intersect}}, \code{\link{gintervals.diff}},
#' \code{\link{gintervals}}, \code{\link{gintervals.2d}}
#' @keywords ~union
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#'
#' intervs1 <- gscreen("dense_track > 0.15 & dense_track < 0.18")
#' intervs2 <- gscreen("dense_track >= 0.18 & dense_track < 0.2")
#'
#' ## 'intervs3' and 'intervs4' are identical
#' intervs3 <- gintervals.union(intervs1, intervs2)
#' intervs4 <- gscreen("dense_track > 0.15 & dense_track < 0.2")
#'
#' @export gintervals.union
gintervals.union <- function(intervals1 = NULL, intervals2 = NULL, intervals.set.out = NULL) {
    if (is.null(intervals1) || is.null(intervals2)) {
        stop("Usage: gintervals.union(intervals1, intervals2, intervals.set.out = NULL)", call. = FALSE)
    }

    if (.gintervals.is_bigset(intervals1) || .gintervals.is_bigset(intervals2) || !is.null(intervals.set.out)) {
        res <- NULL

        FUN <- function(intervals, intervals.set.out, envir) {
            intervals1 <- intervals[[1]]
            intervals2 <- intervals[[2]]
            chrom_res <- .gcall("gintervunion", intervals1, intervals2, .misha_env())
            if (!is.null(chrom_res) && nrow(chrom_res) > 0) {
                if (is.null(intervals.set.out)) {
                    assign("res", c(get("res", envir = envir), list(chrom_res)), envir = envir)
                    .gverify_max_data_size(sum(unlist(lapply(get("res", envir), nrow))), arguments = "intervals.set.out")
                }
            }
            chrom_res
        }

        # use for .gintervals.apply chromosomes from both intervals1 and intervals2
        chroms <- NULL
        if (.gintervals.is1d(intervals1)) {
            chroms <- rbind(gintervals.chrom_sizes(intervals1), gintervals.chrom_sizes(intervals2))
            chroms <- factor(chroms$chrom, levels(chroms$chrom))
            chroms <- unique(chroms)
            chroms <- sort(chroms)
            chroms <- data.frame(chrom = chroms)
        } else {
            chroms <- rbind(gintervals.chrom_sizes(intervals1), gintervals.chrom_sizes(intervals2))
            chroms <- data.frame(chrom1 = chroms$chrom1, chrom2 = chroms$chrom2)
            chroms <- unique(chroms)
            chroms <- chroms[with(chroms, order(chrom1, chrom2)), ]
        }
        .gintervals.apply(chroms, list(intervals1, intervals2), intervals.set.out, FUN, intervals.set.out, environment())

        if (!is.null(res)) {
            res <- do.call(.grbind, res)
        } # much faster than calling rbind incrementally in FUN

        if (is.null(intervals.set.out)) {
            if (!is.null(res) && nrow(res)) {
                res
            } else {
                NULL
            }
        } else {
            retv <- 0
        } # suppress return value
    } else {
        res <- .gcall("gintervunion", intervals1, intervals2, .misha_env())
        res
    }
}
