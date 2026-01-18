# Interval set management (ls, load, save, rm, update)

#' Returns a list of named intervals sets
#'
#' Returns a list of named intervals sets in Genomic Database.
#'
#' This function returns a list of named intervals sets that match the pattern
#' (see 'grep'). If called without any arguments all named intervals sets are
#' returned.
#'
#' When multiple databases are connected, the 'db' parameter can be used to
#' filter intervals to only those from a specific database.
#'
#' @param pattern,ignore.case,perl,fixed,useBytes see 'grep'
#' @param db optional database path to filter intervals. If specified, only
#' interval sets from that database are returned.
#' @return An array that contains the names of intervals sets.
#' @seealso \code{\link{grep}}, \code{\link{gintervals.exists}},
#' \code{\link{gintervals.load}}, \code{\link{gintervals.save}},
#' \code{\link{gintervals.rm}}, \code{\link{gintervals}},
#' \code{\link{gintervals.2d}}, \code{\link{gintervals.dataset}}
#' @keywords ~intervals ~ls
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gintervals.ls()
#' gintervals.ls(pattern = "annot*")
#'
#' @export gintervals.ls
gintervals.ls <- function(pattern = "", db = NULL, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE) {
    .gcheckroot()

    intervals <- get("GINTERVS", envir = .misha)

    # Filter by database if specified
    if (!is.null(db)) {
        db <- normalizePath(db, mustWork = FALSE)
        intervals_db <- get("GINTERVALS_DATASET", envir = .misha)
        if (is.null(intervals_db) || length(intervals) == 0) {
            if (!identical(db, get("GROOT", envir = .misha))) {
                return(character(0))
            }
        }
        db_by_intervals <- intervals_db[intervals]
        intervals <- intervals[!is.na(db_by_intervals) & db_by_intervals == db]
        if (length(intervals) == 0) {
            return(character(0))
        }
    }

    grep(pattern, intervals, value = TRUE, ignore.case = ignore.case, perl = perl, fixed = fixed, useBytes = useBytes)
}

#' Returns the database/dataset path for interval sets
#'
#' Returns the path of the database or dataset containing an interval set.
#'
#' When datasets are loaded, interval sets can come from either the working database
#' or from loaded datasets. This function returns the source path for each interval set.
#'
#' @param intervals interval set name or a vector of interval set names
#' @return A character vector containing the database paths for each interval set.
#' Returns NA for interval sets that don't exist in any connected database.
#' @seealso \code{\link{gintervals.dbs}}, \code{\link{gintervals.exists}},
#' \code{\link{gintervals.ls}}, \code{\link{gdataset.ls}}
#' @keywords ~intervals ~path ~database
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gintervals.dataset("annotations1")
#'
#' @export gintervals.dataset
gintervals.dataset <- function(intervals = NULL) {
    if (is.null(substitute(intervals))) {
        stop("Usage: gintervals.dataset(intervals)", call. = FALSE)
    }
    .gcheckroot()

    intervalsstr <- do.call(.gexpr2str, list(substitute(intervals)), envir = parent.frame())
    if (length(intervalsstr) == 0) {
        return(character(0))
    }

    intervals_db <- get("GINTERVALS_DATASET", envir = .misha)
    if (is.null(intervals_db) || length(intervals_db) == 0) {
        intervals_all <- get("GINTERVS", envir = .misha)
        groot <- get("GROOT", envir = .misha)
        return(ifelse(intervalsstr %in% intervals_all, groot, NA_character_))
    }

    unname(intervals_db[intervalsstr])
}

#' Returns all database paths containing an interval set
#'
#' Returns all database paths that contain a version of an interval set.
#'
#' When datasets are loaded, an interval set may exist in multiple locations.
#' This function computes on-demand and returns all such paths.
#'
#' @param intervals interval set name
#' @param dataframe return a data frame with columns \code{intervals} and \code{db}
#' @return A named character vector of database paths. If \code{dataframe} is TRUE,
#' returns a data frame with columns \code{intervals} and \code{db}.
#' @seealso \code{\link{gintervals.dataset}}, \code{\link{gintervals.ls}},
#' \code{\link{gdataset.ls}}
#' @keywords ~intervals ~path ~database
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gintervals.dbs("annotations1")
#'
#' @export gintervals.dbs
gintervals.dbs <- function(intervals = NULL, dataframe = FALSE) {
    if (is.null(substitute(intervals))) {
        stop("Usage: gintervals.dbs(intervals)", call. = FALSE)
    }
    .gcheckroot()

    intervalsstr <- do.call(.gexpr2str, list(substitute(intervals)), envir = parent.frame())
    .gdb.resource_dbs_impl(intervalsstr, ".interv", "intervals", dataframe, gintervals.dbs)
}


#' Combines several sets of intervals
#'
#' Combines several sets of intervals into one set.
#'
#' This function combines several intervals sets into one set. It works in a
#' similar manner as 'rbind' yet it is faster. Also it supports intervals sets
#' that are stored in files including the big intervals sets.
#'
#' If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
#' set. If the format of the output intervals is set to be "big" (determined
#' implicitly based on the result size and options), the order of the resulted
#' intervals is altered as they are sorted by chromosome (or chromosomes pair -
#' for 2D).
#'
#' @param ... intervals sets to combine
#' @param intervals intervals set
#' @param intervals.set.out intervals set name where the function result is
#' optionally outputted
#' @return If 'intervals.set.out' is 'NULL' a data frame combining intervals
#' sets.
#' @seealso \code{\link{gintervals}}, \code{\link{gintervals.2d}},
#' \code{\link{gintervals.canonic}}
#' @keywords ~rbind
#' @examples
#' \dontshow{
#' options(gmultitasking = FALSE)
#' }
#'
#' gdb.init_examples()
#'
#' intervs1 <- gextract("sparse_track", gintervals(c(1, 2), 1000, 4000))
#' intervs2 <- gextract("sparse_track", gintervals(c(2, "X"), 2000, 5000))
#' gintervals.save("testintervs", intervs2)
#' gintervals.rbind(intervs1, "testintervs")
#' gintervals.rm("testintervs", force = TRUE)
#'
#' @export gintervals.rbind
gintervals.rbind <- function(..., intervals.set.out = NULL) {
    intervals <- list(...)
    if (!length(intervals)) {
        stop("Usage: gintervals.rbind([intervals]+, intervals.set.out = NULL)", call. = FALSE)
    }
    .gcheckroot()

    intervals.set.out <- do.call(.gexpr2str, list(substitute(intervals.set.out)), envir = parent.frame())

    res <- NULL
    if (any(unlist(lapply(intervals, function(intervals) {
        .gintervals.is_bigset(intervals)
    }))) || !is.null(intervals.set.out)) {
        if (is.null(intervals.set.out)) {
            FUN <- function(intervals, intervals.set.out, envir) {
                assign("res", c(get("res", envir = envir), intervals), envir = envir)
                .gverify_max_data_size(sum(unlist(lapply(get("res", envir), nrow))), arguments = "intervals.set.out")
                intervals[[1]]
            }

            # preserve the order of intervals inside the answer
            lapply(intervals, f <- function(intervals) {
                .gintervals.apply(gintervals.chrom_sizes(intervals), intervals, NULL, FUN, NULL, parent.frame(2))
            })
            if (!is.null(res)) {
                res <- do.call(.grbind, res)
            } # much faster than calling rbind incrementally in FUN
        } else {
            FUN <- function(intervals, intervals.set.out, envir) {
                intervals <- do.call(.grbind, intervals)
                intervals
            }

            # use for .gintervals.apply chromosomes from all intervals
            chroms <- NULL
            chroms <- lapply(intervals, gintervals.chrom_sizes)
            chroms <- do.call(rbind, chroms)
            if (.gintervals.is1d(intervals[[1]])) {
                chroms <- factor(chroms$chrom, levels(chroms$chrom))
                chroms <- unique(chroms)
                chroms <- sort(chroms)
                chroms <- data.frame(chrom = chroms)
            } else {
                chroms <- data.frame(chrom1 = chroms$chrom1, chrom2 = chroms$chrom2)
                chroms <- unique(chroms)
                chroms <- chroms[with(chroms, order(chrom1, chrom2)), ]
            }

            .gintervals.apply(chroms, intervals, intervals.set.out, FUN, intervals.set.out, environment())
        }
    } else {
        intervals <- lapply(intervals, .gintervals.load_ext)
        res <- do.call(.grbind, intervals) # much faster than calling rbind incrementally in FUN
    }

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


#' Deletes a named intervals set
#'
#' Deletes a named intervals set.
#'
#' This function deletes a named intervals set from the Genomic Database. By
#' default 'gintervals.rm' requires the user to interactively confirm the
#' deletion. Set 'force' to 'TRUE' to suppress the user prompt.
#'
#' @param intervals.set name of an intervals set
#' @param force if 'TRUE', suppresses user confirmation of a named intervals set
#' removal
#' @param db optional database path. When multiple databases are connected,
#' this specifies which database to delete the intervals set from. If NULL (the
#' default), the intervals set is deleted from the working database (GROOT).
#' @return None.
#' @seealso \code{\link{gintervals.save}}, \code{\link{gintervals.exists}},
#' \code{\link{gintervals.ls}}, \code{\link{gintervals}},
#' \code{\link{gintervals.2d}}, \code{\link{gtrack.rm}}
#' @keywords ~intervals
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' intervs <- gintervals(c(1, 2))
#' gintervals.save("testintervs", intervs)
#' gintervals.ls()
#' gintervals.rm("testintervs", force = TRUE)
#' gintervals.ls()
#'
#' @export gintervals.rm
gintervals.rm <- function(intervals.set = NULL, force = FALSE, db = NULL) {
    if (is.null(substitute(intervals.set))) {
        stop("Usage: gintervals.rm(intervals.set, force = FALSE, db = NULL)", call. = FALSE)
    }
    .gcheckroot()

    intervals.set <- do.call(.gexpr2str, list(substitute(intervals.set)), envir = parent.frame())

    # Determine the file path based on db parameter
    if (!is.null(db)) {
        db <- normalizePath(db, mustWork = FALSE)
        groot <- get("GROOT", envir = .misha)
        gdatasets <- get("GDATASETS", envir = .misha)
        if (is.null(gdatasets)) gdatasets <- character(0)
        groots <- c(groot, gdatasets)
        if (!(db %in% groots)) {
            stop(sprintf("Database %s is not connected", db), call. = FALSE)
        }
        fname <- file.path(db, "tracks", paste0(gsub("\\.", "/", intervals.set), ".interv"))
    } else {
        fname <- sprintf("%s.interv", paste(get("GWD", envir = .misha), gsub("\\.", "/", intervals.set), sep = "/"))
    }

    # check whether intervals.set appears among GINTERVS
    if (!(intervals.set %in% get("GINTERVS", envir = .misha))) {
        if (force) {
            unlink(fname, recursive = TRUE)
            .gdb.rm_intervals.set(intervals.set, db = db)
            return(invisible())
        }
        stop(sprintf("Intervals set %s does not exist", intervals.set), call. = FALSE)
    }

    if (!is.null(db) && !file.exists(fname)) {
        if (force) {
            return(invisible())
        }
        stop(sprintf("Intervals set %s does not exist in database %s", intervals.set, db), call. = FALSE)
    }

    answer <- "N"
    if (force) {
        answer <- "Y"
    } else {
        str <- sprintf("Are you sure you want to delete intervals set %s (Y/N)? ", intervals.set)
        message(str)
        answer <- toupper(readLines(n = 1))
    }

    if (answer == "Y" || answer == "YES") {
        # remove the intervals set
        unlink(fname, recursive = TRUE)

        if (file.exists(fname)) {
            message(sprintf("Failed to delete intervals set %s", intervals.set))
        } else {
            # refresh the list of GINTERVS, etc.
            .gdb.rm_intervals.set(intervals.set, db = db)
        }
    }
}


#' Creates a named intervals set
#'
#' Saves intervals to a named intervals set.
#'
#' This function saves 'intervals' as a named intervals set.
#'
#' @param intervals.set.out name of the new intervals set
#' @param intervals intervals to save
#' @return None.
#' @seealso \code{\link{gintervals.rm}}, \code{\link{gintervals.load}},
#' \code{\link{gintervals.exists}}, \code{\link{gintervals.ls}},
#' \code{\link{gintervals}}, \code{\link{gintervals.2d}}
#' @keywords ~intervals
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' intervs <- gintervals(c(1, 2))
#' gintervals.save("testintervs", intervs)
#' gintervals.ls()
#' gintervals.rm("testintervs", force = TRUE)
#'
#' @export gintervals.save
gintervals.save <- function(intervals.set.out = NULL, intervals = NULL) {
    if (is.null(substitute(intervals.set.out)) || is.null(intervals)) {
        stop("Usage: gintervals.save(intervals.set.out, intervals)", call. = FALSE)
    }
    .gcheckroot()

    intervals <- rescue_ALLGENOME(intervals, as.character(substitute(intervals)))

    intervals.set.out <- do.call(.gexpr2str, list(substitute(intervals.set.out)), envir = parent.frame())
    .gintervals.apply(gintervals.chrom_sizes(intervals), intervals, intervals.set.out, function(intervs, ...) {
        intervs[[1]]
    })
    retv <- NULL
}


#' Updates a named intervals set
#'
#' Updates a named intervals set.
#'
#' This function replaces all intervals of given chromosome (or chromosome
#' pair) within 'intervals.set' with 'intervals'. Chromosome is specified by
#' 'chrom' for 1D intervals set or 'chrom1', 'chrom2' for 2D intervals set.
#'
#' If 'intervals' is 'NULL' all intervals of given chromosome are removed from
#' 'intervals.set'.
#'
#' @param intervals.set name of an intervals set
#' @param intervals intervals or 'NULL'
#' @param chrom chromosome for 1D intervals set
#' @param chrom1 first chromosome for 2D intervals set
#' @param chrom2 second chromosome for 2D intervals set
#' @return None.
#' @seealso \code{\link{gintervals.save}}, \code{\link{gintervals.load}},
#' \code{\link{gintervals.exists}}, \code{\link{gintervals.ls}}
#' @keywords ~intervals
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' intervs <- gscreen(
#'     "sparse_track > 0.2",
#'     gintervals(c(1, 2), 0, 10000)
#' )
#' gintervals.save("testintervs", intervs)
#' gintervals.load("testintervs")
#' gintervals.update("testintervs", intervs[intervs$chrom == "chr2", ][1:5, ], chrom = 2)
#' gintervals.load("testintervs")
#' gintervals.update("testintervs", NULL, chrom = 2)
#' gintervals.load("testintervs")
#' gintervals.rm("testintervs", force = TRUE)
#'
#' @export gintervals.update
gintervals.update <- function(intervals.set = NULL, intervals = "", chrom = NULL, chrom1 = NULL, chrom2 = NULL) {
    if (is.null(substitute(intervals.set)) || identical(intervals, "")) {
        stop("Usage: gintervals.update(intervals.set, intervals, chrom = NULL, chrom1 = NULL, chrom2 = NULL)", call. = FALSE)
    }
    .gcheckroot()

    intervals <- rescue_ALLGENOME(intervals, as.character(substitute(intervals)))

    if (identical(intervals.set, intervals)) {
        return(retv <- NULL)
    }

    if (is.null(chrom) && is.null(chrom1) && is.null(chrom2)) {
        stop("Chromosome must be specified in chrom (for 2D intervals: chrom1, chrom2) parameter", call. = FALSE)
    }

    if (!is.null(chrom)) {
        chrom <- .gchroms(chrom)
        if (length(chrom) > 1) {
            stop("chrom parameter should mark only one chromosome")
        }
    }

    if (!is.null(chrom1)) {
        chrom1 <- .gchroms(chrom1)
        if (length(chrom1) > 1) {
            stop("chrom1 parameter should mark only one chromosome")
        }
    }

    if (!is.null(chrom2)) {
        chrom2 <- .gchroms(chrom2)
        if (length(chrom2) > 1) {
            stop("chrom2 parameter should mark only one chromosome")
        }
    }

    if (!is.null(chrom) && !is.null(chrom1)) {
        stop("Cannot use chrom and chrom1 parameters in the same call", call. = FALSE)
    }

    if (!is.null(chrom) && !is.null(chrom2)) {
        stop("Cannot use chrom and chrom2 parameters in the same call", call. = FALSE)
    }

    if (!is.character(intervals.set) || length(intervals.set) != 1) {
        stop("Invalid format of intervals.set parameter", call. = FALSE)
    }

    if (is.na(match(intervals.set, get("GINTERVS", envir = .misha)))) {
        stop(sprintf("Intervals set %s does not exist", intervals.set), call. = FALSE)
    }

    path <- gsub(".", "/", intervals.set, fixed = TRUE)
    path <- paste(path, ".interv", sep = "")
    fullpath <- paste(get("GWD", envir = .misha), path, sep = "/")

    if (!is.null(intervals)) {
        if (!is.null(chrom)) {
            intervals <- .gintervals.load_ext(intervals, chrom = chrom)
        } else {
            intervals <- .gintervals.load_ext(intervals, chrom1 = chrom1, chrom2 = chrom2)
        }
    }

    # big: update stats (including delete), save chrom (or delete), convert to small if needed
    if (.gintervals.is_bigset(intervals.set)) {
        is1d <- .gintervals.big.is1d(intervals.set)
        meta <- .gintervals.big.meta(intervals.set)
        stats <- meta$stats
        zeroline <- meta$zeroline

        if (!is.null(intervals) && !identical(sapply(intervals, "class"), sapply(zeroline, "class"))) {
            stop(sprintf("Cannot update intervals set %s: columns differ", intervals.set), call. = FALSE)
        }

        if (is1d) {
            if (is.null(chrom)) {
                stop("chrom parameter is not specified", call. = FALSE)
            }
            idx <- which(stats$chrom == chrom)
            if (length(idx) > 0) {
                stats <- stats[-idx, ]
            }
            if (!is.null(intervals) && nrow(intervals)) {
                stat <- .gcall("gintervals_stats", intervals, .misha_env())
                stats <- rbind(stats, data.frame(chrom = chrom, stat))
                stats <- stats[order(stats$chrom), ]
            }
            .gintervals.big.save(fullpath, intervals, chrom = chrom)
        } else {
            if (is.null(chrom1) || is.null(chrom2)) {
                stop("chrom1 and chrom2 parameters must be specified", call. = FALSE)
            }
            idx <- which(stats$chrom1 == chrom1 & stats$chrom2 == chrom2)
            if (length(idx) > 0) {
                stats <- stats[-idx, ]
            }
            if (!is.null(intervals) && nrow(intervals)) {
                stat <- .gcall("gintervals_stats", intervals, .misha_env())
                stats <- rbind(stats, data.frame(chrom1 = chrom1, chrom2 = chrom2, stat))
                stats <- stats[order(stats$chrom1, stats$chrom2), ]
            }
            .gintervals.big.save(fullpath, intervals, chrom1 = chrom1, chrom2 = chrom2)
        }

        if (nrow(stats) > 1) {
            rownames(stats) <- 1:nrow(stats)
        }
        .gintervals.big.save_meta(fullpath, stats, zeroline)

        if (!.gintervals.needs_bigset(intervals.set)) {
            .gintervals.big2small(intervals.set)
        }
    }

    # small: load all, update in place (including delete), save back, convert to big if needed
    else {
        tgt.intervals <- .gintervals.load_ext(intervals.set)
        is1d <- .gintervals.is1d(intervals.set)

        if (!is.null(intervals) && !identical(sapply(intervals, "class"), sapply(tgt.intervals, "class"))) {
            stop(sprintf("Cannot update intervals set %s: columns differ", intervals.set), call. = FALSE)
        }

        if (is1d) {
            if (is.null(chrom)) {
                stop("chrom parameter is not specified", call. = FALSE)
            }
            idx <- which(tgt.intervals$chrom == chrom)
            if (length(idx) > 0) {
                tgt.intervals <- tgt.intervals[-idx, ]
            }
            if (!is.null(intervals) && nrow(intervals)) {
                tgt.intervals <- .grbind(tgt.intervals, intervals)
                tgt.intervals <- tgt.intervals[order(tgt.intervals$chrom), ]
            }
        } else {
            if (is.null(chrom1) || is.null(chrom2)) {
                stop("chrom1 and chrom2 parameters must be specified", call. = FALSE)
            }
            idx <- which(tgt.intervals$chrom1 == chrom1 & tgt.intervals$chrom2 == chrom2)
            if (length(idx) > 0) {
                tgt.intervals <- tgt.intervals[-idx, ]
            }
            if (!is.null(intervals) && nrow(intervals)) {
                tgt.intervals <- .grbind(tgt.intervals, intervals)
                tgt.intervals <- tgt.intervals[order(tgt.intervals$chrom1, tgt.intervals$chrom2), ]
            }
        }
        if (.gintervals.needs_bigset(tgt.intervals)) {
            .gintervals.small2big(intervals.set, tgt.intervals)
        } else {
            .gintervals.save_file(fullpath, tgt.intervals)
        }
    }

    retv <- 0 # suppress return value
}
