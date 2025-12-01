# Interval loading and saving infrastructure

.gintervals.apply <- function(chroms, intervals, intervals.set.out, FUN, ...) {
    if (!is.null(intervals.set.out)) {
        fullpath <- .gintervals.check_new_set(intervals.set.out)
    }

    if (is.data.frame(intervals)) {
        intervals <- list(intervals)
    }

    # sort chroms
    chroms$size <- NULL
    if ("chrom" %in% colnames(chroms)) {
        chroms <- data.frame(chrom = chroms[with(chroms, order(chrom)), ])
    } else {
        chroms <- chroms[with(chroms, order(chrom1, chrom2)), ]
    }


    # let's assume that if any of the input intervals sets is big then intervals.set.out should be big as well
    if (any(unlist(lapply(intervals, function(intervals) {
        .gintervals.is_bigset(intervals) || .gintervals.needs_bigset(intervals)
    })))) {
        stats <- NULL
        zeroline <- NULL
        success <- FALSE
        t <- Sys.time()
        progress.percentage <- -1
        tryCatch(
            {
                # if any of the source intervals sets is big then create the output intervals set big too
                if (!is.null(intervals.set.out)) {
                    dir.create(fullpath, recursive = TRUE, mode = "0777")
                }

                if (.gintervals.is1d(intervals[[1]])) {
                    mapply(function(chrom) {
                        loaded_intervals <- lapply(intervals, function(intervals) {
                            .gintervals.load_ext(intervals, chrom = chrom)
                        })
                        res <- do.call(FUN, list(loaded_intervals, ...))
                        if (!is.null(intervals.set.out) && !is.null(res) && nrow(res) > 0) {
                            zeroline <<- res[0, ]
                            .gintervals.big.save(fullpath, res, chrom = chrom)
                            stat <- .gcall("gintervals_stats", res, .misha_env())
                            stats <<- rbind(stats, data.frame(chrom = chrom, stat))
                        }
                        if (as.integer(difftime(Sys.time(), t, units = "secs")) > 3) {
                            t <<- Sys.time()
                            percentage <- as.integer(100 * match(chrom, chroms$chrom) / nrow(chroms))
                            if (percentage < 100 && progress.percentage != percentage) {
                                message(sprintf("%d%%...", percentage), appendLF = FALSE)
                                progress.percentage <<- percentage
                            }
                        }
                    }, chroms$chrom)
                } else {
                    mapply(function(chrom1, chrom2) {
                        loaded_intervals <- lapply(intervals, function(intervals) {
                            .gintervals.load_ext(intervals, chrom1 = chrom1, chrom2 = chrom2)
                        })
                        res <- do.call(FUN, list(loaded_intervals, ...))
                        if (!is.null(intervals.set.out) && !is.null(res) && nrow(res) > 0) {
                            zeroline <<- res[0, ]
                            .gintervals.big.save(fullpath, res, chrom1 = chrom1, chrom2 = chrom2)
                            stat <- .gcall("gintervals_stats", res, .misha_env())
                            stats <<- rbind(stats, data.frame(chrom1 = chrom1, chrom2 = chrom2, stat))
                        }
                        if (as.integer(difftime(Sys.time(), t, units = "secs")) > 3) {
                            t <<- Sys.time()
                            percentage <- as.integer(100 * which(chroms$chrom1 == chrom1 & chroms$chrom2 == chrom2) / nrow(chroms))
                            if (percentage < 100 && progress.percentage != percentage) {
                                message(sprintf("%d%%...", percentage), appendLF = FALSE)
                                progress.percentage <<- percentage
                            }
                        }
                    }, chroms$chrom1, chroms$chrom2)
                }

                if (!is.null(intervals.set.out)) {
                    if (is.null(stats)) {
                        return(retv <- NULL)
                    }
                    .gintervals.big.save_meta(fullpath, stats, zeroline)
                }

                if (progress.percentage >= 0) {
                    message("100%")
                }

                success <- TRUE

                # check whether the output intervals set needs to remain in big format
                if (!is.null(intervals.set.out) && !.gintervals.needs_bigset(intervals.set.out)) {
                    .gintervals.big2small(intervals.set.out)
                }

                # If indexed format is enabled and output is bigset, convert to indexed format
                if (!is.null(intervals.set.out) && getOption("gmulticontig.indexed_format", FALSE) && .gintervals.is_bigset(intervals.set.out)) {
                    if (.gintervals.is1d(intervals.set.out)) {
                        gintervals.convert_to_indexed(intervals.set.out, remove.old = TRUE)
                    } else {
                        gintervals.2d.convert_to_indexed(intervals.set.out, remove.old = TRUE)
                    }
                }
            },
            finally = {
                if (!success && !is.null(intervals.set.out)) {
                    unlink(fullpath, recursive = TRUE)
                }
            }
        )
    } else {
        loaded_intervals <- lapply(intervals, .gintervals.load_ext)
        res <- do.call(FUN, list(loaded_intervals, ...))
        if (!is.null(intervals.set.out) && !is.null(res) && nrow(res) > 0) {
            if (.gintervals.is1d(res)) {
                res <- res[order(res$chrom), ]
            } else {
                res <- res[order(res$chrom1, res$chrom2), ]
            }
            if (.gintervals.needs_bigset(res)) {
                .gintervals.small2big(intervals.set.out, res)

                # If indexed format is enabled, convert to indexed format
                if (getOption("gmulticontig.indexed_format", FALSE)) {
                    if (.gintervals.is1d(res)) {
                        gintervals.convert_to_indexed(intervals.set.out, remove.old = TRUE)
                    } else {
                        gintervals.2d.convert_to_indexed(intervals.set.out, remove.old = TRUE)
                    }
                }
            } else {
                .gintervals.save_file(fullpath, res)
            }
        } else {
            return(NULL)
        }
    }

    # refresh the list of GINTERVS, etc.
    if (!is.null(intervals.set.out)) {
        .gdb.add_intervals.set(intervals.set.out)
    }
}

.gintervals.check_new_set <- function(intervals.set) {
    if (!is.na(match(intervals.set, get("GINTERVS", envir = .misha)))) {
        stop(sprintf("Intervals set %s already exists", intervals.set), call. = FALSE)
    }

    if (!length(grep("^[A-Za-z][\\w.]*$", intervals.set, perl = TRUE))) {
        stop("Invalid interval name %s. Only alphanumeric characters and _ are allowed in the name.")
    }

    path <- gsub(".", "/", intervals.set, fixed = TRUE)
    path <- paste(path, ".interv", sep = "")
    fullpath <- paste(get("GWD", envir = .misha), path, sep = "/")
    dir <- dirname(path)
    fulldir <- paste(get("GWD", envir = .misha), dir, sep = "/")

    if (!file.exists(fulldir)) {
        stop(sprintf("Directory %s does not exist", dir), call. = FALSE)
    }

    if (file.exists(fullpath)) {
        stop(sprintf("File %s already exists", path), call. = FALSE)
    }

    if (!is.na(match(intervals.set, get("GTRACKS", envir = .misha)))) {
        stop(sprintf("Track %s already exists", intervals.set), call. = FALSE)
    }

    if (!is.na(match(intervals.set, gvtrack.ls()))) {
        stop(sprintf("Virtual track %s already exists", intervals.set), call. = FALSE)
    }

    if (.ggetOption(".gautocompletion", FALSE) && exists(intervals.set, envir = .misha)) {
        stop(sprintf("Variable \"%s\" shadows the name of the new intervals set.\nPlease remove this variable from the environment or switch off autocompletion mode.", intervals.set), call. = FALSE)
    }

    fullpath
}

.gintervals.load_ext <- function(intervals.set = NULL, chrom = NULL, chrom1 = NULL, chrom2 = NULL, progress = FALSE) {
    if (is.null(intervals.set)) {
        stop("Usage: gintervals.load(intervals.set, chrom = NULL, chrom1 = NULL, chrom2 = NULL)", call. = FALSE)
    }
    .gcheckroot()

    if (is.character(intervals.set) && length(intervals.set) == 1 && is.na(match(intervals.set, get("GINTERVS", envir = .misha))) && is.na(match(intervals.set, get("GTRACKS", envir = .misha)))) {
        stop(sprintf("Intervals set %s does not exist", intervals.set), call. = FALSE)
    }

    .gintervals.load(intervals.set, chrom, chrom1, chrom2, progress)
}

.gintervals.load <- function(intervals.set = NULL, chrom = NULL, chrom1 = NULL, chrom2 = NULL, progress = FALSE) {
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

    if (is.character(intervals.set) && length(intervals.set) != 1 || !is.character(intervals.set) && !.gintervals.is1d(intervals.set) && !.gintervals.is2d(intervals.set)) {
        stop("Invalid format of intervals", call. = FALSE)
    }

    res <- NULL
    if (.gintervals.is_bigset(intervals.set)) {
        meta <- .gintervals.big.meta(intervals.set)
        zeroline <- meta$zeroline
        t <- Sys.time()
        progress.percentage <- -1
        if (.gintervals.big.is1d(intervals.set)) {
            if (!is.null(chrom1) || !is.null(chrom2)) {
                stop(sprintf("%s is a 1D big intervals set.\nchrom1 or chrom2 parameters can be applied only to 2D intervals.", intervals.set), call. = FALSE)
            }

            if (!is.null(chrom)) {
                # Convert to character for comparison to handle different factor levels
                meta$stats <- meta$stats[as.character(meta$stats$chrom) == as.character(chrom), ]
            }

            if (!.gintervals.loadable(intervals.set, chrom = chrom)) {
                if (is.null(chrom)) {
                    stop(sprintf(
                        "Cannot load a big intervals set %s: its size (%d) exceeds the limit (%d) controlled by gmax.data.size option.\nFor big intervals sets only one chromosome pair can be loaded at a time.",
                        intervals.set, sum(meta$stats$size), .ggetOption("gmax.data.size", 10000000)
                    ), call. = FALSE)
                } else {
                    stop(sprintf(
                        "Cannot load chromosome %s of an intervals set %s: its size (%d) exceeds the limit (%d) controlled by gmax.data.size option.",
                        chrom, intervals.set, sum(meta$stats$size), .ggetOption("gmax.data.size", 10000000)
                    ), call. = FALSE)
                }
            }
            if (nrow(meta$stats) > 1) {
                res <- list(zeroline)
                lapply(
                    meta$stats$chrom,
                    function(chrom) {
                        loaded_intervs <- .gintervals.load_file(intervals.set, chrom = chrom)
                        if (!identical(sapply(loaded_intervs, "class"), sapply(zeroline, "class"))) {
                            stop(sprintf("Intervals set %s, chrom %s: invalid columns definition", intervals.set, chrom), call. = FALSE)
                        }
                        res <<- c(res, list(loaded_intervs))
                        if (as.integer(difftime(Sys.time(), t, units = "secs")) > 3) {
                            t <<- Sys.time()
                            percentage <- as.integer(100 * match(chrom, meta$stats$chrom) / length(meta$stats$chrom))
                            if (progress && percentage < 100 && progress.percentage != percentage) {
                                message(sprintf("%d%%...", percentage), appendLF = FALSE)
                                progress.percentage <<- percentage
                            }
                        }
                    }
                )
                res <- do.call(.grbind, res) # much faster than calling rbind incrementally in mapply
            } else if (nrow(meta$stats) == 1) {
                res <- .gintervals.load_file(intervals.set, chrom = meta$stat$chrom[1])
                if (!identical(sapply(res, "class"), sapply(zeroline, "class"))) {
                    stop(sprintf("Intervals set %s, chrom %s: invalid columns definition", intervals.set, meta$stat$chrom[1]), call. = FALSE)
                }
            } else {
                res <- meta$zeroline
            }
        } else {
            if (!is.null(chrom)) {
                stop(sprintf("%s is a 2D big intervals set.\nchrom parameter can be applied only to 1D intervals.", intervals.set), call. = FALSE)
            }

            if (!is.null(chrom1)) {
                # Convert to character for comparison to handle different factor levels
                meta$stats <- meta$stats[as.character(meta$stats$chrom1) == as.character(chrom1), ]
            }
            if (!is.null(chrom2)) {
                # Convert to character for comparison to handle different factor levels
                meta$stats <- meta$stats[as.character(meta$stats$chrom2) == as.character(chrom2), ]
            }

            if (!.gintervals.loadable(intervals.set, chrom1 = chrom1, chrom2 = chrom2)) {
                if (!is.null(chrom1) && !is.null(chrom2)) {
                    stop(sprintf(
                        "Cannot load chromosome pair (%s, %s) of an intervals set %s: its size (%d) exceeds the limit (%d) controlled by gmax.data.size option.",
                        chrom1, chrom2, intervals.set, sum(meta$stats$size), .ggetOption("gmax.data.size", 10000000)
                    ), call. = FALSE)
                } else if (!is.null(chrom1)) {
                    stop(sprintf(
                        "Cannot load chromosome %s of an intervals set %s: its size (%d) exceeds the limit (%d) controlled by gmax.data.size option.",
                        chrom1, intervals.set, sum(meta$stats$size), .ggetOption("gmax.data.size", 10000000)
                    ), call. = FALSE)
                } else if (!is.null(chrom2)) {
                    stop(sprintf(
                        "Cannot load chromosome %s of an intervals set %s: its size (%d) exceeds the limit (%d) controlled by gmax.data.size option.",
                        chrom2, intervals.set, sum(meta$stats$size), .ggetOption("gmax.data.size", 10000000)
                    ), call. = FALSE)
                } else {
                    stop(sprintf(
                        "Cannot load a big intervals set %s: its size (%d) exceeds the limit (%d) controlled by gmax.data.size option.\nFor big intervals sets only one chromosome pair can be loaded at a time.",
                        intervals.set, sum(meta$stats$size), .ggetOption("gmax.data.size", 10000000)
                    ), call. = FALSE)
                }
            }

            if (nrow(meta$stats) > 1) {
                res <- list(zeroline)
                mapply(
                    function(chrom1, chrom2) {
                        loaded_intervs <- .gintervals.load_file(intervals.set, chrom1 = chrom1, chrom2 = chrom2)
                        if (!identical(sapply(loaded_intervs, "class"), sapply(zeroline, "class"))) {
                            stop(sprintf("Interval set %s, chrom1 %s, chrom2 %s: invalid columns definition", intervals.set, chrom1, chrom2), call. = FALSE)
                        }
                        res <<- c(res, list(loaded_intervs))
                        if (as.integer(difftime(Sys.time(), t, units = "secs")) > 3) {
                            t <<- Sys.time()
                            percentage <- as.integer(100 * which(meta$stats$chrom1 == chrom1 & meta$stats$chrom2 == chrom2) / nrow(meta$stats))
                            if (progress && percentage < 100 && progress.percentage != percentage) {
                                message(sprintf("%d%%...", percentage), appendLF = FALSE)
                                progress.percentage <<- percentage
                            }
                        }
                    },
                    meta$stats$chrom1, meta$stats$chrom2
                )
                res <- do.call(.grbind, res) # much faster than calling rbind incrementally in mapply
            } else if (nrow(meta$stats) == 1) {
                res <- .gintervals.load_file(intervals.set, chrom1 = meta$stat$chrom1[1], chrom2 = meta$stat$chrom2[1])
                if (!identical(sapply(res, "class"), sapply(zeroline, "class"))) {
                    stop(sprintf("Interval set %s, chrom1 %s, chrom2 %s: invalid columns definition", intervals.set, meta$stat$chrom1[1], meta$stat$chrom2[1]), call. = FALSE)
                }
            } else {
                res <- meta$zeroline
            }
        }

        if (progress.percentage >= 0) {
            message("100%")
        }
    } else {
        if (is.character(intervals.set) && length(intervals.set) == 1) {
            res <- .gintervals.load_file(intervals.set)
        } else {
            res <- intervals.set
        }
        if (!is.null(res)) {
            if (!.gintervals.is1d(res) && !is.null(chrom)) {
                stop("chrom parameter can be applied only to 1D intervals", call. = FALSE)
            }

            if (!.gintervals.is2d(res) && (!is.null(chrom1) || !is.null(chrom2))) {
                stop("chrom1 or chrom2 parameters can be applied only to 2D intervals", call. = FALSE)
            }

            if (nrow(res) > 0) {
                if (!is.null(chrom)) {
                    res <- res[res$chrom == chrom, ]
                    if (nrow(res)) {
                        rownames(res) <- 1:nrow(res)
                    }
                }
                if (!is.null(chrom1)) {
                    res <- res[res$chrom1 == chrom1, ]
                    if (nrow(res)) {
                        rownames(res) <- 1:nrow(res)
                    }
                }
                if (!is.null(chrom2)) {
                    res <- res[res$chrom2 == chrom2, ]
                    if (nrow(res)) {
                        rownames(res) <- 1:nrow(res)
                    }
                }
            }
        }
    }
    res
}

.gintervals.load_file <- function(intervals.set, chrom = NULL, chrom1 = NULL, chrom2 = NULL) {
    intervfname <- sprintf("%s.interv", paste(get("GWD", envir = .misha), gsub("\\.", "/", intervals.set), sep = "/"))
    if (!is.null(chrom)) {
        chrom <- .gchroms(chrom)
        intervfname <- sprintf("%s/%s", intervfname, chrom)
    } else if (!is.null(chrom1) && !is.null(chrom2)) {
        chrom1 <- .gchroms(chrom1)
        chrom2 <- .gchroms(chrom2)
        intervfname <- sprintf("%s/%s-%s", intervfname, chrom1, chrom2)
    }

    if (intervals.set %in% get("GTRACKS", envir = .misha)) {
        .gcall("gtrack_intervals_load", intervals.set, chrom, chrom1, chrom2, .misha_env())
    } else {
        if (file.exists(intervfname) || (is.null(chrom) && is.null(chrom1) && is.null(chrom2))) {
            f <- file(intervfname, "rb")
            intervals.set <- unserialize(f)
            close(f)
            intervals.set
        } else {
            if (.gintervals.is_bigset(intervals.set)) {
                # For big intervals sets, check if we're loading a specific chromosome
                # and use the C++ function that handles both per-chromosome and indexed formats
                if (!is.null(chrom)) {
                    return(.gcall("gbigintervs_load_chrom", intervals.set, chrom, .misha_env()))
                } else if (!is.null(chrom1) && !is.null(chrom2)) {
                    return(.gcall("gbigintervs_load_chrom2d", intervals.set, chrom1, chrom2, .misha_env()))
                } else {
                    .gintervals.big.meta(intervals.set)$zeroline
                }
            } else {
                stop(sprintf("File %s does not exist", intervfname), call. = FALSE)
            }
        }
    }
}

.gintervals.save_file <- function(filename, intervs) {
    if (nrow(intervs)) {
        if (.gintervals.is1d(intervs)) {
            intervs$chrom <- as.factor(intervs$chrom)
            point.intervs <- intervs$start == intervs$end
            intervs[point.intervs, ]$end <- intervs[point.intervs, ]$end + 1
        } else {
            intervs$chrom1 <- as.factor(intervs$chrom1)
            intervs$chrom2 <- as.factor(intervs$chrom2)
            point.intervs <- intervs$start1 == intervs$end1
            intervs[point.intervs, ]$end1 <- intervs[point.intervs, ]$end1 + 1
            point.intervs <- intervs$start2 == intervs$end2
            intervs[point.intervs, ]$end2 <- intervs[point.intervs, ]$end2 + 1
        }
    }

    f <- file(filename, "wb")
    serialize(intervs, f)
    close(f)
    nrow(intervs)
}

# Internal helper: save a computed intervals data.frame to a new intervals set, or return it
.gintervals.save_set_or_return <- function(result, intervals.set.out) {
    intervals.set.out <- do.call(.gexpr2str, list(substitute(intervals.set.out)), envir = parent.frame())

    if (is.null(intervals.set.out)) {
        return(result)
    }

    fullpath <- .gintervals.check_new_set(intervals.set.out)

    # Nothing to save
    if (is.null(result) || nrow(result) == 0) {
        return(NULL)
    }

    success <- FALSE
    tryCatch(
        {
            # Ensure deterministic ordering
            if (.gintervals.is1d(result)) {
                result <- result[order(result$chrom), ]
            } else {
                result <- result[order(result$chrom1, result$chrom2), ]
            }

            # Save as big or small set
            if (.gintervals.needs_bigset(result)) {
                .gintervals.small2big(intervals.set.out, result)
            } else {
                .gintervals.save_file(fullpath, result)
            }
            success <- TRUE
        },
        finally = {
            if (!success) {
                unlink(fullpath, recursive = TRUE)
            }
        }
    )

    # Register the new set and suppress return value
    .gdb.add_intervals.set(intervals.set.out)
    retv <- 0
}
