# Big intervals set handling

.gintervals.is_bigset <- function(intervals.set, err_if_non_exist = TRUE) {
    if (is.character(intervals.set) & length(intervals.set) == 1) {
        if (intervals.set %in% get("GTRACKS", envir = .misha)) {
            intervfname <- sprintf("%s.track", paste(get("GWD", envir = .misha), gsub("\\.", "/", intervals.set), sep = "/"))
        } else {
            intervfname <- sprintf("%s.interv", paste(get("GWD", envir = .misha), gsub("\\.", "/", intervals.set), sep = "/"))
        }
        if (file.exists(intervfname)) {
            if (file.info(intervfname)$isdir) {
                return(TRUE)
            }
        } else if (err_if_non_exist) {
            stop(sprintf("Intervals set %s does not exist", intervals.set), call. = FALSE)
        }
    }
    FALSE
}

.gintervals.needs_bigset <- function(intervals = NULL, size = NULL) {
    if (!is.null(intervals)) {
        chromsizes <- gintervals.chrom_sizes(intervals)
        nrow(chromsizes) && sum(chromsizes$size) > min(.ggetOption("gmax.data.size", 10000000), .ggetOption("gbig.intervals.size", 1000000))
    } else {
        size > min(.ggetOption("gmax.data.size", 10000000), .ggetOption("gbig.intervals.size", 1000000))
    }
}

.gintervals.loadable <- function(intervals = NULL, size = NULL, chrom = NULL, chrom1 = NULL, chrom2 = NULL) {
    if (!is.null(intervals)) {
        chromsizes <- gintervals.chrom_sizes(intervals)
        if (nrow(chromsizes)) {
            if (!is.null(chrom)) {
                chromsizes <- chromsizes[chromsizes$chrom == chrom, ]
            }
            if (!is.null(chrom1)) {
                chromsizes <- chromsizes[chromsizes$chrom1 == chrom1, ]
            }
            if (!is.null(chrom2)) {
                chromsizes <- chromsizes[chromsizes$chrom2 == chrom2, ]
            }
        }
        !nrow(chromsizes) || sum(chromsizes$size) <= .ggetOption("gmax.data.size", 10000000)
    } else {
        size <= .ggetOption("gmax.data.size", 10000000)
    }
}

.gintervals.big.is1d <- function(intervals.set) {
    "chrom" %in% colnames(.gintervals.big.meta(intervals.set)$stats)
}

.gintervals.big.is2d <- function(intervals.set) {
    "chrom1" %in% colnames(.gintervals.big.meta(intervals.set)$stats)
}

.gintervals.big2small <- function(intervals.set) {
    # We assume that writing the intervals might be a lengthy process.
    # During this time the process might get interrupted leaving the intervals set in incomplete state.
    # Even though it's not fully transaction-safe, we prefer to create a temporary file and then move it hoping it's fast enough.
    path <- gsub(".", "/", intervals.set, fixed = TRUE)
    path <- paste(get("GWD", envir = .misha), path, sep = "/")
    path <- paste(path, ".interv", sep = "")

    intervals <- .gintervals.load(intervals.set)
    tmpfilename <- tempfile(".", dirname(path)) # tmpdir = the parent directory of intervals set, otherwise rename might not work
    # (tmpdir might be at another file system)
    file.rename(path, tmpfilename)
    success <- FALSE
    tryCatch(
        {
            .gintervals.save_file(path, intervals)
            success <- TRUE
        },
        finally = {
            if (!success) {
                unlink(path, recursive = TRUE)
                file.rename(tmpfilename, path)
            }
            unlink(tmpfilename, recursive = TRUE)
        }
    )
}

.gintervals.small2big <- function(intervals.set, intervals = NULL) {
    # We assume that writing the intervals might be a lengthy process.
    # During this time the process might get interrupted leaving the intervals set in incomplete state.
    # Even though it's not fully transaction-safe, we prefer to create a temporary file and then move it hoping it's fast enough.
    if (is.null(intervals)) {
        intervals <- .gintervals.load(intervals.set)
    }

    path <- gsub(".", "/", intervals.set, fixed = TRUE)
    path <- paste(get("GWD", envir = .misha), path, sep = "/")
    path <- paste(path, ".interv", sep = "")

    tmpfilename <- tempfile(".", dirname(path)) # tmpdir = the parent directory of intervals set, otherwise rename might not work
    # (tmpdir might be at another file system)
    file.rename(path, tmpfilename)
    gintervs <- get("GINTERVS", envir = .misha)
    assign("GINTERVS", gintervs[gintervs != intervals.set], envir = .misha)
    success <- FALSE
    tryCatch(
        {
            .gcall_noninteractive(gintervals.save, intervals.set, intervals)
            success <- TRUE
        },
        finally = {
            if (!success) {
                unlink(path, recursive = TRUE)
                file.rename(tmpfilename, path)
                gintervs <- c(gintervs, intervals.set)
                gintervs <- unique(gintervs)
                sort(gintervs)
                assign("GINTERVS", gintervs, envir = .misha)
            }
            unlink(tmpfilename, recursive = TRUE)
        }
    )
}

.gintervals.big.meta <- function(intervals.set) {
    metafname <- ""
    if (intervals.set %in% get("GTRACKS", envir = .misha)) {
        metafname <- sprintf("%s.track/.meta", paste(get("GWD", envir = .misha), gsub("\\.", "/", intervals.set), sep = "/"))
        if (!file.exists(metafname)) {
            .gcall("gtrack_create_meta", intervals.set, .misha_env())
        }
    } else {
        metafname <- sprintf("%s.interv/.meta", paste(get("GWD", envir = .misha), gsub("\\.", "/", intervals.set), sep = "/"))
    }
    f <- file(metafname, "rb")
    res <- unserialize(f)
    close(f)
    res
}

.gintervals.big.save_meta <- function(path, stats, zeroline) {
    f <- file(sprintf("%s/.meta", path), "wb")
    serialize(list(stats = stats, zeroline = zeroline), f)
    close(f)
}

.gintervals.big.save <- function(path, intervs, chrom = NULL, chrom1 = NULL, chrom2 = NULL) {
    if (!is.null(chrom)) {
        filename <- sprintf("%s/%s", path, chrom)
    } else {
        filename <- sprintf("%s/%s-%s", path, chrom1, chrom2)
    }

    if (is.null(intervs) || nrow(intervs) == 0) {
        unlink(filename, recursive = TRUE)
    } else {
        if (!is.null(chrom)) {
            intervs <- intervs[intervs$chrom == chrom, ]
        } else {
            intervs <- intervs[intervs$chrom1 == chrom1 & intervs$chrom2 == chrom2, ]
        }
        .gintervals.save_file(filename, intervs)
    }
}
