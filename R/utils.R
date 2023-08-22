.gcall <- function(...) {
    tryCatch(
        {
            res <- .Call(...)
        },
        interrupt = function(interrupt) {
            stop("Command interrupted!", call. = FALSE)
        }
    )
    res
}

.misha_env <- function() {
    e <- new.env(parent = parent.frame(2))
    assign(".misha", .misha, envir = e)
    return(e)
}

.gcall_noninteractive <- function(FUN, ...) {
    .ginteractive <- .ggetOption(".ginteractive")
    tryCatch(
        {
            options(.ginteractive = F)
            do.call(FUN, list(...))
        },
        finally = {
            options(.ginteractive = .ginteractive)
        }
    )
}

.gdefine_autocompletion_vars <- function(tracks, intervals, vtracks, interactive) {
    for (track in tracks) {
        if (exists(track)) {
            stop(sprintf("Variable \"%s\" shadows the name of identically named track.\nPlease remove this variable from the environment or switch off autocompletion mode.", track), call. = F)
        }
    }

    for (interval in intervals) {
        if (exists(interval)) {
            stop(sprintf("Variable \"%s\" shadows the name of identically named intervals set.\nPlease remove this variable from the environment or switch off autocompletion mode.", interval), call. = F)
        }
    }

    for (vtrack in vtracks) {
        if (exists(vtrack)) {
            stop(sprintf("Variable \"%s\" shadows the name of identically named virtual track.\nPlease remove this variable from the environment or switch off autocompletion mode.", vtrack), call. = F)
        }
    }

    if (interactive) {
        # set track to NULL otherwise evaluation of track expression pmin(track, 2) will produce a string "2"
        for (track in tracks) {
            assign(track, NULL, envir = .GlobalEnv)
        }

        for (vtrack in vtracks) {
            assign(vtrack, NULL, envir = .GlobalEnv)
        }
    } else {
        for (track in tracks) {
            assign(track, track, envir = .GlobalEnv)
        }

        for (vtrack in vtracks) {
            assign(vtrack, vtrack, envir = .GlobalEnv)
        }
    }

    for (interval in intervals) {
        assign(interval, interval, envir = .GlobalEnv)
    }
}

.gundefine_autocompletion_vars <- function() {
    options(warn = -1) # disable warnings: some variables might be removed already by the user

    if (exists("GTRACKS", envir = .misha)) {
        remove(list = get("GTRACKS", envir = .misha), envir = .GlobalEnv)
    }

    if (exists("GINTERVS", envir = .misha)) {
        remove(list = get("GINTERVS", envir = .misha), envir = .GlobalEnv)
    }


    vtracks <- gvtrack.ls()
    if (!is.null(vtracks)) {
        remove(list = vtracks, envir = .GlobalEnv)
    }

    options(warn = -1) # restore warnings
}

.ggetOption <- function(x, default = NULL) {
    if (missing(default)) {
        return(options(x)[[1L]])
    }
    if (x %in% names(options())) {
        options(x)[[1L]]
    } else {
        default
    }
}

.gexpr2str <- function(x) {
    if (.ggetOption(".ginteractive", FALSE)) {
        if (is.null(substitute(x))) {
            NULL
        } else {
            str <- deparse(substitute(x), width.cutoff = 500)[1]
            gsub("^\"(.*)\"$", "\\1", str, perl = T)
        }
    } else {
        eval.parent(x)
    }
}

.giterator <- function(iterator) {
    if (typeof(iterator) == "integer" || typeof(iterator) == "double") {
        return(iterator)
    }

    iterator.str <- do.call(.gexpr2str, list(substitute(iterator)), envir = parent.frame())

    if (typeof(iterator.str) == "character") {
        if (!is.na(match(iterator.str, get("GTRACKS", envir = .misha))) || !is.na(match(iterator.str, get("GINTERVS", envir = .misha)))) {
            return(iterator.str)
        }
    }

    iterator
}

.grbind <- function(...) {
    objs <- list(...)

    zerolines <- lapply(objs, function(obj) {
        obj[0, ]
    })

    diffs <- sapply(zerolines, FUN = attr.all.equal, zerolines[[1]])
    if (!all(sapply(diffs, FUN = is.null))) {
        stop("Cannot rbind objects: columns differ", call. = F)
    }

    .gcall("grbind", objs, .misha_env())
}

.gverify_max_data_size <- function(size, data_name = "Result", arguments = NULL) {
    max.data.size <- .ggetOption("gmax.data.size", 10000000)

    if (size > max.data.size) {
        if (is.null(arguments)) {
            stop(sprintf(
                paste("%s size exceeded the maximal allowed (%d).",
                    "Note: the maximum data size is controlled via gmax.data.size option (see options, getOptions).",
                    sep = "\n"
                ),
                data_name, max.data.size
            ), call. = F)
        } else {
            stop(sprintf(
                paste("%s size exceeded the maximal allowed (%d).",
                    "Consider saving the result in a file (use %s argument).",
                    "Note: the maximum data size is controlled via gmax.data.size option (see options, getOptions).",
                    sep = "\n"
                ),
                data_name, max.data.size, paste(arguments, collapse = " or ")
            ), call. = F)
        }
    }
}




#' Sets input mode and auto-completion for track expressions, track names,
#' virtual tracks and interval sets
#'
#' Sets input mode and auto-completion for track expressions, track names,
#' virtual tracks and interval sets.
#'
#' This function enables / disables auto-completion for track names, virtual
#' tracks and interval sets. It also controls whether these objects together
#' with track expressions should be passed as strings or "as is" to the various
#' package functions.
#'
#' If 'autocompletion' is 'TRUE' all the track names, virtual tracks and
#' intervals sets are defined as R variables (auxiliary variables) which allows
#' them to be auto-completed by TAB key. The values of these variables are
#' meaningless for the user and they should not be altered.
#'
#' If 'interactive' is 'TRUE' track names, virtual tracks, interval sets and
#' track expressions are passed to the package functions "as is", i.e.
#' unquoted.
#'
#' 'autocompletion' is required to be switched on if 'interactive' mode is on
#' too.
#'
#' Please beware of the consequences of using interactive mode as it creates a
#' bunch of new variables in R environment. Though collision with the existing
#' variables is checked at the time of the call to 'gset_input_mode', yet
#' nothing prevents the user to modify the value of the auxiliary variables
#' later. This might cause unexpected behaviour in some of the package
#' functions. Also the auxiliary variables are automatically undefined once the
#' interactive mode is switched off. User who mistakenly uses auxiliary
#' variables to store the data might therefore accidentially loose it.
#'
#' @param autocompletion if 'TRUE' auto-completion is switched on
#' @param interactive if 'TRUE' interactive mode is switched on
#'
#' @return None.
#' @keywords ~autocompletion ~interactive ~quoted
#' @examples
#'
#' gdb.init_examples()
#' gset_input_mode(interactive = FALSE)
#' gsummary("dense_track + 10")
#' gset_input_mode(autocompletion = TRUE, interactive = TRUE)
#' gsummary(dense_track + 10)
#'
#' # roll-back to default input mode
#' gset_input_mode()
#'
#' @export gset_input_mode
gset_input_mode <- function(autocompletion = FALSE, interactive = FALSE) {
    if (!is.logical(autocompletion) || !is.logical(interactive) || length(autocompletion) != 1 || length(interactive) != 1) {
        stop("Usage: gset_input_mode(autocompletion = FALSE, interactive = FALSE)", call. = F)
    }

    if (!autocompletion && interactive) {
        stop("Autocompletion must be switched on in interactive mode", call. = F)
    }

    if (.ggetOption(".gautocompletion") != autocompletion) {
        if (autocompletion) {
            tracks <- NULL
            intervals <- NULL
            if (exists("GTRACKS", envir = .misha)) {
                tracks <- get("GTRACKS", envir = .misha)
            }
            if (exists("GINTERVS", envir = .misha)) {
                intervals <- get("GINTERVS", envir = .misha)
            }
            .gdefine_autocompletion_vars(tracks, intervals, gvtrack.ls(), interactive)
        } else {
            .gundefine_autocompletion_vars()
        }
    }

    options(.gautocompletion = autocompletion)
    options(.ginteractive = interactive)
}



#' Prints call stack of the last uncaught error
#'
#' Prints call stack of the last uncaught error in a friendly way.
#'
#' Similarly to 'traceback' this function prints the call stack of the last
#' uncaught error. Yet 'gtraceback' does it in a more friendly way by omitting
#' the calls that occurred inside the library.
#'
#' @param x see 'traceback'
#' @param max.lines see 'traceback'
#' @return See 'traceback'.
#' @seealso \code{\link{traceback}}
#' @keywords ~trace ~error ~exception
#' @examples
#'
#' gdb.init_examples()
#' f <- function() {
#'     gscreen("blablabla")
#' }
#' f()
#' gtraceback()
#'
#' @export gtraceback
gtraceback <- function(x = NULL, max.lines = getOption("deparse.max.lines")) {
    x <- NULL

    if (is.null(x) && (exists(".Traceback", envir = baseenv()))) {
        x <- get(".Traceback", envir = baseenv())
    }

    if (!is.null(x) && length(x) > 0) {
        # get the call stack and concatenate all complex commands together
        x <- sapply(x, paste, collapse = "")

        # extract call stack function names
        fnames <- gsub("^(\\S+)\\s*\\(.*\\)$", "\\1", x, perl = T)

        # get the indices of lib functions
        libindices <- which(fnames %in% get(".GFUNCS", envir = .misha))

        # cut whatever comes after the first lib function
        if (length(libindices) > 0) {
            x <- get(".Traceback", envir = .misha)[libindices[length(libindices)]:length(get(".Traceback", envir = .misha))]
        }
    }

    traceback(x, max.lines)
}




#' Downloads files from FTP server
#'
#' Downloads multiple files from FTP server
#'
#' This function downloads files from FTP server given by 'url'. The address in
#' 'url' can contain wildcards to download more than one file at once. Files
#' are downloaded to a directory given by 'path' argument.  If 'path' is
#' 'NULL', file are downloaded into 'GROOT/downloads'.
#'
#' @param url URL of FTP server
#' @param path directory path where the downloaded files are stored
#' @return An array of file names that have been downloaded.
#' @seealso \code{\link{gtrack.import_set}}
#' @keywords ~ftp
#' @examples
#'
#' gdb.init_examples()
#' gwget("ftp://hgdownload.cse.ucsc.edu/logs/2012/access_log.2012071*")
#'
#' @export gwget
gwget <- function(url = NULL, path = NULL) {
    if (is.null(url)) {
        stop("Usage: gwget(url, path = NULL)", call. = F)
    }

    if (is.null(path)) {
        .gcheckroot()
        path <- paste(get("GROOT", envir = .misha), "/downloads", sep = "")
        dir.create(path, showWarnings = F, recursive = T, mode = "0777")
    }

    if (!length(grep("^ftp\\:\\/\\/(\\w+(\\.\\w+)+)\\/(.+)", url, perl = T))) {
        url <- paste("ftp://", url, sep = "")
    }

    if (!length(grep("^ftp\\:\\/\\/(\\w+(\\.\\w+)+)\\/(.+)", url, perl = T))) {
        stop("Invalid format of URL", call. = F)
    }

    old.files <- dir(path)
    oldwd <- getwd()
    tryCatch(
        {
            setwd(path)

            server <- gsub("^ftp\\:\\/\\/(\\w+(\\.\\w+)+)\\/(.+)", "\\1", url, perl = T)
            rpath <- gsub("^ftp\\:\\/\\/(\\w+(\\.\\w+)+)\\/(.+)", "\\3", url, perl = T)

            rdir <- dirname(rpath)
            rfiles <- basename(rpath)
            if (system(paste(
                "/bin/sh -c \"",
                "ftp -i -n -v",
                server,
                "<< EOF\n",
                "user anonymous anonymous\n",
                "cd", rdir, "\n",
                paste("mget \"", rfiles, "\"", sep = ""), "\n",
                "quit\n",
                "EOF\""
            ))) {
                stop("Command failed", call. = F)
            }
        },
        interrupt = function(interrupt) {
            setwd(oldwd)
        },
        finally = {
            setwd(oldwd)
        }
    )

    new.files <- dir(path)
    setdiff(new.files, old.files)
}




#' Runs R commands on a cluster
#'
#' Runs R commands on a cluster that supports SGE.
#'
#' This function runs R commands on a cluster by distributing them among
#' cluster nodes. It must run on a machine that supports Sun Grid Engine (SGE).
#' The order in which the commands are executed can not be guaranteed,
#' therefore the commands must be inter-independent.
#'
#' Optional flags to 'qsub' command can be passed through 'opt.flags'
#' parameter. Users are strongly recommended to use only '-l' flag as other
#' flags might interfere with those that are already used (-terse, -S, -o, -e,
#' -V). For additional information please refer to the manual of 'qsub'.
#'
#' The maximal number of simultaneously submitted jobs is controlled by
#' 'max.jobs'.
#'
#' Set 'debug' argument to 'TRUE to allow additional report prints.
#'
#' 'gcluster.run' launches R on the cluster nodes to execute the commands. 'R'
#' argument specifies how R executable should be invoked.
#'
#' @param ... R commands
#' @param opt.flags optional flags for qsub command
#' @param max.jobs maximal number of simultaneously submitted jobs
#' @param debug if 'TRUE', additional reports are printed
#' @param R command that launches R
#' @return Return value ('retv') is a list, such that 'retv[[i]]' represents
#' the result of the run of command number 'i'. Each result consists of 4
#' fields that can be accessed by 'retv[[i]]$FIELDNAME':
#'
#' \tabular{ll}{ \emph{FIELDNAME} \tab \emph{DESCRIPTION}\cr exit.status \tab
#' Exit status of the command. Possible values: 'success', 'failure' or
#' 'interrupted'.\cr retv \tab Return value of the command.\cr stdout \tab
#' Standard output of the command.\cr stderr \tab Standard error of the
#' command. }
#' @keywords ~cluster
#' @examples
#' \dontrun{
#' gdb.init_examples()
#' v <- 17
#' gcluster.run(
#'     gsummary("dense_track + v"),
#'     {
#'         intervs <- gscreen("dense_track > 0.1", gintervals(1, 2))
#'         gsummary("sparse_track", intervs)
#'     },
#'     gsummary("rects_track")
#' )
#' }
#'
#' @export gcluster.run
gcluster.run <- function(..., opt.flags = "", max.jobs = 400, debug = FALSE, R = "R") {
    commands <- as.list(substitute(list(...))[-1L])

    if (length(commands) < 1) {
        stop("Usage: gculster.run(..., opt.flags = \"\" max.jobs = 400, debug = FALSE)", call. = F)
    }

    if (!length(system("which qsub", ignore.stderr = T, intern = T))) {
        stop("gcluster.run must run on a host that supports Sun Grid Engine (qsub)", call. = F)
    }

    .gcheckroot()

    tmp.dirname <- ""

    submitted.jobs <- c()

    tryCatch(
        {
            tmp.dirname <- tempfile(pattern = "", tmpdir = paste(get("GROOT", envir = .misha), "/tmp", sep = ""))
            if (!dir.create(tmp.dirname, recursive = T, mode = "0777")) {
                stop(sprintf("Failed to create a directory %s", tmp.dirname), call. = F)
            }

            # save the environment + options
            # parent.frame() is the environment of the caller
            cat("Preparing for distribution...\n")
            save(.misha$.GLIBDIR, file = paste(tmp.dirname, "libdir", sep = "/"))
            vars <- ls(all.names = TRUE, envir = parent.frame())
            envir <- parent.frame()
            while (!identical(envir, .GlobalEnv)) {
                envir <- parent.env(envir)
                vars <- union(vars, ls(all.names = TRUE, envir = envir))
            }
            save(list = vars, file = paste(tmp.dirname, "envir", sep = "/"), envir = parent.frame())
            .GSGECMD <- commands
            save(.GSGECMD, file = paste(tmp.dirname, "commands", sep = "/"))
            opts <- options()
            save(opts, file = paste(tmp.dirname, "opts", sep = "/"))

            cat("Running the commands...\n")
            completed.jobs <- c()
            progress <- -1
            repeat {
                # submit the commands
                num.running.jobs <- length(submitted.jobs) - length(completed.jobs)
                if (length(submitted.jobs) < length(commands) && num.running.jobs < max.jobs) {
                    istart <- length(submitted.jobs) + 1
                    iend <- min(length(commands), istart + (max.jobs - num.running.jobs) - 1)

                    for (i in istart:iend) {
                        out.file <- sprintf("%s/%d.out", tmp.dirname, i)
                        err.file <- sprintf("%s/%d.err", tmp.dirname, i)
                        script <- paste(get(".GLIBDIR", envir = .misha), "exec", "sgjob.sh", sep = "/")
                        command <- sprintf(
                            "unset module; qsub -terse -S /bin/bash -o %s -e %s -V %s %s %d '%s' '%s'",
                            out.file, err.file, opt.flags, script, i, tmp.dirname, R
                        )
                        jobid <- system(command, intern = TRUE)
                        if (length(jobid) != 1) {
                            stop("Failed to run qsub", call. = FALSE)
                        }
                        if (debug) {
                            cat(sprintf("\tSubmitted job %d (id: %s)\n", i, jobid))
                        }
                        submitted.jobs <- c(submitted.jobs, jobid)
                    }
                }

                # monitor progress
                Sys.sleep(3)
                running.jobs <- .gcluster.running.jobs(submitted.jobs)

                old.completed.jobs <- completed.jobs
                completed.jobs <- setdiff(submitted.jobs, running.jobs)
                if (debug) {
                    delta.jobs <- setdiff(completed.jobs, old.completed.jobs)
                    if (length(delta.jobs) > 0) {
                        for (jobid in delta.jobs) {
                            cat(sprintf("\tJob %d (id: %s) completed\n", match(jobid, submitted.jobs), jobid))
                        }
                    }

                    if (!length(running.jobs) && length(submitted.jobs) == length(commands)) {
                        break
                    }

                    new.progress <- length(completed.jobs)
                    if (new.progress != progress) {
                        progress <- new.progress
                        cat(sprintf("\t%d job(s) still in progress\n", length(commands) - progress))
                    }
                } else {
                    if (!length(running.jobs) && length(submitted.jobs) == length(commands)) {
                        break
                    }

                    new.progress <- as.integer(100 * length(completed.jobs) / length(commands))
                    if (new.progress != progress) {
                        progress <- new.progress
                        cat(sprintf("%d%%...", progress))
                    } else {
                        cat(".")
                    }
                }
            }
            if (!debug && progress != -1 && progress != 100) {
                cat("100%\n")
            }
        },
        interrupt = function(interrupt) {
            cat("\n")
            stop("Command interrupted!", call. = FALSE)
        },
        finally = {
            # We should perform now cleanup. If the user presses again Ctr+C "finaly" statement will be interrupted and the cleanup will
            # be incomplete. Unfortunately there is no way to block interrupts up until "finally" is done.
            # The only way is to solve the problem is to run "finally" in a loop.
            cleanup.finished <- FALSE
            while (!cleanup.finished) {
                tryCatch(
                    {
                        # kill still running jobs
                        if (length(submitted.jobs) > 0) {
                            # pack the answer
                            running.jobs <- .gcluster.running.jobs(submitted.jobs)
                            answer <- c()
                            for (i in 1:length(commands)) {
                                res <- list()
                                res$exit.status <- NA
                                res$retv <- NA
                                res$stdout <- NA
                                res$stderr <- NA

                                if (submitted.jobs[i] %in% running.jobs) {
                                    res$exit.status <- "interrupted"
                                } else {
                                    fname <- sprintf("%s/%d.retv", tmp.dirname, i)
                                    if (file.exists(fname)) {
                                        load(fname)
                                        res$exit.status <- "success"
                                        res$retv <- retv
                                    } else {
                                        res$exit.status <- "failure"
                                    }
                                }

                                out.file <- sprintf("%s/%d.out", tmp.dirname, i)
                                if (file.exists(out.file)) {
                                    f <- file(out.file, "rc")
                                    res$stdout <- readChar(f, 1e6)
                                    close(f)
                                }

                                err.file <- sprintf("%s/%d.err", tmp.dirname, i)
                                if (file.exists(err.file)) {
                                    f <- file(err.file, "rc")
                                    res$stderr <- readChar(f, 1e6)
                                    close(f)
                                }

                                answer[[i]] <- res
                            }
                            for (job in running.jobs) {
                                system(sprintf("qdel %s", job), ignore.stderr = T, intern = T)
                            }

                            unlink(tmp.dirname, recursive = TRUE)
                            return(answer)
                        }
                        unlink(tmp.dirname, recursive = TRUE)
                        cleanup.finished <- TRUE
                    },
                    interrupt = function(interrupt) {}
                )
            }
        }
    )
}


.gcluster.running.jobs <- function(jobids) {
    str <- system("qstat | sed 's/^[ ]*//' | cut -f 1 -d\" \"", intern = T)
    if (length(str) > 2) {
        intersect(jobids, str)
    } else {
        c()
    }
}
