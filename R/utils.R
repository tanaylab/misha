.gcall <- function(...) {
    # C++ cannot raise a warning or a message itself: a caller that catches either with
    # an exiting handler (tryCatch, or options(warn = 2)) would longjmp out of the C
    # frame and leave misha's PROTECT counter inflated. Diagnostics are appended to
    # .misha$.GPENDING.DIAGNOSTICS - a list of c(severity, text) pairs, see
    # rdb::add_pending_diagnostic - and raised here, in the order they were produced,
    # once the call has returned.
    #
    # An outer call may have queued diagnostics of its own: gintervals.mapply runs the
    # user's FUN, and any misha call inside it nests under it. Take the outer queue
    # aside and put it back, so this call raises only what it produced. That has to
    # happen on the way out of an error or an interrupt too, or the outer call loses
    # its own diagnostics and inherits whatever the failed inner call had queued.
    outer <- .gtake_pending_diagnostics()
    completed <- FALSE
    on.exit(
        {
            mine <- .gtake_pending_diagnostics()
            if (length(outer)) {
                assign(".GPENDING.DIAGNOSTICS", outer, envir = .misha)
            }
            # A call that did not return produced no results for a diagnostic to
            # qualify, so its queue is dropped rather than raised - both because the
            # error is the thing worth reporting and because a warning raised while an
            # error unwinds would replace it under options(warn = 2).
            if (completed) {
                .graise_diagnostics(mine)
            }
        },
        add = TRUE
    )
    # A multitasking child process must end inside C (rexit(), i.e. a signal), and must
    # never come back here: R would carry on running the caller's code inside the fork and
    # then leave through R's own shutdown, whose R_CleanTempDir() deletes the session
    # tempdir - taking with it any database that lives there (gdb.init_examples(), the
    # isolated test databases). A pid that changed across the .Call means this call forked
    # and we are the child, which is always a bug in the C layer. Die at once and hard:
    # SIGKILL runs no exit handler and removes no file, and the parent reports the child's
    # abnormal death instead of silently dropping its share of the result.
    #
    # pskill/SIGKILL are called unqualified, not as tools::pskill/tools::SIGKILL: a `::`
    # reference loads the namespace on first use if it is not loaded yet, which is file I/O
    # and R evaluation in a process already known to be in a broken state. @importFrom in
    # misha-package.R makes R load `tools` and resolve both bindings while the *parent*
    # loads misha - before any child exists - so the fork inherits them already resolved and
    # this call is nothing but a lookup in misha's own imports environment.
    pid <- Sys.getpid()
    tryCatch(
        {
            res <- .Call(...)
            if (!identical(Sys.getpid(), pid)) {
                pskill(Sys.getpid(), SIGKILL)
            }
        },
        interrupt = function(interrupt) {
            stop("Command interrupted!", call. = FALSE)
        }
    )
    completed <- TRUE
    res
}

# Removes the diagnostics queued by the C++ layer, if any, and returns them: a list of
# c(severity, text) character pairs, oldest first.
.gtake_pending_diagnostics <- function() {
    if (!exists(".GPENDING.DIAGNOSTICS", envir = .misha, inherits = FALSE)) {
        return(list())
    }
    diags <- get(".GPENDING.DIAGNOSTICS", envir = .misha)
    rm(list = ".GPENDING.DIAGNOSTICS", envir = .misha)
    diags
}

.graise_diagnostics <- function(diags) {
    for (diag in diags) {
        if (identical(diag[[1]], "message")) {
            message(diag[[2]])
        } else {
            warning(diag[[2]], call. = FALSE)
        }
    }
}

# RdbInitializer's two lifetime counters, as c(ref_count = , protect_count = ): the number of
# entry-point scopes currently alive and the depth of misha's own protection stack. Both are
# zero between calls, and both are restored by ~RdbInitializer on every exit - error and
# interrupt included. Either one stuck above zero means a .Call left through a longjmp that
# skipped that destructor, which is the failure mode behind the worst defects the package has
# shipped. Internal and unexported: it exists so tests can assert on the invariant itself
# rather than on a side effect of it.
.glifetime_counters <- function() {
    .Call("C_lifetime_counters")
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
            options(.ginteractive = FALSE)
            on.exit(options(.ginteractive = .ginteractive))
            do.call(FUN, list(...))
        },
        finally = {
            options(.ginteractive = .ginteractive)
        }
    )
}


.ggetOption <- function(x, default = NULL) {
    if (missing(default)) {
        val <- options(x)[[1L]]
        # If option not set, use centralized default from .misha_defaults
        if (is.null(val) && exists(".misha_defaults") && x %in% names(.misha_defaults)) {
            return(.misha_defaults[[x]])
        }
        return(val)
    }
    if (x %in% names(options())) {
        options(x)[[1L]]
    } else {
        default
    }
}

# Runs `expr` with misha's file-creation umask in force and restores the
# process umask afterwards.
#
# Why misha needs one at all: a genomic database is shared between the members
# of a group, so every file and directory misha creates in it has to stay
# group-writable. On a host with the common umask 022 default a database
# created without this would come out 644/755 and no collaborator could add a
# track to it.
#
# The default, "0007", produces 660 files and 770 directories - group rwx, no
# world access - which is what misha's C++ layer already sets around every
# write (RdbInitializer, src/rdbutils.cpp) and what the real shared databases
# look like on disk. Set options(gpermissions.umask = NULL) to leave the
# process umask alone, or e.g. "0002" for the world-readable 664/775 that
# misha produced before 5.11.21.
#
# The umask is process-global state, so it is held for the shortest possible
# span and restored by on.exit(), which runs on a normal return, on an error
# and on an interrupt alike. Never set it anywhere that outlives a call - it
# used to be set once in .onLoad() and never restored, which silently changed
# the permissions of everything the R session wrote, misha or not.
#
# getOption() is deliberate here: .ggetOption() falls back to .misha_defaults
# when the option reads back NULL, and R *removes* an option that is set to
# NULL, so .ggetOption() cannot tell "unset" from "explicitly turned off".
# .onLoad() seeds the live option from .misha_defaults instead.
.gwith_umask <- function(expr) {
    um <- getOption("gpermissions.umask")
    if (!is.null(um)) {
        old <- tryCatch(Sys.umask(um), error = function(e) {
            stop(sprintf(
                "Invalid gpermissions.umask option (%s): %s",
                paste(format(um), collapse = " "), conditionMessage(e)
            ), call. = FALSE)
        })
        on.exit(Sys.umask(old), add = TRUE)
    }
    expr
}

.gexpr2str <- function(x) {
    if (.ggetOption(".ginteractive", FALSE)) {
        if (is.null(substitute(x))) {
            NULL
        } else {
            str <- deparse(substitute(x), width.cutoff = 500)[1]
            gsub("^\"(.*)\"$", "\\1", str, perl = TRUE)
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
        stop("Cannot rbind objects: columns differ", call. = FALSE)
    }

    .gcall("grbind", objs, .misha_env())
}

.gverify_max_data_size <- function(size, data_name = "Result", arguments = NULL) {
    max.data.size <- .ggetOption("gmax.data.size")

    if (size > max.data.size) {
        # gmax.data.size is a (multi-GB) double; format as a plain integer string
        # so sprintf doesn't error with "%d" on a non-integer/over-2^31 value.
        max.data.size.str <- format(max.data.size, scientific = FALSE, trim = TRUE)
        if (is.null(arguments)) {
            stop(sprintf(
                paste("%s size exceeded the maximal allowed (%s).",
                    "Note: the maximum data size is controlled via gmax.data.size option (see options, getOptions).",
                    sep = "\n"
                ),
                data_name, max.data.size.str
            ), call. = FALSE)
        } else {
            stop(sprintf(
                paste("%s size exceeded the maximal allowed (%s).",
                    "Consider saving the result in a file (use %s argument).",
                    "Note: the maximum data size is controlled via gmax.data.size option (see options, getOptions).",
                    sep = "\n"
                ),
                data_name, max.data.size.str, paste(arguments, collapse = " or ")
            ), call. = FALSE)
        }
    }
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
#' gdb.init_examples()
#' \donttest{
#' outdir <- tempdir()
#' gwget("ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/md5sum.txt", path = outdir)
#' }
#'
#' @export gwget
gwget <- function(url = NULL, path = NULL) {
    if (is.null(url)) {
        stop("Usage: gwget(url, path = NULL)", call. = FALSE)
    }

    if (is.null(path)) {
        .gcheckroot()
        path <- paste(get("GROOT", envir = .misha), "/downloads", sep = "")
        .gwith_umask(dir.create(path, showWarnings = FALSE, recursive = TRUE, mode = "0777"))
    }

    if (!length(grep("^ftp\\:\\/\\/(\\w+(\\.\\w+)+)\\/(.+)", url, perl = TRUE))) {
        url <- paste("ftp://", url, sep = "")
    }

    if (!length(grep("^ftp\\:\\/\\/(\\w+(\\.\\w+)+)\\/(.+)", url, perl = TRUE))) {
        stop("Invalid format of URL", call. = FALSE)
    }

    gftp_download_glob(url, path, verbose = FALSE, handle = NULL)
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
#' @param control_dir directory where the control files are stored. Note that this
#' directory should be accessible from all nodes. If 'NULL', a temporary directory
#' would be created under the current misha database.
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
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#' \donttest{
#' gdb.init_examples()
#' # Run only on systems with Sun Grid Engine (SGE)
#' if (FALSE) {
#'     v <- 17
#'     gcluster.run(
#'         gsummary("dense_track + v"),
#'         {
#'             intervs <- gscreen("dense_track > 0.1", gintervals(1, 2))
#'             gsummary("sparse_track", intervs)
#'         },
#'         gsummary("rects_track")
#'     )
#' }
#' }
#'
#' @export gcluster.run
gcluster.run <- function(..., opt.flags = "", max.jobs = 400, debug = FALSE, R = "R", control_dir = NULL) {
    commands <- as.list(substitute(list(...))[-1L])

    if (length(commands) < 1) {
        stop("Usage: gcluster.run(..., opt.flags = \"\" max.jobs = 400, debug = FALSE)", call. = FALSE)
    }

    if (!length(system("which qsub", ignore.stderr = TRUE, intern = TRUE))) {
        stop("gcluster.run must run on a host that supports Sun Grid Engine (qsub)", call. = FALSE)
    }

    .gcheckroot()

    tmp.dirname <- ""
    if (is.null(control_dir)) {
        control_dir <- paste(get("GROOT", envir = .misha), "/tmp", sep = "")
    }

    # if the path of control dir starts with '/tmp' warn that the tempdir needs to be shared
    if (grepl("^/tmp", control_dir)) {
        message("Warning: The control directory is in /tmp. Make sure that it is shared between all nodes.")
    }

    submitted.jobs <- c()

    tryCatch(
        {
            tmp.dirname <- tempfile(pattern = "", tmpdir = control_dir)
            if (!.gwith_umask(dir.create(tmp.dirname, recursive = TRUE, mode = "0777"))) {
                stop(sprintf("Failed to create a directory %s", tmp.dirname), call. = FALSE)
            }

            # save the environment + options
            # parent.frame() is the environment of the caller
            message("Preparing for distribution...\n")

            save(.misha, file = paste(tmp.dirname, "misha", sep = "/"))
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

            message("Running the commands...")
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
                            message(sprintf("\tSubmitted job %d (id: %s)", i, jobid))
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
                            message(sprintf("\tJob %d (id: %s) completed", match(jobid, submitted.jobs), jobid))
                        }
                    }

                    if (!length(running.jobs) && length(submitted.jobs) == length(commands)) {
                        break
                    }

                    new.progress <- length(completed.jobs)
                    if (new.progress != progress) {
                        progress <- new.progress
                        message(sprintf("\t%d job(s) still in progress", length(commands) - progress))
                    }
                } else {
                    if (!length(running.jobs) && length(submitted.jobs) == length(commands)) {
                        break
                    }

                    new.progress <- as.integer(100 * length(completed.jobs) / length(commands))
                    if (new.progress != progress) {
                        progress <- new.progress
                        message(sprintf("%d%%...", progress), appendLF = FALSE)
                    } else {
                        message(".", appendLF = FALSE)
                    }
                }
            }
            if (!debug && progress != -1 && progress != 100) {
                message("100%")
            }
        },
        interrupt = function(interrupt) {
            message("\n")
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
                            for (i in seq_along(commands)) {
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
                                system(sprintf("qdel %s", job), ignore.stderr = TRUE, intern = TRUE)
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
    str <- system("qstat | sed 's/^[ ]*//' | cut -f 1 -d\" \"", intern = TRUE)
    if (length(str) > 2) {
        intersect(jobids, str)
    } else {
        c()
    }
}

rescue_ALLGENOME <- function(intervals, intervals_name) {
    if (length(intervals_name) == 0) {
        return(intervals)
    }
    if (intervals_name[1] == "ALLGENOME") {
        warning("ALLGENOME was deprecated at version 4.2.0. Use .misha$ALLGENOME instead.", call. = FALSE)
        intervals <- .misha$ALLGENOME
    }
    if (intervals_name[1] == "ALLGENOME[[1]]") {
        warning("ALLGENOME[[1]] was deprecated at version 4.2.0. Use .misha$ALLGENOME[[1]] instead.", call. = FALSE)
        intervals <- .misha$ALLGENOME[[1]]
    }
    if (intervals_name[1] == "ALLGENOME[[2]]") {
        warning("ALLGENOME[[2]] was deprecated at version 4.2.0. Use .misha$ALLGENOME[[2]] instead.", call. = FALSE)
        intervals <- .misha$ALLGENOME[[2]]

        # Check if 2D is deferred and generate on demand
        if (.is_2d_deferred(intervals)) {
            warning("2D genome structure is deferred for this multi-contig genome. Use gintervals.2d.all() to generate the full structure, or gintervals.2d() for specific pairs.", call. = FALSE)
            # Return deferred structure as-is - users should use gintervals.2d.all()
        }
    }

    # Handle deferred 2D structures when ALLGENOME is passed as a whole.
    if (is.list(intervals) && length(intervals) >= 2 && .is_2d_deferred(intervals[[2]])) {
        # When users pass ALLGENOME directly (the default in many APIs), the 2D
        # component might be a deferred placeholder. Expand it here so C++
        # interval conversions always see a proper data frame.
        intervals[[2]] <- if (identical(intervals, get("ALLGENOME", envir = .misha))) {
            gintervals.2d.all()
        } else {
            mode <- getOption("gmulticontig.2d.mode", "diagonal")
            .generate_2d_on_demand(intervals[[1]], mode)
        }
    }
    return(intervals)
}
