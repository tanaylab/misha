.gcall <- function(...) {
    # The number of misha calls running in this process (this one not counted): above zero, this
    # call is nested in another one, made by R code that the outer call runs (a track expression,
    # the FUN of gintervals.mapply), in the parent or in one of its kids. Kept in R, and restored
    # by on.exit, so that an error raised inside a .Call cannot leave it stuck: such an error is a
    # longjmp that skips ~RdbInitializer, which leaves the C++ counters above zero (and the umask,
    # the signal handlers and the multitasking state in the shape that call left them). From
    # 5.16.0 that turned every later call of the session serial. A call made at depth zero has no
    # misha frame below it - every entry point that evaluates R code is reached through here - so
    # C_lifetime_repair() puts back what the leaked scope did not before the call starts.
    depth <- .gcall_depth()
    if (depth == 0L) {
        .Call("C_lifetime_repair")
    }
    on.exit(assign(".GCALL_DEPTH", depth, envir = .misha), add = TRUE)
    assign(".GCALL_DEPTH", depth + 1L, envir = .misha)

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
    outer_write_error <- .Call("C_take_write_error", NULL)
    completed <- FALSE
    on.exit(
        {
            mine <- .gtake_pending_diagnostics()
            # A file the call failed to finish writing (see BufferedFile::write_error): the
            # final flush happens when C++ closes the file, often in a destructor that cannot
            # raise anything. Taken on every exit, so a failed call does not leave it behind
            # for the next one; an outer call's (this one is nested in it) is put back.
            write_error <- .Call("C_take_write_error", outer_write_error)
            if (length(outer)) {
                assign(".GPENDING.DIAGNOSTICS", outer, envir = .misha)
            }
            # A call that did not return produced no results for a diagnostic to
            # qualify, so its queue is dropped rather than raised - both because the
            # error is the thing worth reporting and because a warning raised while an
            # error unwinds would replace it under options(warn = 2).
            if (completed) {
                if (!is.null(write_error)) {
                    stop(write_error, call. = FALSE)
                }
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

# Closes a connection that wrote `path`. close() writes the last buffer, and if that fails R
# only warns and returns a non-zero status - a truncated file behind a call that seems to have
# succeeded. Removes the file and raises an error naming it instead.
.gclose_written <- function(con, path) {
    why <- ""
    status <- withCallingHandlers(close(con), warning = function(w) {
        if (startsWith(conditionMessage(w), "Problem closing connection")) {
            why <<- trimws(sub("^Problem closing connection:", "", conditionMessage(w)))
            invokeRestart("muffleWarning")
        }
    })
    if (!is.null(status) && status != 0) {
        unlink(path)
        stop(sprintf("Failed to write %s: %s", path, why), call. = FALSE)
    }
    invisible()
}

# Serializes x into a new file at `path`. A write that fails midway, fails when close() writes
# the last buffer, or is interrupted removes the file and raises an error naming it, instead of
# leaving a truncated file (and an open connection) behind.
.gserialize_file <- function(x, path) {
    f <- file(path, "wb")
    tryCatch(serialize(x, f),
        error = function(e) {
            suppressWarnings(close(f))
            unlink(path)
            stop(sprintf("Failed to write %s: %s", path, conditionMessage(e)), call. = FALSE)
        },
        interrupt = function(i) {
            suppressWarnings(close(f))
            unlink(path)
            stop("Command interrupted!", call. = FALSE)
        }
    )
    .gclose_written(f, path)
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

# The number of .gcall()s running in this process (see .gcall)
.gcall_depth <- function() {
    get0(".GCALL_DEPTH", envir = .misha, inherits = FALSE, ifnotfound = 0L)
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
            val <- .misha_defaults[[x]]
        }
    } else if (x %in% names(options())) {
        val <- options(x)[[1L]]
    } else {
        val <- default
    }
    if (identical(x, "gmultitasking")) {
        # if (.ggetOption("gmultitasking")) would stop with "missing value where TRUE/FALSE needed"
        if (length(val) == 1 && is.na(val)) {
            stop("options(gmultitasking) must be TRUE or FALSE, not NA", call. = FALSE)
        }
        # A misha call made by R code that another misha call runs (the FUN of
        # gintervals.mapply, a track expression), in the parent or in one of its kids, runs
        # serially: forking would reuse the outer call's shared memory and kid slots. Told by
        # .gcall()'s depth, not by RdbInitializer's ref_count, which a scope that ended through a
        # longjmp leaves above zero until the next call repairs it.
        if (.gcall_depth() > 0L) {
            val <- FALSE
        }
    }
    val
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
    # A tibble, grouped tibble or data.table is bound as the plain data frame it is; their
    # own attributes (groups, the data.table self reference) would not describe the result.
    objs <- lapply(list(...), function(obj) if (is.data.frame(obj)) as.data.frame(obj) else obj)

    zerolines <- lapply(objs, function(obj) {
        obj[0, ]
    })

    diffs <- sapply(zerolines, FUN = attr.all.equal, zerolines[[1]])
    if (!all(sapply(diffs, FUN = is.null))) {
        stop("Cannot rbind objects: columns differ", call. = FALSE)
    }

    if (length(objs) > 1) {
        objs <- .grbind_unify_cols(objs)
    }
    .gcall("grbind", objs, .misha_env())
}

# The C++ grbind fills each column of the result with the values of every object as they are
# stored, in a column made like the first object's: a factor's level codes would be read with the
# first object's levels (rows silently got another chromosome, or another label), and a column of
# another type ended in an R error inside the .Call. So each column is given one type here first,
# the one base rbind() gives it:
#   - factors, and characters bound to a factor, become one factor whose levels are those of the
#     factors and the values of the characters, in the order met (NA is a level only when a factor
#     has it);
#   - factors bound to a character column become characters;
#   - logical, integer and double columns become the widest of them.
# Other columns whose types differ (a number and a string, a Date and a number) are an error. As in
# rbind(), objects with no rows take no part in the choice when another object has rows; they get
# an empty column of the chosen type.
.grbind_unify_cols <- function(objs) {
    nrows <- vapply(objs, function(obj) .row_names_info(obj, 2L), integer(1))
    donors <- if (any(nrows > 0L)) which(nrows > 0L) else seq_along(objs)
    for (j in seq_along(objs[[1]])) {
        cols <- lapply(objs, .subset2, j)
        if (.grbind_same_type(cols)) {
            next
        }
        unified <- .grbind_unify_col(cols[donors], donors, names(objs[[1]])[j])
        for (i in setdiff(seq_along(objs), donors)) {
            objs[[i]][[j]] <- unified[[1]][0]
        }
        for (k in seq_along(donors)) {
            objs[[donors[k]]][[j]] <- unified[[k]]
        }
    }
    objs
}

# TRUE when the columns have one type, one class and (factors) one set of levels
.grbind_same_type <- function(cols) {
    first <- cols[[1]]
    type <- typeof(first)
    cls <- oldClass(first)
    levs <- if (is.factor(first)) levels(first)
    for (x in cols[-1]) {
        if (typeof(x) != type || !identical(oldClass(x), cls) || (!is.null(levs) && !identical(levels(x), levs))) {
            return(FALSE)
        }
    }
    TRUE
}

# cols: the column of each object that has rows, idx: the positions of these objects
.grbind_unify_col <- function(cols, idx, name) {
    if (.grbind_same_type(cols)) {
        return(cols)
    }
    first <- cols[[1]]
    is_fac <- vapply(cols, is.factor, logical(1))
    is_chr <- vapply(cols, is.character, logical(1))

    if (is.factor(first) && all(is_fac | is_chr)) {
        levs <- unique(unlist(lapply(cols, function(x) if (is.factor(x)) levels(x) else x), use.names = FALSE))
        na_level <- any(vapply(cols[is_fac], function(x) anyNA(levels(x)), logical(1)))
        if (!na_level) {
            levs <- levs[!is.na(levs)]
        }
        cls <- if (all(vapply(cols[is_fac], is.ordered, logical(1)))) c("ordered", "factor") else "factor"
        return(lapply(cols, function(x) {
            if (is.factor(x)) {
                codes <- match(levels(x), levs)[unclass(x)]
                if (na_level) {
                    codes[is.na(unclass(x))] <- match(NA_character_, levs)
                }
            } else {
                codes <- match(x, levs)
            }
            structure(codes, levels = levs, class = cls)
        }))
    }

    if (is.character(first) && all(is_fac | is_chr)) {
        return(lapply(cols, function(x) if (is.factor(x)) as.character(x) else x))
    }

    num_types <- c("logical", "integer", "double")
    types <- vapply(cols, typeof, character(1))
    if (all(types %in% num_types) && all(vapply(cols, function(x) is.null(oldClass(x)), logical(1)))) {
        type <- num_types[max(match(types, num_types))]
        return(lapply(cols, function(x) {
            storage.mode(x) <- type
            x
        }))
    }

    describe <- function(x) if (is.factor(x)) "a factor" else if (is.object(x)) sprintf("%s (%s)", class(x)[1], typeof(x)) else class(x)[1]
    other <- which(!vapply(cols, function(x) identical(class(x), class(first)) && typeof(x) == typeof(first), logical(1)))[1]
    stop(sprintf(
        "Cannot rbind objects: column '%s' is %s in object %d and %s in object %d",
        name, describe(first), idx[1], describe(cols[[other]]), idx[other]
    ), call. = FALSE)
}

# options(gmax.data.size) as the C++ side reads it (ConfigurationManager::get_max_data_size): unset
# (NULL, documented as no limit) or not a number is Inf, NA is an error
.gmax_data_size <- function() {
    size <- .ggetOption("gmax.data.size")
    if (length(size) && is.na(size[1])) {
        stop("options(gmax.data.size) must be a number, not NA", call. = FALSE)
    }
    if (!is.numeric(size) || !length(size)) Inf else size[1]
}

# A number (a coordinate, a count, gmax.data.size) for a message, as written: sprintf("%d") refuses
# doubles above 2^31, and %g prints 61331566 as 6.13316e+07
.gformat_num <- function(x) {
    format(x, scientific = FALSE, trim = TRUE)
}

.gverify_max_data_size <- function(size, data_name = "Result", arguments = NULL) {
    max.data.size <- .gmax_data_size()

    if (size > max.data.size) {
        max.data.size.str <- .gformat_num(max.data.size)
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


# Called by exec/sgjob.sh to root a job where the caller was. `file` is the
# copy of .misha that gcluster.run saved; loading it only yields a detached
# environment that misha never reads, so it has to be replayed through the API.
.gcluster.restore_db <- function(file) {
    saved <- new.env()
    load(file, envir = saved)
    state <- saved$.misha
    # The example db that gdb.init_examples() extracted into the caller's tempdir,
    # which other hosts do not see: keep the job's own copy.
    if (!dir.exists(state$GROOT) && grepl("/trackdb/test$", state$GROOT)) {
        return(invisible())
    }
    gsetroot(state$GROOT)
    # not gsetroot(dir = GWD): it refuses the db root itself, which gdir.cd("..") allows
    if (!identical(state$GWD, .misha$GWD)) {
        .gdir.cd(state$GWD, FALSE)
    }
    for (dataset in state$GDATASETS) {
        gdataset.load(dataset, force = TRUE)
    }
    if (!is.null(state$GVTRACKS)) {
        assign("GVTRACKS", state$GVTRACKS, envir = .misha)
    }
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

    # The whole-genome scope passed as a whole (.misha$ALLGENOME) keeps its 2D part as it is, the
    # placeholder above gmulticontig.2d.threshold contigs. The C++ conversion builds it in
    # options(gmulticontig.2d.mode) when 2D intervals are asked for, or stops without the option; for a
    # computation of the track expression scanner (.gscope_lazy_2d), only when it is 2D.
    return(intervals)
}

# The message of a computation over the deferred whole-genome 2D scope (as DEFERRED_2D_SCOPE_MSG in C++)
.DEFERRED_2D_SCOPE_MSG <- paste(
    "The whole-genome 2D scope of a database with more than getOption(\"gmulticontig.2d.threshold\") contigs is not built",
    "by default, and a 2D computation over it would miss the trans chromosome pairs. Set",
    "options(gmulticontig.2d.mode = \"full\") for all chromosome pairs, options(gmulticontig.2d.mode = \"diagonal\")",
    "for the cis pairs only, or give the 2D intervals explicitly."
)

# The whole-genome scope (.misha$ALLGENOME), or its deferred 2D part alone, for an interval operation
# done without the track expression scanner (in R or in C++): the 2D part, deferred above
# gmulticontig.2d.threshold contigs, built in options(gmulticontig.2d.mode), or an error without it.
# Anything else is returned as it is.
.gscope_expand_2d <- function(intervals) {
    placeholder <- function(x) is.list(x) && !is.data.frame(x) && !is.null(attr(x, "deferred"))
    whole <- is.list(intervals) && !is.data.frame(intervals) && length(intervals) == 2 && placeholder(intervals[[2]])
    if (!whole && !placeholder(intervals)) {
        return(intervals)
    }
    intervals2d <- .gintervals.2d.whole_genome()
    if (is.null(intervals2d)) {
        stop(.DEFERRED_2D_SCOPE_MSG, call. = FALSE)
    }
    if (whole) {
        intervals[[2]] <- intervals2d
        intervals
    } else {
        intervals2d
    }
}

# The whole-genome scope (.misha$ALLGENOME) for a computation of the track expression scanner: its
# deferred 2D part (the placeholder above gmulticontig.2d.threshold contigs) is marked "lazy", so that
# the C++ conversion leaves it to the scanner, which builds it only for a 2D computation
# (TrackExprScanner::build_deferred_2d_scope). Unmarked, a conversion that asks for 2D intervals
# builds it at once or stops. The deferred 2D part alone (.misha$ALLGENOME[[2]]) is a 2D scope, so it is
# built at once (.gscope_expand_2d). Anything else is returned as it is.
.gscope_lazy_2d <- function(intervals) {
    if (is.list(intervals) && !is.data.frame(intervals) && !is.null(attr(intervals, "deferred"))) {
        return(.gscope_expand_2d(intervals))
    }
    if (is.list(intervals) && !is.data.frame(intervals) && length(intervals) == 2 &&
        is.list(intervals[[2]]) && !is.data.frame(intervals[[2]]) && !is.null(attr(intervals[[2]], "deferred"))) {
        attr(intervals[[2]], "lazy") <- TRUE
    }
    intervals
}
