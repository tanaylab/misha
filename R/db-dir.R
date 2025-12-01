# Directory management functions

#' Changes current working directory in Genomic Database
#'
#' Changes current working directory in Genomic Database.
#'
#' This function changes the current working directory in Genomic Database (not
#' to be confused with shell's current working directory). The list of database
#' objects - tracks, intervals, track variables - is rescanned recursively
#' under 'dir'. Object names are updated with the respect to the new current
#' working directory. Example: a track named 'subdir.dense' will be referred as
#' 'dense' once current working directory is set to 'subdir'. All virtual
#' tracks are removed.
#'
#' @param dir directory path
#' @return None.
#' @seealso \code{\link{gdb.init}}, \code{\link{gdir.cwd}},
#' \code{\link{gdir.create}}, \code{\link{gdir.rm}}
#' @keywords ~db ~data ~database ~cd ~dir ~directory ~folder
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gdir.cd("subdir")
#' gtrack.ls()
#' gdir.cd("..")
#' gtrack.ls()
#'
#' @export gdir.cd
gdir.cd <- function(dir = NULL) {
    if (is.null(dir)) {
        stop("Usage: gdir.cd(dir)", call. = FALSE)
    }

    success <- FALSE
    oldgwd <- get("GWD", envir = .misha)

    tryCatch(
        {
            .gdir.cd(dir, TRUE)
            success <- TRUE
        },
        finally = {
            if (!success) {
                .gdir.cd(oldgwd, TRUE)
            }
        }
    )
}


#' Creates a new directory in Genomic Database
#'
#' Creates a new directory in Genomic Database.
#'
#' This function creates a new directory in Genomic Database. Creates only the
#' last element in the specified path.
#'
#' @param dir directory path
#' @param showWarnings see 'dir.create'
#' @param mode see 'dir.create'
#' @return None.
#' @note A new directory cannot be created within an existing track directory.
#' @seealso \code{\link{dir.create}}, \code{\link{gdb.init}},
#' \code{\link{gdir.cwd}}, \code{\link{gdir.rm}}
#' @keywords ~db ~data ~database ~dir ~directory ~folder ~create
#' @export gdir.create
gdir.create <- function(dir = NULL, showWarnings = TRUE, mode = "0777") {
    if (is.null(dir)) {
        stop("Usage: gdir.create(dir, showWarnings = TRUE, mode = \"0777\")", call. = FALSE)
    }

    oldwd <- getwd()
    on.exit(setwd(oldwd), add = TRUE)
    setwd(get("GWD", envir = .misha))
    tryCatch(
        {
            d <- dirname(dir)

            if (!file.exists(d)) {
                stop(sprintf("Path %s does not exist.\nNote: recursive directory creation is forbidden.", d), call. = FALSE)
            }

            t <- .gfindtrackinpath(d)
            if (!is.null(t)) {
                stop(sprintf("Cannot create a directory within a track %s", t), call. = FALSE)
            }

            if (length(grep("\\.track$", basename(dir))) > 0) {
                stop("gdir.create cannot create track directories", call. = FALSE)
            }

            dir.create(dir, showWarnings = showWarnings, recursive = FALSE, mode = mode)
        },
        interrupt = function(interrupt) {
            setwd(oldwd)
        },
        finally = {
            setwd(oldwd)
        }
    )
}

#' Create directories needed for track creation
#'
#' @description This function creates the directories needed for track creation.
#' For example, if the track name is 'proj.sample.my_track', this function
#' creates the directories 'proj' and 'sample'. Use this function with caution -
#' a long track name may create a deep directory structure.
#'
#' @param track name of the track
#'
#' @inheritParams gdir.create
#'
#' @return None.
#' @examples
#'
#' gdb.init_examples()
#'
#' # This creates the directories 'proj' and 'sample'
#' gtrack.create_dirs("proj.sample.my_track")
#'
#' @export
gtrack.create_dirs <- function(track, mode = "0777") {
    # split the track name into directories
    dirs <- dirname(gsub("\\.", "/", track))
    dirs <- strsplit(dirs, "/")[[1]]
    dir <- dirs[1]
    for (i in 1:length(dirs)) {
        if (i > 1) {
            dir <- paste(dir, dirs[i], sep = "/")
        }
        gdir.create(dir, mode = mode)
    }
}


#' Returns the current working directory in Genomic Database
#'
#' Returns the absolute path of the current working directory in Genomic
#' Database.
#'
#' This function returns the absolute path of the current working directory in
#' Genomic Database (not to be confused with shell's current working
#' directory).
#'
#' @return A character string of the path.
#' @seealso \code{\link{gdb.init}}, \code{\link{gdir.cd}},
#' \code{\link{gdir.create}}, \code{\link{gdir.rm}}
#' @keywords ~db ~data ~database ~cwd ~pwd ~dir ~directory ~folder
#' @export gdir.cwd
gdir.cwd <- function() {
    .gcheckroot()
    get("GWD", envir = .misha)
}


#' Deletes a directory from Genomic Database
#'
#' Deletes a directory from Genomic Database.
#'
#' This function deletes a directory from Genomic Database. If 'recursive' is
#' 'TRUE', the directory is deleted with all the files/directories it contains.
#' If the directory contains tracks or intervals, the user is prompted to
#' confirm the deletion. Set 'force' to 'TRUE' to suppress the prompt.
#'
#' @param dir directory path
#' @param recursive if 'TRUE', the directory is deleted recursively
#' @param force if 'TRUE', suppresses user confirmation of tracks/intervals
#' removal
#' @return None.
#' @seealso \code{\link{gdb.init}}, \code{\link{gdir.create}},
#' \code{\link{gdir.cd}}, \code{\link{gdir.cwd}}
#' @keywords ~db ~data ~database ~dir ~directory ~folder ~rm
#' @export gdir.rm
gdir.rm <- function(dir = NULL, recursive = FALSE, force = FALSE) {
    if (is.null(dir)) {
        stop("Usage: gdir.rm(dir, recursive = FALSE, force = FALSE)", call. = FALSE)
    }

    oldwd <- getwd()
    on.exit(setwd(oldwd), add = TRUE)
    setwd(get("GWD", envir = .misha))
    tryCatch(
        {
            if (!file.exists(dir)) {
                if (force) {
                    return(invisible())
                }
                stop(sprintf("Directory %s does not exist", dir), call. = FALSE)
            }

            r <- file.info(dir)
            if (r[names(r) == "isdir"] != 1) {
                stop(sprintf("%s is not a directory", dir), call. = FALSE)
            }

            t <- .gfindtrackinpath(dir)
            if (!is.null(t)) {
                stop(sprintf("Directory %s belongs to track %s", dir, t), call. = FALSE)
            }

            answer <- "Y"

            if (recursive && !force) {
                res <- .gcall("gfind_tracks_n_intervals", dir, .misha_env())
                tracks <- res[[1]]
                intervals <- res[[2]]

                if (!force && length(tracks) + length(intervals) > 0) {
                    message(sprintf("Directory %s contains tracks or intervals. Are you still sure you want to delete it (Y/N)? ", dir))
                    answer <- toupper(readLines(n = 1))
                }
            }

            if (answer == "Y" || answer == "YES") {
                if (recursive) {
                    unlink(dir, recursive)
                } else {
                    file.remove(dir)
                }

                if (file.exists(dir)) {
                    stop("Failed to remove the directory", call. = FALSE)
                }
            }
            gdb.reload()
        },
        interrupt = function(interrupt) {
            setwd(oldwd)
        },
        finally = {
            setwd(oldwd)
        }
    )
}


#' Sets read-only track attributes
#'
#' Sets read-only track attributes.
#'
#' This function sets the list of read-only track attributes. The specified
#' attributes may or may not already exist in the tracks.
#'
#' If 'attrs' is 'NULL' the list of read-only attributes is emptied.
#'
#' @param attrs a vector of read-only attributes names or 'NULL'
#' @return None.
#' @seealso \code{\link{gdb.get_readonly_attrs}},
#' \code{\link{gtrack.attr.get}}, \code{\link{gtrack.attr.set}}
#' @keywords ~attr ~attribute
#' @export gdb.set_readonly_attrs
