.gchroms <- function(chroms) {
    if (!is.character(chroms)) {
        chroms <- as.character(chroms)
    }

    idx <- substr(chroms, 1, 3) != "chr"
    chroms[idx] <- paste("chr", chroms[idx], sep = "")

    indices <- match(chroms, get("ALLGENOME")[[1]]$chrom)
    err.chroms <- chroms[is.na(indices)]
    if (length(err.chroms) > 0) {
        stop(sprintf("Chromosome %s does not exist in the database", err.chroms[1]))
    }
    get("ALLGENOME")[[1]]$chrom[indices] # return factor
}


.gcheckroot <- function() {
    if (!exists("GROOT", envir = .GlobalEnv) || !exists("ALLGENOME", envir = .GlobalEnv) || is.null(get("GROOT")) || is.null(get("ALLGENOME"))) {
        stop("Database root directory is not set. Please call gdb.init().", call. = F)
    }
}


.gdir.cd <- function(dir, rescan) {
    oldwd <- getwd()
    setwd(get("GWD"))
    tryCatch(
        {
            t <- .gfindtrackinpath(dir)
            if (!is.null(t)) {
                stop(sprintf("Directory %s belongs to track %s", dir, t), call. = F)
            }

            setwd(dir)
            newwd <- getwd()

            if (.ggetOption(".gautocompletion", FALSE)) {
                .gundefine_autocompletion_vars()
            }
            assign("GWD", newwd, envir = .GlobalEnv)
            setwd(oldwd)
            gdb.reload(rescan)
        },
        interrupt = function(interrupt) {
            setwd(oldwd)
        },
        finally = {
            setwd(oldwd)
        }
    )
}

#' @rdname gdb.init
#' @export
gsetroot <- function(groot = NULL, dir = NULL, rescan = FALSE) {
    if (is.null(groot)) {
        stop("Usage: gsetroot(groot, dir = NULL, rescan = FALSE)", call. = F)
    }

    groot <- normalizePath(groot)

    if (exists("GROOT") && exists("ALLGENOME") && !is.null(get("GROOT")) && !is.null(get("ALLGENOME")) && .ggetOption(".gautocompletion", FALSE)) {
        .gundefine_autocompletion_vars()
    }

    assign("ALLGENOME", NULL, envir = .GlobalEnv)
    assign("GROOT", NULL, envir = .GlobalEnv)

    chromsizes <- read.csv(paste(groot, "chrom_sizes.txt", sep = "/"), sep = "\t", header = F)
    colnames(chromsizes) <- c("chrom", "size")
    intervals <- data.frame(
        chrom = as.factor(paste("chr", as.character(chromsizes$chrom), sep = "")),
        start = 0, end = as.numeric(chromsizes$size)
    )

    if (nrow(intervals) == 0) {
        stop("chrom_sizes.txt file does not contain any chromosomes", call. = F)
    }

    for (chrom in intervals$chrom) {
        if (length(grep(sprintf("^%s$", chrom), intervals$chrom)) > 1) {
            stop(sprintf("Chromosome \"%s\" appears more than once in chrom_sizes.txt", chrom))
        }
    }
    intervals <- intervals[order(intervals$chrom), ]
    rownames(intervals) <- 1:nrow(intervals)

    cartesian <- expand.grid(1:nrow(intervals), 1:nrow(intervals))
    intervals2d <- cbind(intervals[cartesian[, 2], ], intervals[cartesian[, 1], ])
    names(intervals2d) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2")
    rownames(intervals2d) <- 1:nrow(intervals2d)

    assign("ALLGENOME", list(intervals, intervals2d), envir = .GlobalEnv)
    assign("GROOT", groot, envir = .GlobalEnv)
    assign("GWD", groot, envir = .GlobalEnv)

    success <- F
    tryCatch(
        {
            if (is.null(dir)) {
                .gdir.cd(paste(groot, "tracks", sep = "/"), rescan)
            } else {
                if (nchar(dir) < 1) {
                    stop("dir argument is an empty string")
                }

                c <- substr(dir, 1, 1)
                if (c == "~" || c == "/") {
                    .gdir.cd(dir, rescan)
                } else {
                    .gdir.cd(paste(groot, dir, sep = "/"), rescan)
                }
            }
            success <- T
        },
        finally = {
            if (!success) {
                assign("ALLGENOME", NULL, envir = .GlobalEnv)
                assign("GROOT", NULL, envir = .GlobalEnv)
                assign("GWD", NULL, envir = .GlobalEnv)
            }
        }
    )
}




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
#' \code{\link{gdir.create}}, \code{\link{gdir.rm}},
#' \code{\link{gset_input_mode}}
#' @keywords ~db ~data ~database ~cd ~dir ~directory ~folder
#' @examples
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
        stop("Usage: gdir.cd(dir)", call. = F)
    }

    success <- FALSE
    oldgwd <- get("GWD")

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
        stop("Usage: gdir.create(dir, showWarnings = TRUE, mode = \"0777\")", call. = F)
    }

    oldwd <- getwd()
    setwd(get("GWD"))
    tryCatch(
        {
            d <- dirname(dir)

            if (!file.exists(d)) {
                stop(sprintf("Path %s does not exist.\nNote: recursive directory creation is forbidden.", d), call. = F)
            }

            t <- .gfindtrackinpath(d)
            if (!is.null(t)) {
                stop(sprintf("Cannot create a directory within a track %s", t), call. = F)
            }

            if (length(grep("\\.track$", basename(dir))) > 0) {
                stop("gdir.create cannot create track directories", call. = F)
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
    get("GWD")
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
#' @param force if 'TRUE', supresses user confirmation of tracks/intervals
#' removal
#' @return None.
#' @seealso \code{\link{gdb.init}}, \code{\link{gdir.create}},
#' \code{\link{gdir.cd}}, \code{\link{gdir.cwd}}
#' @keywords ~db ~data ~database ~dir ~directory ~folder ~rm
#' @export gdir.rm
gdir.rm <- function(dir = NULL, recursive = FALSE, force = FALSE) {
    if (is.null(dir)) {
        stop("Usage: gdir.rm(dir, recursive = FALSE, force = FALSE)", call. = F)
    }

    oldwd <- getwd()
    setwd(get("GWD"))
    tryCatch(
        {
            if (!file.exists(dir)) {
                if (force) {
                    return(invisible())
                }
                stop(sprintf("Directory %s does not exist", dir), call. = F)
            }

            r <- file.info(dir)
            if (r[names(r) == "isdir"] != 1) {
                stop(sprintf("%s is not a directory", dir), call. = F)
            }

            t <- .gfindtrackinpath(dir)
            if (!is.null(t)) {
                stop(sprintf("Directory %s belongs to track %s", dir, t), call. = F)
            }

            answer <- "Y"

            if (recursive && !force) {
                res <- .gcall("gfind_tracks_n_intervals", dir, new.env(parent = parent.frame()), silent = TRUE)
                tracks <- res[[1]]
                intervals <- res[[2]]

                if (!force && length(tracks) + length(intervals) > 0) {
                    cat(sprintf("Directory %s contains tracks or intervals. Are you still sure you want to delete it (Y/N)? ", dir))
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
                    stop("Failed to remove the directory", call. = F)
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
gdb.set_readonly_attrs <- function(attrs) {
    .gcheckroot()

    filename <- paste(get("GROOT"), ".ro_attributes", sep = "/")

    if (is.null(attrs)) {
        unlink(filename)
    } else {
        attrs <- as.character(attrs)

        idx <- which(duplicated(attrs))[1]
        if (!is.na(idx)) {
            stop(sprintf("Attribute %s appears more than once", attrs[idx]), call. = F)
        }

        idx <- which(attrs == "")[1]
        if (!is.na(idx)) {
            stop("Attribute name cannot be an empty string", call. = F)
        }

        f <- file(filename, "wb")
        serialize(attrs, f)
        close(f)
    }
    retv <- 0 # suppress return value
}




#' Creates a new Genomic Database
#'
#' Creates a new Genomic Database.
#'
#' This function creates a new Genomic Database at the location specified by
#' 'groot'. FASTA files are converted to 'Seq' format and appropriate
#' 'chrom_sizes.txt' file is generated (see "User Manual" for more details).
#'
#' If 'genes.file' is not 'NULL' four sets of intervals are created in the
#' database: \code{tss}, \code{exons}, \code{utr3} and \code{utr5}. See
#' \link{gintervals.import_genes} for more details about importing genes
#' intervals.
#'
#' 'fasta', 'genes.file' and 'annots.file' can be either a file path or URL in
#' a form of 'ftp://[address]/[file]'. 'fasta' can also contain wildcards to
#' indicate multiple files. Files that these arguments point to can be zipped
#' or unzipped.
#'
#' @param groot path to newly created database
#' @param fasta an array of names or URLs of FASTA files. Can contain wildcards
#' for multiple files
#' @param genes.file name or URL of file that contains genes. If 'NULL' no
#' genes are imported
#' @param annots.file name of URL file that contains annotations. If 'NULL' no
#' annotations are imported
#' @param annots.names annotations names
#' @return None.
#' @seealso \code{\link{gdb.init}}, \code{\link{gdb.reload}},
#' \code{\link{gintervals.import_genes}}
#' @keywords ~database ~create ~genes
#' @examples
#' \dontrun{
#' ftp <- "ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19"
#' gdb.create(
#'     "mydb",
#'     c(
#'         paste(ftp, "chromosomes/chr2*", sep = "/"),
#'         paste(ftp, "chromosomes/chr7.*", sep = "/")
#'     ),
#'     paste(ftp, "database/knownGene.txt.gz", sep = "/"),
#'     paste(ftp, "database/kgXref.txt.gz", sep = "/"),
#'     c(
#'         "kgID", "mRNA", "spID", "spDisplayID", "geneSymbol",
#'         "refseq", "protAcc", "description", "rfamAcc",
#'         "tRnaName"
#'     )
#' )
#' gdb.init("mydb")
#' gintervals.ls()
#' gintervals.all()
#' }
#'
#' @export gdb.create
gdb.create <- function(groot = NULL, fasta = NULL, genes.file = NULL, annots.file = NULL, annots.names = NULL) {
    if (is.null(groot) || is.null(fasta)) {
        stop("Usage: gdb.create(groot, fasta, genes.file = NULL, annots.file = NULL, annots.names = NULL)", call. = F)
    }

    if (file.exists(groot)) {
        stop(sprintf("Directory %s already exists", groot), call. = F)
    }

    success <- FALSE
    allgenome.old <- NULL
    groot.old <- NULL
    if (exists("ALLGENOME")) {
        allgenome.old <- get("ALLGENOME")
    }
    if (exists("GROOT")) {
        groot.old <- get("GROOT")
    }

    tryCatch(
        {
            dir.create(groot, showWarnings = F, recursive = TRUE, mode = "0777")
            dir.create(paste(groot, "pssms", sep = "/"), showWarnings = F, recursive = TRUE, mode = "0777")
            dir.create(paste(groot, "seq", sep = "/"), showWarnings = F, recursive = TRUE, mode = "0777")
            dir.create(paste(groot, "tracks", sep = "/"), showWarnings = F, recursive = TRUE, mode = "0777")

            chroms <- .gseq.import(groot, fasta)

            if (!length(chroms)) {
                stop("No FASTA files were imported", call. = F)
            }

            seq.files <- paste("chr", chroms, ".seq", sep = "")
            seq.files <- paste(paste(groot, "seq", sep = "/"), seq.files, sep = "/")
            chrom.sizes <- data.frame(chrom = chroms, size = file.info(seq.files)$size)
            write.table(chrom.sizes, paste(groot, "chrom_sizes.txt", sep = "/"), quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)

            # before calling gintervals.import_genes new ALLGENOME must be set
            intervals <- data.frame(
                chrom = as.factor(paste("chr", as.character(chrom.sizes$chrom), sep = "")),
                start = 0, end = as.numeric(chrom.sizes$size)
            )
            intervals <- intervals[order(intervals$chrom), ]
            rownames(intervals) <- 1:nrow(intervals)

            cartesian <- expand.grid(1:nrow(intervals), 1:nrow(intervals))
            intervals2d <- cbind(intervals[cartesian[, 2], ], intervals[cartesian[, 1], ])
            names(intervals2d) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2")
            rownames(intervals2d) <- 1:nrow(intervals2d)

            assign("ALLGENOME", list(intervals, intervals2d), envir = .GlobalEnv)
            assign("GROOT", groot, envir = .GlobalEnv)

            if (!is.null(genes.file)) {
                intervs <- gintervals.import_genes(genes.file, annots.file, annots.names)
                if (!is.null(intervs)) {
                    for (i in 1:length(intervs)) {
                        if (!is.null(intervs$tss)) {
                            .gcall_noninteractive(.gintervals.save_file, sprintf("%s/tracks/%s.interv", groot, names(intervs)[i]), intervs[[i]])
                        }
                    }
                }
            }

            # write read-only attributes
            f <- file(paste(groot, ".ro_attributes", sep = "/"), "wb")
            serialize(c("created.by", "created.date"), f)
            close(f)

            cat("Database was successfully created\n")
            success <- TRUE
        },
        finally = {
            assign("ALLGENOME", allgenome.old, envir = .GlobalEnv)
            assign("GROOT", groot.old, envir = .GlobalEnv)
            if (!success) {
                unlink(groot, recursive = TRUE)
            }
        }
    )
    retv <- 0 # suppress return value
}



#' Returns a list of read-only track attributes
#'
#' Returns a list of read-only track attributes.
#'
#' This function returns a list of read-only track attributes. These attributes
#' are not allowed to be modified or deleted.
#'
#' If no attributes are marked as read-only a 'NULL' is returned.
#'
#' @return A list of read-only track attributes.
#' @seealso \code{\link{gdb.set_readonly_attrs}},
#' \code{\link{gtrack.attr.get}}, \code{\link{gtrack.attr.set}}
#' @keywords ~attr ~attribute
#' @export gdb.get_readonly_attrs
gdb.get_readonly_attrs <- function() {
    .gcheckroot()

    filename <- paste(get("GROOT"), ".ro_attributes", sep = "/")
    attrs <- NULL
    if (file.exists(filename)) {
        f <- file(filename, "rb")
        attrs <- unserialize(f)
        close(f)
        if (!is.character(attrs)) {
            stop(sprintf("Invalid format of read-only atrributes file %s", filename), call. = F)
        }

        attrs <- unique(attrs)
        attrs <- attrs[attrs != ""]
    }
    attrs
}



#' Initializes connection with Genomic Database
#'
#' Initializes connection with Genomic Database: loads the list of tracks,
#' intervals, etc.
#'
#' 'gdb.init' initializes the connection with the Genomic Database. It is
#' typically called first prior to any other function. When the package is
#' attached it internally calls to 'gdb.init.examples' which opens the
#' connection with the database located at 'PKGDIR/trackdb/test' directory,
#' where 'PKGDIR' is the directory where the package is installed.
#'
#' The current working directory inside the Genomic Database is set to 'dir'.
#' If 'dir' is 'NULL', the current working directory is set to 'GROOT/tracks'.
#'
#' If 'rescan' is 'TRUE', the list of tracks and intervals is achieved by
#' rescanning directory structure under the current current working directory.
#' Otherwise 'gdb.init' attempts to use the cached list that resides in
#' 'groot/.db.cache' file.
#'
#' Upon completion the connection is established with the database. If
#' auto-completion mode is switched on (see 'gset_input_method') the list of
#' tracks and intervals sets is loaded and added as variables to the global
#' environment allowing auto-completion of object names with <TAB> key. Also a
#' few global variables are defined. These variables should not be modified by
#' user.
#'
#' \tabular{ll}{ GROOT \tab Root directory of Genomic Database\cr GWD \tab
#' Current working directory inside Genomic Database\cr GTRACKS \tab List of
#' all available tracks\cr GINTERVS \tab List of all available intervals\cr
#' GVTRACKS \tab List of all available virtual tracks\cr ALLGENOME \tab List of
#' all chromosomes and their sizes\cr GITERATOR.INTERVALS \tab A set of
#' iterator intervals for which the track expression is evaluated\cr }
#'
#' @aliases gdb.init gdb.init.examples gsetroot
#' @param groot the root directory of the Genomic Database
#' @param dir the current working directory inside the Genomic Database
#' @param rescan indicates whether the file structure should be rescanned
#' @return None.
#' @seealso \code{\link{gdb.reload}}, \code{\link{gdb.create}},
#' \code{\link{gdir.cd}}, \code{\link{gtrack.ls}}, \code{\link{gintervals.ls}},
#' \code{\link{gvtrack.ls}}, \code{\link{gset_input_mode}}
#' @keywords ~db ~data ~database
#' @export gdb.init
gdb.init <- function(groot = NULL, dir = NULL, rescan = FALSE) {
    if (is.null(groot)) {
        stop("Usage: gdb.init(groot, dir = NULL, rescan = FALSE)", call. = F)
    }
    gsetroot(groot, dir, rescan)
}

#' @rdname gdb.init
#' @export
gdb.init_examples <- function() {
    gsetroot(system.file("trackdb/test", package = "misha"))
}



#' Reloads database from the disk
#'
#' Reloads database from disk: list of tracks, intervals, etc.
#'
#' Reloads Genomic Database from disk: list of tracks, intervals, etc. Use this
#' function if you manually add tracks or if for any reason the database
#' becomes corrupted. If 'rescan' is 'TRUE', the list of tracks and intervals
#' is achieved by rescanning directory structure under the current current
#' working directory. Otherwise 'gdb.reload' attempts to use the cached list
#' that resides in 'GROOT/.db.cache' file.
#'
#' @param rescan indicates whether the file structure should be rescanned
#' @seealso \code{\link{gdb.init}}, \code{\link{gdb.create}},
#' \code{\link{gdir.cd}}, \code{\link{gset_input_mode}}
#' @keywords ~db
#' @export gdb.reload
gdb.reload <- function(rescan = TRUE) {
    if (!exists("GROOT")) {
        stop("gdb.init() must be called beforehand.", call. = F)
    }

    if (.ggetOption(".gautocompletion", FALSE)) {
        .gundefine_autocompletion_vars()
    }

    assign("GTRACKS", NULL, envir = .GlobalEnv)
    assign("GINTERVS", NULL, envir = .GlobalEnv)

    dir <- get("GWD")

    res <- ""

    if (get("GWD") != paste(get("GROOT"), "tracks", sep = "/")) {
        rescan <- TRUE
    }

    db.filename <- paste(get("GROOT"), ".db.cache", sep = "/")

    options(warn = -1) # disable warnings since dir() on non dir or non existing dir produces warnings
    if (!rescan) {
        retv <- try(
            {
                f <- file(db.filename, "rb")
                res <- unserialize(f)
                close(f)
            },
            silent = T
        )

        if (inherits(retv, "try-error")) {
            rescan <- TRUE
        }
    }

    if (rescan) {
        res <- .gcall("gfind_tracks_n_intervals", dir, new.env(parent = parent.frame()), silent = TRUE)
        if (get("GWD") == paste(get("GROOT"), "tracks", sep = "/")) {
            try(
                {
                    f <- file(db.filename, "wb")
                    serialize(res, f)
                    close(f)
                },
                silent = T
            )
        } else {
            unlink(db.filename, recursive = TRUE)
        }
    }
    options(warn = 0) # restore the warning behavior

    tracks <- res[[1]]
    intervals <- res[[2]]

    tracks <- sort(tracks)
    intervals <- sort(intervals)

    res <- intersect(tracks, intervals)
    if (length(res) > 0) {
        stop("The following tracks exist also as intervals: ", paste(res, collapse = " "))
    }

    if (.ggetOption(".gautocompletion", FALSE)) {
        .gdefine_autocompletion_vars(tracks, intervals, gvtrack.ls(), .ggetOption(".ginteractive", FALSE))
    }

    assign("GTRACKS", tracks, envir = .GlobalEnv)
    assign("GINTERVS", intervals, envir = .GlobalEnv)
}


.gdb.convert_attrs <- function() {
    .gcheckroot()

    ro_attrs <- c("created.by", "created.date")
    .gcall_noninteractive(gdb.set_readonly_attrs, ro_attrs)

    for (track in GTRACKS) {
        for (attr in ro_attrs) {
            try(
                {
                    if (.gcall_noninteractive(.gtrack.var.exists, track, attr)) {
                        .gcall_noninteractive(.gtrack.attr.set, track, attr, as.character(.gtrack.var.get(track, attr))[1], T)
                        .gcall_noninteractive(gtrack.var.rm, track, attr)
                    }
                },
                silent = TRUE
            )
        }
        cat(sprintf("%s\n", track))
    }
}

.gdb.convert_tracks <- function() {
    .gcheckroot()

    for (track in GTRACKS) {
        try(
            {
                retv <- try(.gcall_noninteractive(gtrack.info, track), silent = TRUE)
                if (inherits(retv, "try-error") & length(grep("obsolete", retv)) > 0) {
                    cat(sprintf("Converting track %s\n", track))
                    .gcall_noninteractive(gtrack.convert, track)
                }
            },
            silent = TRUE
        )
    }
}


.gconfirmtrackcreate <- function(track) {
    if (!is.na(match(track, get("GTRACKS")))) {
        stop(sprintf("Track %s already exists", track), call. = F)
    }

    path <- gsub(".", "/", track, fixed = T)
    dir <- dirname(path)
    fulldir <- paste(get("GWD"), dir, sep = "/")
    fullpath <- sprintf("%s.track", paste(get("GWD"), path, sep = "/"))

    if (!file.exists(fulldir)) {
        stop(sprintf("Directory %s does not exist", dir), call. = F)
    }

    if (file.exists(fullpath)) {
        stop(sprintf("File %s already exists", path), call. = F)
    }

    if (!is.na(match(track, get("GINTERVS")))) {
        stop(sprintf("Interval %s already exists", track), call. = F)
    }

    if (!is.na(match(track, gvtrack.ls()))) {
        stop(sprintf("Virtual track %s already exists", track), call. = F)
    }

    if (.ggetOption(".gautocompletion", FALSE) && exists(track)) {
        stop(sprintf("Variable \"%s\" shadows the name of the new track.\nPlease remove this variable from the environment or switch off autocompletion mode.", track), call. = F)
    }
}

.gdb.add_track <- function(track) {
    .gcheckroot()

    trackdir <- sprintf("%s.track", paste(get("GWD"), gsub("\\.", "/", track), sep = "/"))
    if (file.exists(trackdir)) {
        tracks <- sort(c(get("GTRACKS"), track))
        intervals <- sort(get("GINTERVS"))

        res <- intersect(tracks, intervals)
        if (length(res) > 0) {
            stop("The following tracks exist also as intervals: ", paste(res, collapse = " "))
        }

        if (.ggetOption(".gautocompletion", FALSE)) {
            if (exists(track, envir = .GlobalEnv)) {
                stop(sprintf("Variable \"%s\" shadows the name of identically named track.\nPlease remove this variable from the environment or switch off autocompletion mode.", track), call. = F)
            }

            if (.ggetOption(".ginteractive", FALSE)) { # set track to NULL otherwise evaluation of track expression pmin(track, 2) will produce a string "2"
                assign(track, NULL, envir = .GlobalEnv)
            } else {
                assign(track, track, envir = .GlobalEnv)
            }
        }

        assign("GTRACKS", tracks, envir = .GlobalEnv)
    }
}

.gdb.rm_track <- function(track) {
    .gcheckroot()

    trackdir <- sprintf("%s.track", paste(get("GWD"), gsub("\\.", "/", track), sep = "/"))
    if (!file.exists(trackdir)) {
        if (.ggetOption(".gautocompletion", FALSE)) {
            if (exists(track, envir = .GlobalEnv)) {
                remove(list = track, envir = .GlobalEnv)
            }
        }

        tracks <- get("GTRACKS")
        tracks <- tracks[tracks != track]
        assign("GTRACKS", tracks, envir = .GlobalEnv)
    }
}

.gdb.add_intervals.set <- function(intervals.set) {
    .gcheckroot()

    fname <- sprintf("%s.interv", paste(get("GWD"), gsub("\\.", "/", intervals.set), sep = "/"))
    if (file.exists(fname)) {
        tracks <- get("GTRACKS")
        intervals <- sort(c(get("GINTERVS"), intervals.set))

        res <- intersect(tracks, intervals)
        if (length(res) > 0) {
            stop("The following tracks exist also as intervals: ", paste(res, collapse = " "))
        }

        if (.ggetOption(".gautocompletion", FALSE)) {
            if (exists(intervals.set, envir = .GlobalEnv)) {
                stop(sprintf("Variable \"%s\" shadows the name of identically named intervals set.\nPlease remove this variable from the environment or switch off autocompletion mode.", intervals.set), call. = F)
            }

            assign(intervals.set, intervals.set, envir = .GlobalEnv)
        }

        assign("GINTERVS", intervals, envir = .GlobalEnv)
    }
}

.gdb.rm_intervals.set <- function(intervals.set) {
    .gcheckroot()

    fname <- sprintf("%s.interv", paste(get("GWD"), gsub("\\.", "/", intervals.set), sep = "/"))
    if (!file.exists(fname)) {
        if (.ggetOption(".gautocompletion", FALSE)) {
            if (exists(intervals.set, envir = .GlobalEnv)) {
                remove(list = intervals.set, envir = .GlobalEnv)
            }
        }

        intervals <- get("GINTERVS")
        intervals <- intervals[intervals != intervals.set]
        assign("GINTERVS", intervals, envir = .GlobalEnv)
    }
}