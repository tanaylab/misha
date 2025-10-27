.gchroms <- function(chroms) {
    if (!is.character(chroms)) {
        chroms <- as.character(chroms)
    }


    if (exists("CHROM_ALIAS", envir = .misha, inherits = FALSE)) {
        alias_map <- get("CHROM_ALIAS", envir = .misha)
        if (length(alias_map)) {
            mapped <- alias_map[chroms]
            matched <- !is.na(mapped)
            if (any(matched)) {
                chroms[matched] <- mapped[matched]
            }
        }
    }

    allchroms <- get("ALLGENOME", envir = .misha)[[1]]$chrom
    indices <- match(chroms, allchroms)

    err.chroms <- chroms[is.na(indices)]
    if (length(err.chroms) > 0) {
        stop(sprintf("Chromosome %s does not exist in the database", err.chroms[1]))
    }
    allchroms[indices]
}

.gnormalize_chrom_names <- function(intervals) {
    if (!is.data.frame(intervals)) {
        return(intervals)
    }

    normalize_col <- function(df, col) {
        if (col %in% colnames(df)) {
            df[[col]] <- .gchroms(as.character(df[[col]]))
        }
        df
    }

    intervals <- normalize_col(intervals, "chrom")
    intervals <- normalize_col(intervals, "chrom1")
    intervals <- normalize_col(intervals, "chrom2")
    intervals
}

.is_per_chromosome_db <- function(groot, chromsizes) {
    if (file.exists(file.path(groot, "seq", "genome.idx"))) {
        return(FALSE)
    }

    if (!nrow(chromsizes)) {
        return(FALSE)
    }

    names_chr <- chromsizes$chrom
    no_prefix <- names_chr[!startsWith(names_chr, "chr")]

    if (!length(no_prefix)) {
        return(FALSE)
    }

    if (length(no_prefix) < 0.8 * length(names_chr)) {
        return(FALSE)
    }

    seq_dir <- file.path(groot, "seq")
    sample_names <- head(no_prefix, 5)
    prefixed_exist <- vapply(sample_names, function(name) {
        file.exists(file.path(seq_dir, paste0("chr", name, ".seq")))
    }, logical(1))

    if (!all(prefixed_exist)) {
        return(FALSE)
    }

    unprefixed_exist <- vapply(sample_names, function(name) {
        file.exists(file.path(seq_dir, paste0(name, ".seq")))
    }, logical(1))

    if (any(unprefixed_exist)) {
        return(FALSE)
    }

    TRUE
}

.compute_chrom_aliases <- function(chroms) {
    chroms <- unique(as.character(chroms))
    alias <- stats::setNames(chroms, chroms)

    maybe_add <- function(name, target) {
        if (!nzchar(name)) {
            return()
        }
        if (name %in% names(alias)) {
            return()
        }
        if (name %in% chroms) {
            return()
        }
        alias[[name]] <<- target
    }

    for (chrom in chroms) {
        unprefixed <- sub("^chr", "", chrom, perl = TRUE)
        prefixed <- paste0("chr", unprefixed)

        if (!identical(unprefixed, chrom)) {
            maybe_add(unprefixed, chrom)
        }

        if (!identical(prefixed, chrom)) {
            maybe_add(prefixed, chrom)
        }

        upper_unprefixed <- toupper(unprefixed)
        upper_chrom <- toupper(chrom)

        if (upper_unprefixed %in% c("M", "MT") || upper_chrom %in% c("CHRM", "CHRMT")) {
            maybe_add("M", chrom)
            maybe_add("MT", chrom)
            maybe_add("chrM", chrom)
        }
    }

    alias
}

.store_chrom_aliases <- function(chroms) {
    alias_map <- .compute_chrom_aliases(chroms)
    assign("CHROM_ALIAS", alias_map, envir = .misha)
}

# Helper function for lazy 2D genome generation
.generate_2d_on_demand <- function(intervals, mode = "full") {
    if (mode == "diagonal") {
        # Only intra-chromosomal pairs (chrom1 == chrom2)
        intervals2d <- cbind(intervals, intervals)
        names(intervals2d) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2")
    } else if (mode == "full") {
        # Full cartesian product of all chromosome pairs
        cartesian <- expand.grid(1:nrow(intervals), 1:nrow(intervals))
        intervals2d <- cbind(intervals[cartesian[, 2], ], intervals[cartesian[, 1], ])
        names(intervals2d) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2")
    } else {
        stop("Unknown 2D generation mode: ", mode, ". Must be 'diagonal' or 'full'")
    }

    rownames(intervals2d) <- 1:nrow(intervals2d)
    intervals2d
}

# Helper function to create deferred 2D placeholder
.create_deferred_2d <- function(n_contigs) {
    intervals2d <- list() # Use empty list instead of NULL so we can set attributes
    attr(intervals2d, "deferred") <- TRUE
    attr(intervals2d, "n_contigs") <- n_contigs
    intervals2d
}

# Helper function to check if 2D is deferred
.is_2d_deferred <- function(intervals2d) {
    is.null(intervals2d) || !is.null(attr(intervals2d, "deferred"))
}


.gcheckroot <- function() {
    if (!exists("GROOT", envir = .misha) || !exists("ALLGENOME", envir = .misha) || is.null(get("GROOT", envir = .misha)) || is.null(get("ALLGENOME", envir = .misha))) {
        stop("Database root directory is not set. Please call gdb.init().", call. = FALSE)
    }
}


.gdir.cd <- function(dir, rescan) {
    oldwd <- getwd()
    on.exit(setwd(oldwd), add = TRUE)
    setwd(get("GWD", envir = .misha))
    tryCatch(
        {
            t <- .gfindtrackinpath(dir)
            if (!is.null(t)) {
                stop(sprintf("Directory %s belongs to track %s", dir, t), call. = FALSE)
            }

            setwd(dir)
            newwd <- getwd()

            assign("GWD", newwd, envir = .misha)
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
        stop("Usage: gsetroot(groot, dir = NULL, rescan = FALSE)", call. = FALSE)
    }

    groot <- normalizePath(groot)

    assign("ALLGENOME", NULL, envir = .misha)
    assign("GROOT", NULL, envir = .misha)
    assign("CHROM_ALIAS", NULL, envir = .misha)

    chromsizes <- read.csv(
        paste(groot, "chrom_sizes.txt", sep = "/"),
        sep = "\t",
        header = FALSE,
        col.names = c("chrom", "size"),
        colClasses = c("character", "numeric")
    )

    is_per_chromosome <- .is_per_chromosome_db(groot, chromsizes)
    assign("DB_IS_PER_CHROMOSOME", is_per_chromosome, envir = .misha)

    canonical_names <- chromsizes$chrom

    # For per-chromosome databases, add "chr" prefix to match seq file names
    if (is_per_chromosome) {
        # Add chr prefix to names that don't have it
        needs_prefix <- !startsWith(canonical_names, "chr")
        canonical_names[needs_prefix] <- paste0("chr", canonical_names[needs_prefix])
    }

    # Always compute chromosome aliases for better usability
    # This allows users to use both "chr1" and "1" interchangeably
    # and ensures aliases remain available after database conversion
    alias_map <- .compute_chrom_aliases(canonical_names)

    intervals <- data.frame(
        chrom = canonical_names,
        start = 0,
        end = as.numeric(chromsizes$size)
    )
    intervals$chrom <- as.factor(intervals$chrom)

    # Validate genome.idx if it exists (indexed format)
    idx_path <- file.path(groot, "seq", "genome.idx")
    if (file.exists(idx_path)) {
        .gcall("gseq_validate_index", file.path(groot, "seq"), .misha_env())
    }

    if (nrow(intervals) == 0) {
        stop("chrom_sizes.txt file does not contain any chromosomes", call. = FALSE)
    }

    for (chrom in intervals$chrom) {
        if (length(grep(sprintf("^%s$", chrom), intervals$chrom)) > 1) {
            stop(sprintf("Chromosome \"%s\" appears more than once in chrom_sizes.txt", chrom))
        }
    }
    intervals <- intervals[order(intervals$chrom), ]
    rownames(intervals) <- 1:nrow(intervals)

    # Lazy 2D generation: only materialize for small genomes
    contig_threshold <- getOption("gmulticontig.2d.threshold", 100)

    if (nrow(intervals) <= contig_threshold) {
        # Small genome: materialize full 2D grid
        cartesian <- expand.grid(1:nrow(intervals), 1:nrow(intervals))
        intervals2d <- cbind(intervals[cartesian[, 2], ], intervals[cartesian[, 1], ])
        names(intervals2d) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2")
        rownames(intervals2d) <- 1:nrow(intervals2d)
    } else {
        # Large genome: defer 2D generation
        intervals2d <- .create_deferred_2d(nrow(intervals))
        message(sprintf(
            "Deferring 2D genome generation (%d contigs > threshold %d). Use gintervals.2d() for specific pairs.",
            nrow(intervals), contig_threshold
        ))
    }

    assign("ALLGENOME", list(intervals, intervals2d), envir = .misha)
    assign("GROOT", groot, envir = .misha)
    assign("GWD", groot, envir = .misha)
    assign("CHROM_ALIAS", alias_map, envir = .misha)

    success <- FALSE
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
            success <- TRUE
        },
        finally = {
            if (!success) {
                assign("ALLGENOME", NULL, envir = .misha)
                assign("GROOT", NULL, envir = .misha)
                assign("GWD", NULL, envir = .misha)
                assign("CHROM_ALIAS", NULL, envir = .misha)
                assign("DB_IS_PER_CHROMOSOME", NULL, envir = .misha)
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
gdb.set_readonly_attrs <- function(attrs) {
    .gcheckroot()

    filename <- paste(get("GROOT", envir = .misha), ".ro_attributes", sep = "/")

    if (is.null(attrs)) {
        unlink(filename)
    } else {
        attrs <- as.character(attrs)

        idx <- which(duplicated(attrs))[1]
        if (!is.na(idx)) {
            stop(sprintf("Attribute %s appears more than once", attrs[idx]), call. = FALSE)
        }

        idx <- which(attrs == "")[1]
        if (!is.na(idx)) {
            stop("Attribute name cannot be an empty string", call. = FALSE)
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
#' Two database formats are supported:
#' \itemize{
#'   \item \strong{indexed}: Single genome.seq + genome.idx (default). Recommended for
#'         genomes with many contigs. Provides better performance and scalability.
#'   \item \strong{per-chromosome}: Separate .seq file per contig.
#' }
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
#' See the 'Genomes' vignette for details on how to create a database from common
#' genome sources.
#'
#' @param groot path to newly created database
#' @param fasta an array of names or URLs of FASTA files. Can contain wildcards
#' for multiple files
#' @param genes.file name or URL of file that contains genes. If 'NULL' no
#' genes are imported
#' @param annots.file name of URL file that contains annotations. If 'NULL' no
#' annotations are imported
#' @param annots.names annotations names
#' @param format database format: "indexed" (default, single genome.seq + genome.idx)
#' or "per-chromosome" (separate .seq file per contig). If NULL, uses the value from
#' \code{getOption("gmulticontig.indexed_format", TRUE)}
#' @param verbose if TRUE, prints verbose messages
#' @return None.
#' @seealso \code{\link{gdb.init}}, \code{\link{gdb.reload}},
#' \code{\link{gintervals.import_genes}}
#' @keywords ~database ~create ~genes
#' @examples
#' \donttest{
#' ftp <- "ftp://hgdownload.soe.ucsc.edu/goldenPath/mm10"
#' mm10_dir <- file.path(tempdir(), "mm10")
#' # only a single chromosome is loaded in this example
#' # see "Genomes" vignette how to download all of them and how
#' # to download other genomes
#' gdb.create(
#'     mm10_dir,
#'     paste(ftp, "chromosomes", paste0(
#'         "chr", c("X"),
#'         ".fa.gz"
#'     ), sep = "/"),
#'     paste(ftp, "database/knownGene.txt.gz", sep = "/"),
#'     paste(ftp, "database/kgXref.txt.gz", sep = "/"),
#'     c(
#'         "kgID", "mRNA", "spID", "spDisplayID", "geneSymbol",
#'         "refseq", "protAcc", "description", "rfamAcc",
#'         "tRnaName"
#'     )
#' )
#' gdb.init(mm10_dir)
#' gintervals.ls()
#' gintervals.all()
#' }
#'
#' @export gdb.create
gdb.create <- function(groot = NULL, fasta = NULL, genes.file = NULL, annots.file = NULL, annots.names = NULL, format = NULL, verbose = FALSE) {
    if (is.null(groot) || is.null(fasta)) {
        stop("Usage: gdb.create(groot, fasta, genes.file = NULL, annots.file = NULL, annots.names = NULL, format = 'indexed', verbose = FALSE)", call. = FALSE)
    }

    # Determine format: parameter overrides option
    if (is.null(format)) {
        # Check for per-chromosome option for backward compatibility
        use_indexed_format <- getOption("gmulticontig.indexed_format", TRUE)
        format <- if (use_indexed_format) "indexed" else "per-chromosome"
    } else {
        # Validate format parameter
        format <- match.arg(format, c("indexed", "per-chromosome"))
    }
    use_indexed_format <- (format == "indexed")

    if (file.exists(groot)) {
        stop(sprintf("Directory %s already exists", groot), call. = FALSE)
    }

    success <- FALSE
    allgenome.old <- NULL
    groot.old <- NULL
    chrom_alias.old <- NULL
    db_per_chromosome.old <- NULL
    if (exists("ALLGENOME", envir = .misha)) {
        allgenome.old <- get("ALLGENOME", envir = .misha)
    }
    if (exists("GROOT", envir = .misha)) {
        groot.old <- get("GROOT", envir = .misha)
    }
    if (exists("CHROM_ALIAS", envir = .misha)) {
        chrom_alias.old <- get("CHROM_ALIAS", envir = .misha)
    }
    if (exists("DB_IS_PER_CHROMOSOME", envir = .misha)) {
        db_per_chromosome.old <- get("DB_IS_PER_CHROMOSOME", envir = .misha)
    }

    tryCatch(
        {
            assign("CHROM_ALIAS", NULL, envir = .misha)
            assign("DB_IS_PER_CHROMOSOME", FALSE, envir = .misha)
            dir.create(groot, showWarnings = FALSE, recursive = TRUE, mode = "0777")
            dir.create(paste(groot, "pssms", sep = "/"), showWarnings = FALSE, recursive = TRUE, mode = "0777")
            dir.create(paste(groot, "seq", sep = "/"), showWarnings = FALSE, recursive = TRUE, mode = "0777")
            dir.create(paste(groot, "tracks", sep = "/"), showWarnings = FALSE, recursive = TRUE, mode = "0777")

            # Detect import mode: multi-FASTA (indexed) vs per-chromosome
            is_single_file <- (length(fasta) == 1 && file.exists(fasta))

            chroms <- NULL
            chrom.sizes <- NULL

            if (use_indexed_format && is_single_file) {
                # Multi-FASTA import to indexed format
                # C++ function parses FASTA, creates genome.seq and genome.idx,
                # and returns a data frame with columns: name, size (sorted alphabetically)
                if (verbose) message("Creating indexed genome format...")
                contig_info <- .gseq.import_multifasta(groot, fasta, verbose)

                if (nrow(contig_info) == 0) {
                    stop("No contigs were imported from multi-FASTA file", call. = FALSE)
                }

                # Use the data frame directly (already sorted by name, matching ALLGENOME order)
                chrom.sizes <- data.frame(chrom = contig_info$name, size = contig_info$size)
            } else {
                # Per-chromosome import
                if (verbose) message("Creating per-chromosome genome format...")
                chroms <- .gseq.import(groot, fasta)

                if (!length(chroms)) {
                    stop("No FASTA files were imported", call. = FALSE)
                }

                # Use chromosome names as-is (no chr prefix)
                seq.files <- paste(chroms, ".seq", sep = "")
                seq.files <- paste(paste(groot, "seq", sep = "/"), seq.files, sep = "/")
                chrom.sizes <- data.frame(chrom = chroms, size = file.info(seq.files)$size)
            }

            utils::write.table(chrom.sizes, paste(groot, "chrom_sizes.txt", sep = "/"), quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)

            # before calling gintervals.import_genes new ALLGENOME must be set
            intervals <- data.frame(
                chrom = as.factor(as.character(chrom.sizes$chrom)), # Preserve original names
                start = 0, end = as.numeric(chrom.sizes$size)
            )
            intervals <- intervals[order(intervals$chrom), ]
            rownames(intervals) <- 1:nrow(intervals)

            # Lazy 2D generation: only materialize for small genomes
            contig_threshold <- getOption("gmulticontig.2d.threshold", 100)

            if (nrow(intervals) <= contig_threshold) {
                # Small genome: materialize full 2D grid
                cartesian <- expand.grid(1:nrow(intervals), 1:nrow(intervals))
                intervals2d <- cbind(intervals[cartesian[, 2], ], intervals[cartesian[, 1], ])
                names(intervals2d) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2")
                rownames(intervals2d) <- 1:nrow(intervals2d)
            } else {
                # Large genome: defer 2D generation
                intervals2d <- .create_deferred_2d(nrow(intervals))
                if (verbose) {
                    message(sprintf(
                        "Deferring 2D genome generation (%d contigs > threshold %d). Use gintervals.2d() for specific pairs.",
                        nrow(intervals), contig_threshold
                    ))
                }
            }

            assign("ALLGENOME", list(intervals, intervals2d), envir = .misha)
            assign("GROOT", groot, envir = .misha)
            .store_chrom_aliases(levels(intervals$chrom))

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
            serialize(c("created.by", "created.date", "created.user"), f)
            close(f)

            if (verbose) message("Database was successfully created")
            success <- TRUE
        },
        finally = {
            assign("ALLGENOME", allgenome.old, envir = .misha)
            assign("GROOT", groot.old, envir = .misha)
            assign("CHROM_ALIAS", chrom_alias.old, envir = .misha)
            assign("DB_IS_PER_CHROMOSOME", db_per_chromosome.old, envir = .misha)
            if (!success) {
                unlink(groot, recursive = TRUE)
            }
        }
    )
    retv <- 0 # suppress return value
}

#' Create and Load a Genome Database
#'
#' This function downloads, extracts, and loads a misha genome database for the specified genome.
#'
#' @param genome A character string specifying the genome to download. Supported genomes are "mm9", "mm10", "mm39", "hg19", and "hg38".
#' @param path A character string specifying the directory where the genome will be extracted. Defaults to genome name (e.g. "mm10") in the current working directory.
#' @param tmpdir A character string specifying the directory for storing temporary files. This is used for storing the downloaded genome file.
#'
#' @details
#' The function checks if the specified genome is available. If tmpdir, it constructs the download URL, downloads the genome file,
#' extracts it to the specified directory, and loads the genome database using \code{gsetroot}. The function also calls \code{gdb.reload} to reload the genome database.
#'
#' @return None.
#'
#' @examples
#' \donttest{
#' mm10_dir <- tempdir()
#' gdb.create_genome("mm10", path = mm10_dir)
#' list.files(file.path(mm10_dir, "mm10"))
#' gsetroot(file.path(mm10_dir, "mm10"))
#' gintervals.ls()
#' }
#'
#' @export
gdb.create_genome <- function(genome, path = getwd(), tmpdir = tempdir()) {
    supported_genomes <- c("mm9", "mm10", "mm39", "hg19", "hg38")

    if (!genome %in% supported_genomes) {
        stop(paste("The genome", genome, "is not available yet. Available genomes are:", paste(supported_genomes, collapse = ", ")))
    }

    base_url <- "https://misha-genome.s3.eu-west-1.amazonaws.com/"
    url <- paste0(base_url, genome, ".tar.gz")

    temp_file <- tempfile(fileext = ".tar.gz", tmpdir = tmpdir)

    withr::local_options(list(timeout = 60 * 60 * 2))
    message("Downloading ", genome, " genome...")
    utils::download.file(url, temp_file, mode = "wb")

    dir.create(path, showWarnings = FALSE, recursive = TRUE)

    message("Extracting ", genome, " genome...")
    utils::untar(temp_file, exdir = path)

    unlink(temp_file)

    message("Loading misha root...")
    gdb.init(file.path(path, genome), rescan = TRUE)

    message(genome, " genome has been successfully downloaded and extracted to ", file.path(path, genome))
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

    filename <- paste(get("GROOT", envir = .misha), ".ro_attributes", sep = "/")
    attrs <- NULL
    if (file.exists(filename)) {
        f <- file(filename, "rb")
        attrs <- unserialize(f)
        close(f)
        if (!is.character(attrs)) {
            stop(sprintf("Invalid format of read-only atrributes file %s", filename), call. = FALSE)
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
#' few variables are defined at an environment called \code{.misha}, and can be
#' accessed using \code{.misha$variable}, e.g. \code{.misha$ALLGENOME}.
#' These variables should not be modified by user.
#'
#' \tabular{ll}{ GROOT \tab Root directory of Genomic Database\cr GWD \tab
#' Current working directory inside Genomic Database\cr GTRACKS \tab List of
#' all available tracks\cr GINTERVS \tab List of all available intervals\cr
#' GVTRACKS \tab List of all available virtual tracks\cr ALLGENOME \tab List of
#' all chromosomes and their sizes\cr GITERATOR.INTERVALS \tab A set of
#' iterator intervals for which the track expression is evaluated\cr }
#'
#' When option 'gmulticontig.indexed_format' is set to TRUE, the function
#' loads a database with "indexed" track format.
#'
#' @aliases gdb.init gdb.init.examples gsetroot
#' @param groot the root directory of the Genomic Database
#' @param dir the current working directory inside the Genomic Database
#' @param rescan indicates whether the file structure should be rescanned
#' @return None.
#' @seealso \code{\link{gdb.reload}}, \code{\link{gdb.create}},
#' \code{\link{gdir.cd}}, \code{\link{gtrack.ls}}, \code{\link{gintervals.ls}},
#' \code{\link{gvtrack.ls}}
#' @keywords ~db ~data ~database
#' @export gdb.init
gdb.init <- function(groot = NULL, dir = NULL, rescan = FALSE) {
    if (is.null(groot)) {
        stop("Usage: gdb.init(groot, dir = NULL, rescan = FALSE)", call. = FALSE)
    }
    gsetroot(groot, dir, rescan)
}

#' @rdname gdb.init
#' @export
gdb.init_examples <- function() {
    db_dir <- tempdir()
    test_path <- file.path(db_dir, "trackdb/test")
    unlink(test_path, recursive = TRUE)
    utils::untar(system.file("testdb.tar.gz", package = "misha"), exdir = db_dir)
    gsetroot(test_path)
    if (getOption("gmulticontig.indexed_format", FALSE)) {
        gdb.convert_to_indexed(test_path, convert_tracks = TRUE, remove_old_files = TRUE, convert_intervals = TRUE, verbose = FALSE, force = TRUE)
        gsetroot(test_path)
    }
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
#' \code{\link{gdir.cd}},
#' @keywords ~db
#' @return No return value, called for side effects.
#' @export gdb.reload
gdb.reload <- function(rescan = TRUE) {
    if (!exists("GROOT", envir = .misha)) {
        stop("gdb.init() must be called beforehand.", call. = FALSE)
    }

    assign("GTRACKS", NULL, envir = .misha)
    assign("GINTERVS", NULL, envir = .misha)

    dir <- get("GWD", envir = .misha)

    res <- ""

    if (get("GWD", envir = .misha) != paste(get("GROOT", envir = .misha), "tracks", sep = "/")) {
        rescan <- TRUE
    }

    db.filename <- paste(get("GROOT", envir = .misha), ".db.cache", sep = "/")

    suppressWarnings({ # disable warnings since dir() on non dir or non existing dir produces warnings
        if (!rescan) {
            retv <- try(
                {
                    f <- file(db.filename, "rb")
                    res <- unserialize(f)
                    close(f)
                },
                silent = TRUE
            )

            if (inherits(retv, "try-error")) {
                rescan <- TRUE
            }
        }

        if (rescan) {
            res <- .gcall("gfind_tracks_n_intervals", dir, .misha_env())
            if (get("GWD", envir = .misha) == paste(get("GROOT", envir = .misha), "tracks", sep = "/")) {
                try(
                    {
                        f <- file(db.filename, "wb")
                        serialize(res, f)
                        close(f)
                    },
                    silent = TRUE
                )
            } else {
                unlink(db.filename, recursive = TRUE)
            }
        }
    })

    tracks <- res[[1]]
    intervals <- res[[2]]

    tracks <- sort(tracks)
    intervals <- sort(intervals)

    res <- intersect(tracks, intervals)
    if (length(res) > 0) {
        stop("The following tracks exist also as intervals: ", paste(res, collapse = " "))
    }

    assign("GTRACKS", tracks, envir = .misha)
    assign("GINTERVS", intervals, envir = .misha)
}


.gdb.convert_attrs <- function() {
    .gcheckroot()

    ro_attrs <- c("created.by", "created.date", "created.user")
    .gcall_noninteractive(gdb.set_readonly_attrs, ro_attrs)

    for (track in .misha$GTRACKS) {
        for (attr in ro_attrs) {
            try(
                {
                    if (.gcall_noninteractive(.gtrack.var.exists, track, attr)) {
                        .gcall_noninteractive(.gtrack.attr.set, track, attr, as.character(.gtrack.var.get(track, attr))[1], TRUE)
                        .gcall_noninteractive(gtrack.var.rm, track, attr)
                    }
                },
                silent = TRUE
            )
        }
        message(track)
    }
}

.gdb.convert_tracks <- function() {
    .gcheckroot()

    for (track in .misha$GTRACKS) {
        try(
            {
                retv <- try(.gcall_noninteractive(gtrack.info, track))
                if (inherits(retv, "try-error") & length(grep("obsolete", retv)) > 0) {
                    message(sprintf("Converting track %s", track))
                    .gcall_noninteractive(gtrack.convert, track)
                }
            },
            silent = TRUE
        )
    }
}


.gconfirmtrackcreate <- function(track) {
    if (!is.na(match(track, get("GTRACKS", envir = .misha)))) {
        stop(sprintf("Track %s already exists", track), call. = FALSE)
    }

    path <- gsub(".", "/", track, fixed = TRUE)
    dir <- dirname(path)
    fulldir <- paste(get("GWD", envir = .misha), dir, sep = "/")
    fullpath <- sprintf("%s.track", paste(get("GWD", envir = .misha), path, sep = "/"))

    if (!file.exists(fulldir)) {
        stop(sprintf("Directory %s does not exist", dir), call. = FALSE)
    }

    if (file.exists(fullpath)) {
        stop(sprintf("File %s already exists", path), call. = FALSE)
    }

    if (!is.na(match(track, get("GINTERVS", envir = .misha)))) {
        stop(sprintf("Interval %s already exists", track), call. = FALSE)
    }

    if (!is.na(match(track, gvtrack.ls()))) {
        stop(sprintf("Virtual track %s already exists", track), call. = FALSE)
    }

    if (.ggetOption(".gautocompletion", FALSE) && exists(track)) {
        stop(sprintf("Variable \"%s\" shadows the name of the new track.\nPlease remove this variable from the environment or switch off autocompletion mode.", track), call. = FALSE)
    }
}

.gdb.add_track <- function(track) {
    .gcheckroot()

    trackdir <- sprintf("%s.track", paste(get("GWD", envir = .misha), gsub("\\.", "/", track), sep = "/"))
    if (file.exists(trackdir)) {
        tracks <- sort(c(get("GTRACKS", envir = .misha), track))
        intervals <- sort(get("GINTERVS", envir = .misha))

        res <- intersect(tracks, intervals)
        if (length(res) > 0) {
            stop("The following tracks exist also as intervals: ", paste(res, collapse = " "))
        }

        if (.ggetOption(".gautocompletion", FALSE)) {
            if (exists(track, envir = .misha)) {
                stop(sprintf("Variable \"%s\" shadows the name of identically named track.\nPlease remove this variable from the environment or switch off autocompletion mode.", track), call. = FALSE)
            }

            if (.ggetOption(".ginteractive", FALSE)) { # set track to NULL otherwise evaluation of track expression pmin(track, 2) will produce a string "2"
                assign(track, NULL, envir = .misha)
            } else {
                assign(track, track, envir = .misha)
            }
        }

        assign("GTRACKS", tracks, envir = .misha)
    }
}

.gdb.rm_track <- function(track) {
    .gcheckroot()

    trackdir <- sprintf("%s.track", paste(get("GWD", envir = .misha), gsub("\\.", "/", track), sep = "/"))
    if (!file.exists(trackdir)) {
        if (.ggetOption(".gautocompletion", FALSE)) {
            if (exists(track, envir = .misha)) {
                remove(list = track, envir = .misha)
            }
        }

        tracks <- get("GTRACKS", envir = .misha)
        tracks <- tracks[tracks != track]
        assign("GTRACKS", tracks, envir = .misha)
    }
}

.gdb.add_intervals.set <- function(intervals.set) {
    .gcheckroot()

    fname <- sprintf("%s.interv", paste(get("GWD", envir = .misha), gsub("\\.", "/", intervals.set), sep = "/"))
    if (file.exists(fname)) {
        tracks <- get("GTRACKS", envir = .misha)
        intervals <- sort(c(get("GINTERVS", envir = .misha), intervals.set))

        res <- intersect(tracks, intervals)
        if (length(res) > 0) {
            stop("The following tracks exist also as intervals: ", paste(res, collapse = " "))
        }

        if (.ggetOption(".gautocompletion", FALSE)) {
            if (exists(intervals.set, envir = .misha)) {
                stop(sprintf("Variable \"%s\" shadows the name of identically named intervals set.\nPlease remove this variable from the environment or switch off autocompletion mode.", intervals.set), call. = FALSE)
            }

            assign(intervals.set, intervals.set, envir = .misha)
        }

        assign("GINTERVS", intervals, envir = .misha)
    }
}

.gdb.rm_intervals.set <- function(intervals.set) {
    .gcheckroot()

    fname <- sprintf("%s.interv", paste(get("GWD", envir = .misha), gsub("\\.", "/", intervals.set), sep = "/"))
    if (!file.exists(fname)) {
        if (.ggetOption(".gautocompletion", FALSE)) {
            if (exists(intervals.set, envir = .misha)) {
                remove(list = intervals.set, envir = .misha)
            }
        }

        intervals <- get("GINTERVS", envir = .misha)
        intervals <- intervals[intervals != intervals.set]
        assign("GINTERVS", intervals, envir = .misha)
    }
}
