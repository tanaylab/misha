# Root and initialization functions

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

    # Include both canonical names and aliases in factor levels
    # This allows gintervals.all() to work with both forms
    all_chrom_names <- unique(c(canonical_names, names(alias_map)))

    # Validate genome.idx if it exists (indexed format)
    idx_path <- file.path(groot, "seq", "genome.idx")
    if (file.exists(idx_path)) {
        .gcall("gseq_validate_index", file.path(groot, "seq"), .misha_env())
    }

    intervals <- data.frame(
        chrom = canonical_names,
        start = 0,
        end = as.numeric(chromsizes$size)
    )

    # For indexed databases, preserve the order from chrom_sizes.txt to match genome.idx chromid assignments
    # For per-chromosome databases, sort alphabetically for backward compatibility with existing test snapshots
    if (is_per_chromosome) {
        intervals <- intervals[order(canonical_names), ]
        canonical_names <- sort(canonical_names)
    }

    intervals$chrom <- factor(intervals$chrom, levels = canonical_names)

    if (nrow(intervals) == 0) {
        stop("chrom_sizes.txt file does not contain any chromosomes", call. = FALSE)
    }

    # Check for duplicate chromosomes
    dupes <- intervals$chrom[duplicated(intervals$chrom)]
    if (length(dupes) > 0) {
        stop(sprintf("Chromosome \"%s\" appears more than once in chrom_sizes.txt", dupes[1]), call. = FALSE)
    }

    rownames(intervals) <- 1:nrow(intervals)

    # Lazy 2D generation: only materialize for small genomes
    contig_threshold <- getOption("gmulticontig.2d.threshold", 100)

    if (nrow(intervals) <= contig_threshold) {
        # Small genome: materialize full 2D grid
        cartesian <- expand.grid(1:nrow(intervals), 1:nrow(intervals))
        intervals2d <- cbind(intervals[cartesian[, 2], ], intervals[cartesian[, 1], ])
        names(intervals2d) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2")
        # Ensure chrom1 and chrom2 have the same factor levels as intervals$chrom
        # Use canonical_names only (not aliases) for consistency with 1D intervals
        intervals2d$chrom1 <- factor(intervals2d$chrom1, levels = canonical_names)
        intervals2d$chrom2 <- factor(intervals2d$chrom2, levels = canonical_names)
        rownames(intervals2d) <- 1:nrow(intervals2d)
    } else {
        # Large genome: defer 2D generation
        intervals2d <- .create_deferred_2d(nrow(intervals))
    }

    assign("ALLGENOME", list(intervals, intervals2d), envir = .misha)
    assign("GROOT", groot, envir = .misha)
    assign("GWD", groot, envir = .misha)
    assign("CHROM_ALIAS", alias_map, envir = .misha)

    if (.gdb.cache_is_dirty(groot)) {
        rescan <- TRUE
    }

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
