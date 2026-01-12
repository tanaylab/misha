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

            # Check if directory exists before trying to setwd
            if (!dir.exists(dir)) {
                stop(sprintf("Directory does not exist: %s", dir), call. = FALSE)
            }

            tryCatch(
                {
                    setwd(dir)
                },
                error = function(e) {
                    stop(sprintf("Cannot change to directory '%s': %s", dir, e$message), call. = FALSE)
                }
            )

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

    # Support multiple database paths (multi-db feature)
    # Validate paths exist before calling normalizePath
    for (i in seq_along(groot)) {
        if (!dir.exists(groot[i])) {
            stop(sprintf("Database directory does not exist: %s", groot[i]), call. = FALSE)
        }

        # Check for required directories
        tracks_dir <- file.path(groot[i], "tracks")
        seq_dir <- file.path(groot[i], "seq")

        if (!dir.exists(tracks_dir)) {
            stop(sprintf("Database directory '%s' does not contain a 'tracks' subdirectory. This does not appear to be a valid misha database.", groot[i]), call. = FALSE)
        }

        if (!dir.exists(seq_dir)) {
            stop(sprintf("Database directory '%s' does not contain a 'seq' subdirectory. This does not appear to be a valid misha database.", groot[i]), call. = FALSE)
        }
    }

    groots <- normalizePath(groot, mustWork = TRUE)

    # Clear all state
    assign("ALLGENOME", NULL, envir = .misha)
    assign("GROOT", NULL, envir = .misha)
    assign("GROOTS", NULL, envir = .misha)
    assign("CHROM_ALIAS", NULL, envir = .misha)
    assign("GTRACK_DB", NULL, envir = .misha)
    assign("GINTERVALS_DB", NULL, envir = .misha)

    # Read and validate chrom_sizes from first database
    chrom_sizes_path <- file.path(groots[1], "chrom_sizes.txt")
    if (!file.exists(chrom_sizes_path)) {
        stop(sprintf("Database directory '%s' does not contain a chrom_sizes.txt file. This does not appear to be a valid misha database.", groots[1]), call. = FALSE)
    }

    chromsizes <- tryCatch(
        read.csv(
            chrom_sizes_path,
            sep = "\t",
            header = FALSE,
            col.names = c("chrom", "size"),
            colClasses = c("character", "numeric")
        ),
        error = function(e) {
            stop(sprintf("Failed to read chrom_sizes.txt from '%s': %s", groots[1], e$message), call. = FALSE)
        }
    )

    # Validate all databases have identical chrom_sizes.txt
    if (length(groots) > 1) {
        for (g in groots[-1]) {
            cs_path <- file.path(g, "chrom_sizes.txt")
            if (!file.exists(cs_path)) {
                stop(sprintf("Database directory '%s' does not contain a chrom_sizes.txt file. This does not appear to be a valid misha database.", g), call. = FALSE)
            }

            cs <- tryCatch(
                read.csv(
                    cs_path,
                    sep = "\t",
                    header = FALSE,
                    col.names = c("chrom", "size"),
                    colClasses = c("character", "numeric")
                ),
                error = function(e) {
                    stop(sprintf("Failed to read chrom_sizes.txt from '%s': %s", g, e$message), call. = FALSE)
                }
            )

            if (!identical(chromsizes, cs)) {
                stop(sprintf("All databases must have identical chrom_sizes.txt. Database '%s' has different chromosome sizes.", g), call. = FALSE)
            }
        }
    }

    # Use first database for shared resources (seq, chrom_sizes)
    groot <- groots[1]

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
    assign("GROOTS", groots, envir = .misha)
    assign("GWD", groot, envir = .misha)
    assign("CHROM_ALIAS", alias_map, envir = .misha)

    # Check if any database cache is dirty
    for (g in groots) {
        if (.gdb.cache_is_dirty(g)) {
            rescan <- TRUE
            break
        }
    }

    success <- FALSE
    tryCatch(
        {
            if (is.null(dir)) {
                # Default to last database's tracks (most likely writable for user)
                .gdir.cd(file.path(groots[length(groots)], "tracks"), rescan)
            } else {
                if (nchar(dir) < 1) {
                    stop("dir argument is an empty string", call. = FALSE)
                }

                # Check if user is trying to pass a second database path instead of a directory
                # This is a common mistake when migrating from single-db to multi-db
                if (dir.exists(dir)) {
                    has_tracks <- dir.exists(file.path(dir, "tracks"))
                    has_seq <- dir.exists(file.path(dir, "seq"))
                    has_chrom_sizes <- file.exists(file.path(dir, "chrom_sizes.txt"))

                    if ((has_tracks || has_seq) && has_chrom_sizes) {
                        stop(sprintf(
                            "The 'dir' parameter ('%s') looks like a misha database path.\nTo connect multiple databases, use: gsetroot(c(\"%s\", \"%s\"))",
                            dir, groot[1], dir
                        ), call. = FALSE)
                    }
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
                assign("GROOTS", NULL, envir = .misha)
                assign("GWD", NULL, envir = .misha)
                assign("CHROM_ALIAS", NULL, envir = .misha)
                assign("DB_IS_PER_CHROMOSOME", NULL, envir = .misha)
                assign("GTRACK_DB", NULL, envir = .misha)
                assign("GINTERVALS_DB", NULL, envir = .misha)
            }
        }
    )
}

#' List connected databases
#'
#' Returns a list of all database paths currently connected to the session.
#'
#' When multiple databases are connected using \code{gsetroot(c(db1, db2, ...))},
#' this function returns all database paths in connection order. The first
#' database is used for shared resources (seq, chrom_sizes), and later databases
#' take precedence for track resolution ("last wins").
#'
#' @return Character vector of database paths in connection order.
#' @seealso \code{\link{gsetroot}}, \code{\link{gdb.info}}, \code{\link{gtrack.db}}
#' @keywords ~db ~database
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gdb.ls()
#' @export
gdb.ls <- function() {
    .gcheckroot()
    get("GROOTS", envir = .misha)
}


#' Get summary of connected databases
#'
#' Returns summary information about all connected databases.
#'
#' This function returns a data frame with information about each connected
#' database, including the number of tracks, number of interval sets, and
#' whether the database is writable. Use this when multiple databases are
#' connected to get an overview of track distribution.
#'
#' Note: This is different from \code{\link{gdb.info}} which provides detailed
#' information about a single database's format and structure.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{db}{Character: The database path}
#'   \item{tracks}{Integer: Number of tracks from this database}
#'   \item{intervals}{Integer: Number of interval sets from this database}
#'   \item{writable}{Logical: Whether the database's tracks directory is writable}
#' }
#' @seealso \code{\link{gdb.ls}}, \code{\link{gsetroot}}, \code{\link{gtrack.db}}, \code{\link{gdb.info}}
#' @keywords ~db ~database ~info
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gdb.summary()
#'
#' @export gdb.summary
gdb.summary <- function() {
    .gcheckroot()
    groots <- get("GROOTS", envir = .misha)
    track_db <- get("GTRACK_DB", envir = .misha)
    intervals_db <- get("GINTERVALS_DB", envir = .misha)

    # Count tracks and intervals per database
    track_vals <- unlist(track_db, use.names = FALSE)
    if (length(track_vals)) {
        track_counts <- as.integer(table(factor(track_vals, levels = groots)))
    } else {
        track_counts <- rep(0L, length(groots))
    }

    interval_vals <- unlist(intervals_db, use.names = FALSE)
    if (length(interval_vals)) {
        interval_counts <- as.integer(table(factor(interval_vals, levels = groots)))
    } else {
        interval_counts <- rep(0L, length(groots))
    }

    # Check write permissions
    writable <- vapply(groots, function(g) {
        file.access(file.path(g, "tracks"), 2) == 0
    }, logical(1))

    data.frame(
        db = groots,
        tracks = track_counts,
        intervals = interval_counts,
        writable = writable,
        stringsAsFactors = FALSE,
        row.names = NULL
    )
}


#' Create a linked database with symlinks to a parent database
#'
#' Creates a new database directory structure with symbolic links to the
#' parent database's seq/ directory and chrom_sizes.txt file.
#'
#' This is useful for creating a writable database that shares sequence
#' data with a read-only main database. The new database can then be used
#' together with the parent database via \code{gsetroot(c(parent, linked_db))}.
#'
#' @param path Path for the new linked database
#' @param parent Path to the parent database (with seq and chrom_sizes.txt)
#' @return Invisible TRUE on success
#' @seealso \code{\link{gsetroot}}, \code{\link{gdb.create}}, \code{\link{gdb.ls}}
#' @keywords ~db ~database ~create
#' @examples
#' \dontrun{
#' # Create linked database sharing sequence data with main database
#' gdb.create_linked("~/my_tracks", parent = "/shared/genomics/hg38")
#'
#' # Connect to both databases
#' gsetroot(c("/shared/genomics/hg38", "~/my_tracks"))
#' }
#'
#' @export gdb.create_linked
gdb.create_linked <- function(path, parent) {
    if (missing(path) || missing(parent)) {
        stop("Usage: gdb.create_linked(path, parent)", call. = FALSE)
    }

    parent <- normalizePath(parent, mustWork = TRUE)
    path <- path.expand(path)

    # Validate parent database
    parent_chrom_sizes <- file.path(parent, "chrom_sizes.txt")
    parent_seq <- file.path(parent, "seq")
    if (!file.exists(parent_chrom_sizes)) {
        stop("Parent database missing chrom_sizes.txt", call. = FALSE)
    }
    if (!dir.exists(parent_seq)) {
        stop("Parent database missing seq/ directory", call. = FALSE)
    }

    if (dir.exists(path)) {
        stop(sprintf("Directory '%s' already exists", path), call. = FALSE)
    }

    # Create user database structure
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(path, "tracks"), showWarnings = FALSE)

    # Create symlinks to parent database
    chrom_link <- file.path(path, "chrom_sizes.txt")
    seq_link <- file.path(path, "seq")

    if (!file.symlink(parent_chrom_sizes, chrom_link)) {
        unlink(path, recursive = TRUE)
        stop("Failed to create symbolic link for chrom_sizes.txt", call. = FALSE)
    }

    if (!file.symlink(parent_seq, seq_link)) {
        unlink(path, recursive = TRUE)
        stop("Failed to create symbolic link for seq/", call. = FALSE)
    }

    message(sprintf("Created linked database at %s (linked to %s)", path, parent))
    invisible(TRUE)
}
