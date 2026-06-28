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

    # Dataset API: only accept single database path
    if (length(groot) > 1) {
        stop("gsetroot() accepts a single database path. Use gdataset.load() to load additional datasets.", call. = FALSE)
    }

    # Validate path exists before calling normalizePath
    if (!dir.exists(groot)) {
        stop(sprintf("Database directory does not exist: %s", groot), call. = FALSE)
    }

    # Check for required directories
    tracks_dir <- file.path(groot, "tracks")
    seq_dir <- file.path(groot, "seq")

    if (!dir.exists(tracks_dir)) {
        stop(sprintf("Database directory '%s' does not contain a 'tracks' subdirectory. This does not appear to be a valid misha database.", groot), call. = FALSE)
    }

    if (!dir.exists(seq_dir)) {
        stop(sprintf("Database directory '%s' does not contain a 'seq' subdirectory. This does not appear to be a valid misha database.", groot), call. = FALSE)
    }

    groot <- normalizePath(groot, mustWork = TRUE)

    # Read and validate chrom_sizes BEFORE clearing the current session state, so
    # a missing/malformed chrom_sizes.txt leaves the previously-loaded database
    # intact instead of unloading it on a failed gsetroot.
    chrom_sizes_path <- file.path(groot, "chrom_sizes.txt")
    if (!file.exists(chrom_sizes_path)) {
        stop(sprintf("Database directory '%s' does not contain a chrom_sizes.txt file. This does not appear to be a valid misha database.", groot), call. = FALSE)
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
            stop(sprintf("Failed to read chrom_sizes.txt from '%s': %s", groot, e$message), call. = FALSE)
        }
    )

    if (nrow(chromsizes) == 0) {
        stop("chrom_sizes.txt file does not contain any chromosomes", call. = FALSE)
    }

    # Clear all state
    assign("ALLGENOME", NULL, envir = .misha)
    assign("GROOT", NULL, envir = .misha)
    assign("CHROM_ALIAS", NULL, envir = .misha)
    .refresh_chrom_alias_env()
    assign("GTRACK_DATASET", NULL, envir = .misha)
    assign("GINTERVALS_DATASET", NULL, envir = .misha)
    assign("GDATASETS", character(0), envir = .misha)

    is_per_chromosome <- .is_per_chromosome_db(groot, chromsizes)
    assign("DB_IS_PER_CHROMOSOME", is_per_chromosome, envir = .misha)

    canonical_names <- chromsizes$chrom

    # For per-chromosome databases, add "chr" prefix to match seq file names
    if (is_per_chromosome) {
        # Add chr prefix to names that don't have it
        needs_prefix <- !startsWith(canonical_names, "chr")
        canonical_names[needs_prefix] <- paste0("chr", canonical_names[needs_prefix])
    }

    # Always compute chromosome aliases for better usability.
    # This allows users to use both "chr1" and "1" interchangeably and
    # ensures aliases remain available after database conversion. When a
    # <groot>/chrom_aliases.tsv exists (written by gdb.build_genome for
    # NCBI sources), refseq/GenBank/sequenceName aliases are added too.
    alias_map <- .compute_chrom_aliases(canonical_names, groot = groot)

    # Include both canonical names and aliases in factor levels.
    # `names(alias_map)` always starts with canonical_names by construction
    # (see .compute_chrom_aliases: basic_names <- chroms ...) and is already
    # de-duplicated, so this is equivalent to unique(c(canonical_names,
    # names(alias_map))) without the ~1 s cost of allocating + uniqueing
    # ~5M strings on a million-contig database.
    all_chrom_names <- names(alias_map)

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

    rownames(intervals) <- seq_len(nrow(intervals))

    # Lazy 2D generation: only materialize for small genomes
    contig_threshold <- getOption("gmulticontig.2d.threshold", 100)

    if (nrow(intervals) <= contig_threshold) {
        # Small genome: materialize full 2D grid
        cartesian <- expand.grid(seq_len(nrow(intervals)), seq_len(nrow(intervals)))
        intervals2d <- cbind(intervals[cartesian[, 2], ], intervals[cartesian[, 1], ])
        names(intervals2d) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2")
        # Ensure chrom1 and chrom2 have the same factor levels as intervals$chrom
        # Use canonical_names only (not aliases) for consistency with 1D intervals
        intervals2d$chrom1 <- factor(intervals2d$chrom1, levels = canonical_names)
        intervals2d$chrom2 <- factor(intervals2d$chrom2, levels = canonical_names)
        rownames(intervals2d) <- seq_len(nrow(intervals2d))
    } else {
        # Large genome: defer 2D generation
        intervals2d <- .create_deferred_2d(nrow(intervals))
    }

    assign("ALLGENOME", list(intervals, intervals2d), envir = .misha)
    assign("GROOT", groot, envir = .misha)
    assign("GWD", groot, envir = .misha)
    # Split the alias map: only the basic toggles (chr-prefix + MT) go into
    # CHROM_ALIAS for the C++ side; the env mirror gets everything (including
    # any TSV-derived NCBI/GenBank/sequenceName aliases) for the R-side
    # .gchroms() to consult at query time.
    n_basic <- attr(alias_map, "n_basic")
    if (is.null(n_basic) || n_basic > length(alias_map)) n_basic <- length(alias_map)
    cpp_map <- if (n_basic > 0L) alias_map[seq_len(n_basic)] else stats::setNames(character(0), character(0))
    attr(cpp_map, "n_basic") <- NULL
    assign("CHROM_ALIAS", cpp_map, envir = .misha)
    .refresh_chrom_alias_env(full_alias_map = alias_map)

    # Check if database cache is dirty
    if (.gdb.cache_is_dirty(groot)) {
        rescan <- TRUE
    }

    # Best-effort cleanup of leftover trash dirs from prior sessions. Only the
    # top-level tracks/ directory is swept; per-track parents are handled
    # lazily on next operation in that namespace.
    try(
        .gdb.trash_sweep_old(tracks_dir, max_age_hours = 24),
        silent = TRUE
    )

    success <- FALSE
    tryCatch(
        {
            if (is.null(dir)) {
                # Default to tracks directory
                .gdir.cd(file.path(groot, "tracks"), rescan)
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
                            "The 'dir' parameter ('%s') looks like a misha database path.\nTo connect multiple databases, use: gsetroot(\"%s\"); gdataset.load(\"%s\")",
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
                assign("GWD", NULL, envir = .misha)
                assign("CHROM_ALIAS", NULL, envir = .misha)
                .refresh_chrom_alias_env()
                assign("DB_IS_PER_CHROMOSOME", NULL, envir = .misha)
                assign("GTRACK_DATASET", NULL, envir = .misha)
                assign("GINTERVALS_DATASET", NULL, envir = .misha)
                assign("GDATASETS", character(0), envir = .misha)
            }
        }
    )
}


#' Unloads the genome database
#'
#' Unloads the currently active genome database and clears the session state.
#'
#' Resets the misha session to an uninitialized state: the database root, the
#' working directory, the genome intervals, the virtual tracks, the loaded
#' datasets, the chromosome aliases and the track / intervals-set
#' auto-completion symbols are all cleared. After calling this function
#' \code{\link{gdb.init}} (or \code{\link{gsetroot}}) must be called again
#' before any other misha function is used.
#'
#' Package internals are preserved, so the package itself stays loaded; this is
#' not the same as detaching the package.
#'
#' @return None.
#' @seealso \code{\link{gdb.init}}, \code{\link{gsetroot}},
#' \code{\link{gdb.reload}}
#' @keywords ~database
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gdb.unload()
#' gdb.init_examples()
#'
#' @export gdb.unload
gdb.unload <- function() {
    # Remove the per-name auto-completion symbols first, then their registries.
    for (listvar in c("GTRACKS", "GINTERVS")) {
        if (exists(listvar, envir = .misha, inherits = FALSE)) {
            syms <- get(listvar, envir = .misha)
            syms <- syms[vapply(syms, function(s) exists(s, envir = .misha, inherits = FALSE), logical(1))]
            if (length(syms)) {
                remove(list = syms, envir = .misha)
            }
            remove(list = listvar, envir = .misha)
        }
    }

    session_vars <- c(
        "GROOT", "GWD", "ALLGENOME", "GINTERVID", "GITERATOR.INTERVALS",
        "GVTRACKS", "GDATASETS", "GTRACK_DATASET", "GINTERVALS_DATASET",
        "CHROM_ALIAS", "DB_IS_PER_CHROMOSOME"
    )
    for (v in session_vars) {
        if (exists(v, envir = .misha, inherits = FALSE)) {
            remove(list = v, envir = .misha)
        }
    }

    invisible(NULL)
}

#' Create a linked database with symlinks to a parent database
#'
#' Creates a new database directory structure with symbolic links to the
#' parent database's seq/ directory and chrom_sizes.txt file.
#'
#' This is useful for creating a writable database that shares sequence
#' data with a read-only main database. The new database can be set as the
#' working database via \code{gsetroot()}, and then the parent database
#' can be loaded as a dataset via \code{gdataset.load()}.
#'
#' @param path Path for the new linked database
#' @param parent Path to the parent database (with seq and chrom_sizes.txt)
#' @return Invisible TRUE on success
#' @seealso \code{\link{gsetroot}}, \code{\link{gdb.create}},
#'   \code{\link{gdataset.load}}, \code{\link{gdataset.ls}}
#' @keywords ~db ~database ~create
#' @examples
#' \dontrun{
#' # Create linked database sharing sequence data with main database
#' gdb.create_linked("~/my_tracks", parent = "/shared/genomics/hg38")
#'
#' # Set linked database as working database and load parent as dataset
#' gsetroot("~/my_tracks")
#' gdataset.load("/shared/genomics/hg38")
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
