# Dataset management functions

#' Load a dataset into the namespace
#'
#' Loads tracks and intervals from a dataset directory, making them available
#' for analysis alongside the working database.
#'
#' @param path Path to a dataset or misha database directory
#' @param force If TRUE, ignore name collisions (working db wins; for dataset-to-dataset, later-loaded wins)
#' @param verbose If TRUE, print loaded track/interval names and summary counts
#'
#' @return Invisibly returns a list with:
#'   \item{tracks}{Number of visible tracks loaded}
#'   \item{intervals}{Number of visible intervals loaded}
#'   \item{shadowed_tracks}{Number of tracks shadowed by collisions}
#'   \item{shadowed_intervals}{Number of intervals shadowed by collisions}
#'
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' dataset_path <- gdataset.example_path()
#' gdataset.load(dataset_path)
#' gdataset.unload(dataset_path)
#'
#' @seealso \code{\link{gdataset.unload}}, \code{\link{gdataset.save}}, \code{\link{gdataset.ls}}
#' @export
gdataset.load <- function(path, force = FALSE, verbose = FALSE) {
    # Check working database exists
    if (!exists("GROOT", envir = .misha) || is.null(get("GROOT", envir = .misha))) {
        stop("No working database. Call gsetroot() first.", call. = FALSE)
    }

    # Check path exists before normalizing
    if (!dir.exists(path)) {
        stop(sprintf("Dataset path '%s' does not exist", path), call. = FALSE)
    }

    # Normalize path
    path_norm <- normalizePath(path, mustWork = TRUE)

    groot <- get("GROOT", envir = .misha)

    # If path equals working database, silently return zeros (no-op)
    if (path_norm == groot) {
        return(invisible(list(tracks = 0, intervals = 0, shadowed_tracks = 0, shadowed_intervals = 0)))
    }

    # If already loaded, unload first to support reload
    gdatasets <- get("GDATASETS", envir = .misha)
    if (path_norm %in% gdatasets) {
        gdataset.unload(path_norm)
    }

    # Validate path has tracks/ directory
    if (!dir.exists(file.path(path_norm, "tracks"))) {
        stop(sprintf("Path '%s' does not contain a 'tracks' directory", path), call. = FALSE)
    }

    # Validate chrom_sizes.txt exists and matches working database
    cs_path <- file.path(path_norm, "chrom_sizes.txt")
    if (!file.exists(cs_path)) {
        stop(sprintf("Path '%s' does not contain a chrom_sizes.txt file", path), call. = FALSE)
    }

    # Compare genomes using MD5 hash
    groot_cs <- file.path(groot, "chrom_sizes.txt")
    groot_hash <- tools::md5sum(groot_cs)
    dataset_hash <- tools::md5sum(cs_path)

    if (groot_hash != dataset_hash) {
        stop(sprintf("Cannot load dataset '%s': genome does not match working database", path), call. = FALSE)
    }

    # Scan for tracks and intervals
    dataset_tracks <- .gdb.scan_tracks(path_norm)
    dataset_intervals <- .gdb.scan_intervals(path_norm)

    # Get current state
    .gdb.ensure_dataset_maps()
    gtrack_dataset <- get("GTRACK_DATASET", envir = .misha)
    gintervals_dataset <- get("GINTERVALS_DATASET", envir = .misha)

    # Check for collisions
    track_collisions_working <- intersect(dataset_tracks, names(gtrack_dataset)[gtrack_dataset == groot])
    interval_collisions_working <- intersect(dataset_intervals, names(gintervals_dataset)[gintervals_dataset == groot])

    track_collisions_datasets <- intersect(dataset_tracks, names(gtrack_dataset)[gtrack_dataset != groot])
    interval_collisions_datasets <- intersect(dataset_intervals, names(gintervals_dataset)[gintervals_dataset != groot])

    has_collisions <- length(track_collisions_working) > 0 || length(interval_collisions_working) > 0 ||
        length(track_collisions_datasets) > 0 || length(interval_collisions_datasets) > 0

    if (has_collisions && !force) {
        # Build error message
        msgs <- character()

        if (length(track_collisions_working) > 0) {
            tracks_str <- paste0("'", track_collisions_working, "'", collapse = ", ")
            msgs <- c(msgs, sprintf("tracks %s already exist in working database '%s'", tracks_str, groot))
        }

        if (length(interval_collisions_working) > 0) {
            intervals_str <- paste0("'", interval_collisions_working, "'", collapse = ", ")
            msgs <- c(msgs, sprintf("interval sets %s already exist in working database '%s'", intervals_str, groot))
        }

        if (length(track_collisions_datasets) > 0) {
            # Find which dataset(s) have these tracks
            collision_dbs <- unique(gtrack_dataset[track_collisions_datasets])
            tracks_str <- paste0("'", track_collisions_datasets, "'", collapse = ", ")
            db_str <- paste0("'", collision_dbs, "'", collapse = ", ")
            msgs <- c(msgs, sprintf("tracks %s already exist in loaded dataset %s", tracks_str, db_str))
        }

        if (length(interval_collisions_datasets) > 0) {
            collision_dbs <- unique(gintervals_dataset[interval_collisions_datasets])
            intervals_str <- paste0("'", interval_collisions_datasets, "'", collapse = ", ")
            db_str <- paste0("'", collision_dbs, "'", collapse = ", ")
            msgs <- c(msgs, sprintf("interval sets %s already exist in loaded dataset %s", intervals_str, db_str))
        }

        error_msg <- sprintf("Cannot load dataset '%s':\n  - %s\nUse force=TRUE to override.", path, paste(msgs, collapse = "\n  - "))
        stop(error_msg, call. = FALSE)
    }

    # Load tracks and intervals
    shadowed_tracks <- 0
    shadowed_intervals <- 0
    visible_tracks <- 0
    visible_intervals <- 0

    for (track in dataset_tracks) {
        if (track %in% names(gtrack_dataset)) {
            # Collision - skip if working db, otherwise shadow
            if (gtrack_dataset[[track]] == groot) {
                # Working db wins - skip dataset track
                shadowed_tracks <- shadowed_tracks + 1
            } else {
                # Dataset collision - new dataset wins
                gtrack_dataset[[track]] <- path_norm
                visible_tracks <- visible_tracks + 1
            }
        } else {
            # No collision - add track
            gtrack_dataset[[track]] <- path_norm
            visible_tracks <- visible_tracks + 1
        }
    }

    for (interval in dataset_intervals) {
        if (interval %in% names(gintervals_dataset)) {
            if (gintervals_dataset[[interval]] == groot) {
                shadowed_intervals <- shadowed_intervals + 1
            } else {
                gintervals_dataset[[interval]] <- path_norm
                visible_intervals <- visible_intervals + 1
            }
        } else {
            gintervals_dataset[[interval]] <- path_norm
            visible_intervals <- visible_intervals + 1
        }
    }

    # Add to loaded datasets first (needed for gdb.reload)
    gdatasets <- c(gdatasets, path_norm)
    assign("GDATASETS", gdatasets, envir = .misha)

    # Reload the database to properly scan all tracks/intervals
    gdb.reload(rescan = TRUE)

    if (verbose) {
        message(sprintf("Loaded dataset '%s':", path))
        message(sprintf("  Tracks: %d visible, %d shadowed", visible_tracks, shadowed_tracks))
        message(sprintf("  Intervals: %d visible, %d shadowed", visible_intervals, shadowed_intervals))
    }

    invisible(list(
        tracks = visible_tracks,
        intervals = visible_intervals,
        shadowed_tracks = shadowed_tracks,
        shadowed_intervals = shadowed_intervals
    ))
}

#' Unload a dataset from the namespace
#'
#' Removes all tracks and intervals from a previously loaded dataset.
#' If a track was shadowing another, the shadowed track becomes visible again.
#'
#' @param path Path to a previously loaded dataset
#' @param validate If TRUE, error if path is not currently loaded; otherwise silently no-op
#'
#' @return Invisible NULL
#'
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' dataset_path <- gdataset.example_path()
#' gdataset.load(dataset_path)
#' gdataset.unload(dataset_path, validate = TRUE)
#'
#' @seealso \code{\link{gdataset.load}}, \code{\link{gdataset.ls}}
#' @export
gdataset.unload <- function(path, validate = FALSE) {
    # Normalize path
    path_norm <- normalizePath(path, mustWork = FALSE)

    gdatasets <- get("GDATASETS", envir = .misha)

    # Check if loaded
    if (!path_norm %in% gdatasets) {
        if (validate) {
            stop(sprintf("Dataset '%s' is not loaded", path), call. = FALSE)
        }
        return(invisible(NULL))
    }

    # Remove from loaded datasets
    gdatasets <- setdiff(gdatasets, path_norm)
    assign("GDATASETS", gdatasets, envir = .misha)

    # Reload the database to properly rescan all tracks/intervals
    gdb.reload(rescan = TRUE)

    invisible(NULL)
}

#' Save a dataset
#'
#' Creates a new dataset directory containing selected tracks and/or intervals
#' from the working database.
#'
#' @param path Destination directory (must not exist)
#' @param description Required description for metadata
#' @param tracks Character vector of track names to include
#' @param intervals Character vector of interval set names to include
#' @param symlinks If TRUE, create symlinks to tracks/intervals instead of copying
#' @param copy_seq If TRUE, copy seq/ directory instead of symlinking
#'
#' @return Invisible path
#'
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' example_intervs <- gintervals(1, 0, 10000)
#' gintervals.save("example_dataset_intervals", example_intervs)
#' gtrack.create(
#'     "example_dataset_track",
#'     "Example dataset track",
#'     "dense_track",
#'     iterator = "example_dataset_intervals"
#' )
#' dataset_path <- tempfile("misha_dataset_")
#' gdataset.save(
#'     path = dataset_path,
#'     description = "Example dataset",
#'     tracks = "example_dataset_track",
#'     intervals = "example_dataset_intervals"
#' )
#' gtrack.rm("example_dataset_track", force = TRUE)
#' gintervals.rm("example_dataset_intervals", force = TRUE)
#'
#' @seealso \code{\link{gdataset.load}}, \code{\link{gdataset.info}}
#' @export
gdataset.save <- function(path, description, tracks = NULL, intervals = NULL,
                          symlinks = FALSE, copy_seq = FALSE) {
    .gcheckroot()
    .gdb.ensure_dataset_maps()

    # Validate at least one of tracks or intervals specified
    if (is.null(tracks) && is.null(intervals)) {
        stop("At least one of 'tracks' or 'intervals' must be specified", call. = FALSE)
    }

    # Validate path doesn't exist
    if (file.exists(path) || dir.exists(path)) {
        stop(sprintf("Path '%s' already exists", path), call. = FALSE)
    }

    groot <- get("GROOT", envir = .misha)
    gtrack_dataset <- get("GTRACK_DATASET", envir = .misha)
    gintervals_dataset <- get("GINTERVALS_DATASET", envir = .misha)

    # Validate tracks exist
    if (!is.null(tracks)) {
        for (track in tracks) {
            if (!track %in% names(gtrack_dataset)) {
                stop(sprintf("Track '%s' does not exist", track), call. = FALSE)
            }
        }
    }

    # Validate intervals exist
    if (!is.null(intervals)) {
        for (interval in intervals) {
            if (!interval %in% names(gintervals_dataset)) {
                stop(sprintf("Interval set '%s' does not exist", interval), call. = FALSE)
            }
        }
    }

    # Create directory structure
    dir.create(path, recursive = TRUE)
    dir.create(file.path(path, "tracks"), recursive = TRUE)

    # Copy chrom_sizes.txt
    file.copy(file.path(groot, "chrom_sizes.txt"), file.path(path, "chrom_sizes.txt"))

    # Handle seq/ directory
    if (copy_seq) {
        # Copy the seq directory contents directly into path/seq
        # Note: file.copy with recursive=TRUE copies the contents into an existing directory
        seq_dest <- file.path(path, "seq")
        dir.create(seq_dest, showWarnings = FALSE)
        seq_files <- list.files(file.path(groot, "seq"), full.names = TRUE)
        file.copy(seq_files, seq_dest, recursive = TRUE)
    } else {
        file.symlink(file.path(groot, "seq"), file.path(path, "seq"))
    }

    # Copy/link tracks
    if (!is.null(tracks)) {
        for (track in tracks) {
            source_db <- gtrack_dataset[[track]]
            source_path <- file.path(source_db, "tracks", paste0(track, ".track"))
            dest_path <- file.path(path, "tracks", paste0(track, ".track"))

            if (symlinks) {
                file.symlink(source_path, dest_path)
            } else {
                file.copy(source_path, file.path(path, "tracks"), recursive = TRUE)
            }
        }
    }

    # Copy/link intervals
    if (!is.null(intervals)) {
        for (interval in intervals) {
            source_db <- gintervals_dataset[[interval]]
            # Intervals can be files or directories
            source_file <- file.path(source_db, "tracks", paste0(interval, ".interv"))
            source_dir <- source_file # Same path - could be file or dir

            if (file.exists(source_file) || dir.exists(source_dir)) {
                dest_path <- file.path(path, "tracks", paste0(interval, ".interv"))

                if (symlinks) {
                    file.symlink(source_file, dest_path)
                } else {
                    file.copy(source_file, file.path(path, "tracks"), recursive = TRUE)
                }
            }
        }
    }

    # Generate misha.yaml
    yaml_data <- list(
        description = description,
        author = Sys.info()["user"],
        created = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ"),
        original_db = groot,
        misha_version = as.character(utils::packageVersion("misha")),
        track_count = length(tracks),
        interval_count = length(intervals),
        genome = as.character(tools::md5sum(file.path(groot, "chrom_sizes.txt")))
    )

    yaml::write_yaml(yaml_data, file.path(path, "misha.yaml"))

    invisible(path)
}

#' List working database and loaded datasets
#'
#' Returns a list of the working database and all loaded datasets.
#'
#' @param dataframe If FALSE, return character vector; if TRUE, return data frame
#'
#' @return Character vector of paths or data frame with detailed information
#'
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' dataset_path <- gdataset.example_path()
#' gdataset.load(dataset_path)
#' gdataset.ls()
#' gdataset.unload(dataset_path)
#'
#' @seealso \code{\link{gdataset.load}}, \code{\link{gdataset.info}}
#' @export
gdataset.ls <- function(dataframe = FALSE) {
    .gcheckroot()

    groot <- get("GROOT", envir = .misha)
    gdatasets <- get("GDATASETS", envir = .misha)

    all_paths <- c(groot, gdatasets)

    if (!dataframe) {
        return(all_paths)
    }

    # Build data frame
    gtrack_dataset <- get("GTRACK_DATASET", envir = .misha)
    gintervals_dataset <- get("GINTERVALS_DATASET", envir = .misha)

    result <- data.frame(
        path = all_paths,
        tracks_total = integer(length(all_paths)),
        tracks_visible = integer(length(all_paths)),
        intervals_total = integer(length(all_paths)),
        intervals_visible = integer(length(all_paths)),
        has_metadata = logical(length(all_paths)),
        writable = logical(length(all_paths)),
        stringsAsFactors = FALSE
    )

    for (i in seq_along(all_paths)) {
        db <- all_paths[i]

        # Scan filesystem
        result$tracks_total[i] <- length(.gdb.scan_tracks(db))
        result$intervals_total[i] <- length(.gdb.scan_intervals(db))

        # Count visible
        if (is.null(gtrack_dataset)) {
            result$tracks_visible[i] <- if (db == groot) length(get("GTRACKS", envir = .misha)) else 0L
        } else {
            result$tracks_visible[i] <- sum(gtrack_dataset == db, na.rm = TRUE)
        }

        if (is.null(gintervals_dataset)) {
            result$intervals_visible[i] <- if (db == groot) length(get("GINTERVS", envir = .misha)) else 0L
        } else {
            result$intervals_visible[i] <- sum(gintervals_dataset == db, na.rm = TRUE)
        }

        result$has_metadata[i] <- file.exists(file.path(db, "misha.yaml"))
        result$writable[i] <- (db == groot)
    }

    result
}

#' Get dataset information
#'
#' Returns metadata and contents of a dataset.
#'
#' @param path Path to any dataset (loaded or not)
#'
#' @return List with dataset information
#'
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' dataset_path <- gdataset.example_path()
#' gdataset.info(dataset_path)
#'
#' @seealso \code{\link{gdataset.ls}}, \code{\link{gdataset.load}}
#' @export
gdataset.info <- function(path) {
    path_norm <- normalizePath(path, mustWork = TRUE)

    # Read misha.yaml if exists
    yaml_path <- file.path(path_norm, "misha.yaml")
    if (file.exists(yaml_path)) {
        yaml_data <- yaml::read_yaml(yaml_path)
    } else {
        yaml_data <- list()
    }

    # Scan tracks and intervals
    tracks <- .gdb.scan_tracks(path_norm)
    intervals <- .gdb.scan_intervals(path_norm)

    # Check if loaded
    gdatasets <- get("GDATASETS", envir = .misha)
    is_loaded <- path_norm %in% gdatasets

    # Compute genome hash
    cs_path <- file.path(path_norm, "chrom_sizes.txt")
    genome_hash <- if (file.exists(cs_path)) as.character(tools::md5sum(cs_path)) else NA

    list(
        description = yaml_data$description,
        author = yaml_data$author,
        created = yaml_data$created,
        original_db = yaml_data$original_db,
        misha_version = yaml_data$misha_version,
        track_count = length(tracks),
        interval_count = length(intervals),
        genome = genome_hash,
        is_loaded = is_loaded
    )
}

#' Create an example dataset on the fly
#'
#' Creates a small dataset in a temporary directory using the built-in
#' example database. The working database is reset via
#' \code{\link{gdb.init_examples}}.
#' Temporary tracks and intervals are removed after saving to avoid
#' name collisions.
#'
#' @return Path to the created dataset directory
#'
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' dataset_path <- gdataset.example_path()
#' gdataset.load(dataset_path)
#' gdataset.unload(dataset_path)
#'
#' @seealso \code{\link{gdataset.save}}, \code{\link{gdataset.load}},
#'   \code{\link{gdb.init_examples}}
#' @export
gdataset.example_path <- function() {
    gdb.init_examples()

    gtrack.rm("example_dataset_track", force = TRUE)
    gintervals.rm("example_dataset_intervals", force = TRUE)

    example_intervs <- gintervals(1, 0, 10000)
    gintervals.save("example_dataset_intervals", example_intervs)

    gtrack.create(
        "example_dataset_track",
        "Example dataset track",
        "dense_track",
        iterator = "example_dataset_intervals"
    )

    dataset_path <- tempfile("misha_dataset_")
    gdataset.save(
        path = dataset_path,
        description = "Example dataset created on the fly",
        tracks = "example_dataset_track",
        intervals = "example_dataset_intervals"
    )

    gtrack.rm("example_dataset_track", force = TRUE)
    gintervals.rm("example_dataset_intervals", force = TRUE)

    dataset_path
}

# Helper functions

.gdb.scan_tracks <- function(db) {
    tracks_dir <- file.path(db, "tracks")
    if (!dir.exists(tracks_dir)) {
        return(character(0))
    }

    # List all .track directories recursively
    all_files <- list.files(tracks_dir, full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
    track_dirs <- all_files[grepl("\\.track$", all_files) & dir.exists(all_files)]

    if (length(track_dirs) == 0) {
        return(character(0))
    }

    # Extract track names (relative to tracks_dir, without .track extension)
    track_names <- sub("\\.track$", "", sub(paste0("^", tracks_dir, "/"), "", track_dirs))

    track_names
}

.gdb.scan_intervals <- function(db) {
    tracks_dir <- file.path(db, "tracks")
    if (!dir.exists(tracks_dir)) {
        return(character(0))
    }

    # List all .interv files/directories
    all_files <- list.files(tracks_dir, full.names = TRUE, recursive = TRUE)
    interv_paths <- all_files[grepl("\\.interv$", all_files)]

    if (length(interv_paths) == 0) {
        return(character(0))
    }

    # Extract interval names
    interval_names <- sub("\\.interv$", "", sub(paste0("^", tracks_dir, "/"), "", interv_paths))

    # Remove duplicates (in case both file and directory exist)
    unique(interval_names)
}
