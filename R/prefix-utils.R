# Prefix parsing and resolution utilities for multi-database support

# The prefix separator character
.GPREFIX_SEP <- "@"

# Parse a prefixed name into prefix and base name components
#
# @param name A track or interval name, possibly with prefix (e.g., "at@track.name")
# @return List with prefix (NULL if no prefix) and name (the base name)
# @examples
# .gparse_prefixed_name("at@track.name")  # list(prefix = "at", name = "track.name")
# .gparse_prefixed_name("track.name")     # list(prefix = NULL, name = "track.name")
.gparse_prefixed_name <- function(name) {
    if (is.null(name) || length(name) == 0 || is.na(name)) {
        return(list(prefix = NULL, name = name))
    }

    # If ignoring prefixes, treat @ as literal character
    if (isTRUE(getOption("misha.dangerously_ignore_db_prefixes", FALSE))) {
        return(list(prefix = NULL, name = name))
    }

    # Find first @ character
    at_pos <- regexpr(.GPREFIX_SEP, name, fixed = TRUE)

    if (at_pos == -1) {
        # No prefix
        return(list(prefix = NULL, name = name))
    }

    prefix <- substr(name, 1, at_pos - 1)
    base_name <- substr(name, at_pos + 1, nchar(name))

    # Validate prefix is not empty and base name is not empty
    if (nchar(prefix) == 0 || nchar(base_name) == 0) {
        stop(sprintf(
            "Invalid prefixed name '%s': both prefix and name must be non-empty",
            name
        ), call. = FALSE)
    }

    list(prefix = prefix, name = base_name)
}


# Check if a name has a prefix
#
# @param name A track or interval name
# @return TRUE if the name contains @, FALSE otherwise
.ghas_prefix <- function(name) {
    if (isTRUE(getOption("misha.dangerously_ignore_db_prefixes", FALSE))) {
        return(FALSE)
    }
    !is.null(name) && length(name) > 0 && !is.na(name) &&
        grepl(.GPREFIX_SEP, name, fixed = TRUE)
}


# Look up database path for a prefix
#
# @param prefix A database prefix string
# @return Database path or NULL if prefix not found
.gprefix_to_db <- function(prefix) {
    if (is.null(prefix)) {
        return(NULL)
    }

    prefix_map <- get("GPREFIX_MAP", envir = .misha)
    if (is.null(prefix_map) || !(prefix %in% names(prefix_map))) {
        return(NULL)
    }

    prefix_map[[prefix]]
}


# Get the prefix for a database path
#
# @param db_path A database path
# @return Prefix string or NULL if database has no prefix
.gdb_to_prefix <- function(db_path) {
    if (is.null(db_path)) {
        return(NULL)
    }

    configs <- get("GDB_CONFIGS", envir = .misha)
    if (is.null(configs) || !(db_path %in% names(configs))) {
        return(NULL)
    }

    config <- configs[[db_path]]
    if (is.null(config)) {
        return(NULL)
    }

    config$prefix
}


# Check if a database is global (no prefix)
#
# @param db_path A database path
# @return TRUE if database has no prefix, FALSE if it has a prefix
.gis_global_db <- function(db_path) {
    is.null(.gdb_to_prefix(db_path))
}


# Qualify a name with a database prefix
#
# @param name Base track/interval name (without prefix)
# @param db_path Database path
# @return Prefixed name if database has prefix, otherwise name as-is
.gqualify_name <- function(name, db_path) {
    prefix <- .gdb_to_prefix(db_path)
    if (is.null(prefix)) {
        return(name)
    }
    paste0(prefix, .GPREFIX_SEP, name)
}


# Strip prefix from a name
#
# @param name A track or interval name, possibly with prefix
# @return The base name without prefix
.gunqualify_name <- function(name) {
    parsed <- .gparse_prefixed_name(name)
    parsed$name
}


# Stop with error for unknown database prefix
#
# @param prefix The unknown prefix string
.gstop_unknown_prefix <- function(prefix) {
    stop(sprintf("Unknown database prefix '%s'", prefix), call. = FALSE)
}


# Internal helper to resolve a name to its database path
#
# @param name A track or interval name (possibly prefixed)
# @param db_map_name Name of the database map variable ("GTRACK_DB" or "GINTERVALS_DB")
# @return Database path or NULL if not found
.gresolve_name_db <- function(name, db_map_name) {
    parsed <- .gparse_prefixed_name(name)

    if (!is.null(parsed$prefix)) {
        db_path <- .gprefix_to_db(parsed$prefix)
        if (is.null(db_path)) {
            .gstop_unknown_prefix(parsed$prefix)
        }
        return(db_path)
    }

    db_map <- get(db_map_name, envir = .misha)
    if (is.null(db_map) || !(name %in% names(db_map))) {
        return(NULL)
    }

    db_map[[name]]
}


# Resolve a track name to its database path
#
# @param trackname A track name (possibly prefixed)
# @return Database path or NULL if not found
.gresolve_track_db <- function(trackname) {
    .gresolve_name_db(trackname, "GTRACK_DB")
}


# Resolve an interval name to its database path
#
# @param intervname An interval name (possibly prefixed)
# @return Database path or NULL if not found
.gresolve_intervals_db <- function(intervname) {
    .gresolve_name_db(intervname, "GINTERVALS_DB")
}


# Internal helper to suggest prefixed names when an unprefixed name is not found
#
# @param name An unprefixed track or interval name that was not found
# @param extension File extension (".track" or ".interv")
# @param check_fn Function to check existence (dir.exists for tracks, file.exists for intervals)
# @return Suggestion string or NULL if not found in any prefixed database
.gsuggest_prefixed_name <- function(name, extension, check_fn) {
    if (.ghas_prefix(name)) {
        return(NULL)
    }

    prefix_map <- get("GPREFIX_MAP", envir = .misha)
    if (is.null(prefix_map) || length(prefix_map) == 0) {
        return(NULL)
    }

    found_in <- character(0)
    for (prefix in names(prefix_map)) {
        db_path <- prefix_map[[prefix]]
        item_path <- file.path(db_path, "tracks", paste0(gsub("\\.", "/", name), extension))
        if (check_fn(item_path)) {
            found_in <- c(found_in, paste0(prefix, .GPREFIX_SEP, name))
        }
    }

    if (length(found_in) == 0) {
        return(NULL)
    }

    if (length(found_in) == 1) {
        sprintf("Did you mean: %s?", found_in[1])
    } else {
        sprintf("Did you mean: %s?", paste(found_in, collapse = " or "))
    }
}


# Suggest prefixed track names when an unprefixed track is not found
#
# @param trackname An unprefixed track name that was not found
# @return Suggestion string or NULL if not found in any prefixed database
.gsuggest_prefixed_track <- function(trackname) {
    .gsuggest_prefixed_name(trackname, ".track", dir.exists)
}


# Suggest prefixed interval names when an unprefixed interval is not found
#
# @param intervname An unprefixed interval name that was not found
# @return Suggestion string or NULL if not found in any prefixed database
.gsuggest_prefixed_intervals <- function(intervname) {
    .gsuggest_prefixed_name(intervname, ".interv", file.exists)
}


# Resolve target database for track creation from a prefixed name
#
# @param trackname Track name (possibly prefixed)
# @return List with db_path (target database) and base_name (name without prefix)
.gresolve_creation_target <- function(trackname) {
    parsed <- .gparse_prefixed_name(trackname)

    if (!is.null(parsed$prefix)) {
        db_path <- .gprefix_to_db(parsed$prefix)
        if (is.null(db_path)) {
            .gstop_unknown_prefix(parsed$prefix)
        }
        return(list(db_path = db_path, base_name = parsed$name))
    }

    gwd <- get("GWD", envir = .misha)
    db_path <- .gdb.resolve_db_for_path(gwd)

    list(db_path = db_path, base_name = trackname)
}


# Execute a function with GWD temporarily set to a specific database
#
# @param db_path Target database path
# @param fn Function to execute
# @return Result of fn()
.gwith_db_context <- function(db_path, fn) {
    if (is.null(db_path)) {
        return(fn())
    }

    target_gwd <- file.path(db_path, "tracks")
    current_gwd <- get("GWD", envir = .misha)

    if (identical(target_gwd, current_gwd)) {
        return(fn())
    }

    assign("GWD", target_gwd, envir = .misha)
    on.exit(assign("GWD", current_gwd, envir = .misha), add = TRUE)

    fn()
}


# Stop with an error message about a non-existent track, including suggestions
#
# @param trackname The track name that was not found
# @param msg Optional custom message prefix
.gstop_track_not_found <- function(trackname, msg = NULL) {
    if (is.null(msg)) {
        msg <- sprintf("Track %s does not exist", trackname)
    }

    # Try to suggest prefixed alternatives
    suggestion <- .gsuggest_prefixed_track(trackname)
    if (!is.null(suggestion)) {
        msg <- paste0(msg, ". ", suggestion)
    }

    stop(msg, call. = FALSE)
}


# Stop with an error message about a non-existent intervals set, including suggestions
#
# @param intervname The intervals set name that was not found
# @param msg Optional custom message prefix
.gstop_intervals_not_found <- function(intervname, msg = NULL) {
    if (is.null(msg)) {
        msg <- sprintf("Intervals set %s does not exist", intervname)
    }

    # Try to suggest prefixed alternatives
    suggestion <- .gsuggest_prefixed_intervals(intervname)
    if (!is.null(suggestion)) {
        msg <- paste0(msg, ". ", suggestion)
    }

    stop(msg, call. = FALSE)
}
