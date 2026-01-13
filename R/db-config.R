# Database configuration and prefix management functions

# Read .misha config file from a database
#
# @param db_path Path to the database root directory
# @return List with config fields (prefix, description, author, version, created_at).
#         Returns list with NULL values if no .misha file or fields not present.
.gdb.read_config <- function(db_path) {
    config_path <- file.path(db_path, ".misha")
    default_config <- list(
        prefix = NULL,
        global = NULL,
        description = NULL,
        author = NULL,
        version = NULL,
        created_at = NULL
    )

    if (!file.exists(config_path)) {
        return(default_config)
    }

    config <- tryCatch(
        {
            yaml::read_yaml(config_path)
        },
        error = function(e) {
            warning(sprintf("Failed to read .misha config from '%s': %s", db_path, e$message), call. = FALSE)
            return(NULL)
        }
    )

    if (is.null(config)) {
        return(default_config)
    }

    # Merge with defaults
    for (field in names(default_config)) {
        if (!is.null(config[[field]])) {
            default_config[[field]] <- config[[field]]
        }
    }

    # Validate prefix if present
    if (!is.null(default_config$prefix)) {
        .gdb.validate_prefix(default_config$prefix, db_path)
    }

    default_config
}

# Validate a database prefix
#
# @param prefix The prefix string to validate
# @param db_path Path to the database (for error messages)
# @return TRUE if valid, errors otherwise
.gdb.validate_prefix <- function(prefix, db_path = NULL) {
    db_info <- if (!is.null(db_path)) sprintf(" in database '%s'", db_path) else ""

    # Must be a single string
    if (!is.character(prefix) || length(prefix) != 1 || is.na(prefix)) {
        stop(sprintf("Invalid prefix%s: prefix must be a single string", db_info), call. = FALSE)
    }

    # Must not be empty
    if (nchar(prefix) == 0) {
        stop(sprintf("Invalid prefix%s: prefix cannot be empty", db_info), call. = FALSE)
    }

    # Must match pattern: start with letter, alphanumeric + underscore
    if (!grepl("^[a-zA-Z][a-zA-Z0-9_]*$", prefix)) {
        stop(sprintf(
            "Invalid prefix '%s'%s: prefix must start with a letter and contain only alphanumeric characters and underscores",
            prefix, db_info
        ), call. = FALSE)
    }

    # Must not contain @
    if (grepl("@", prefix, fixed = TRUE)) {
        stop(sprintf("Invalid prefix '%s'%s: prefix cannot contain '@' character", prefix, db_info), call. = FALSE)
    }

    # Warn if longer than recommended
    if (nchar(prefix) > 8) {
        warning(sprintf(
            "Prefix '%s'%s is longer than recommended (2-8 characters)",
            prefix, db_info
        ), call. = FALSE)
    }

    invisible(TRUE)
}

# Build prefix map from connected databases
#
# @param groots Vector of database paths
# @param configs Named list of database configs (db_path -> config)
# @return Named list mapping prefix -> db_path
.gdb.build_prefix_map <- function(groots, configs) {
    prefix_map <- list()

    for (groot in groots) {
        config <- configs[[groot]]
        if (!is.null(config) && !is.null(config$prefix)) {
            prefix <- config$prefix

            # Check for duplicate prefixes
            if (prefix %in% names(prefix_map)) {
                stop(sprintf(
                    "Duplicate prefix '%s' in databases '%s' and '%s'",
                    prefix, prefix_map[[prefix]], groot
                ), call. = FALSE)
            }

            prefix_map[[prefix]] <- groot
        }
    }

    prefix_map
}


#' Get database configuration
#'
#' Returns the configuration from a database's `.misha` file.
#'
#' Database configuration is read from a YAML file named `.misha` in the
#' database root directory. This file can contain a `prefix` for namespacing
#' tracks, along with optional metadata like `description`, `author`, etc.
#'
#' @param db Database identifier: NULL for current GWD database, a path to a
#'   database directory, or a prefix string (e.g., "at") for a prefixed database.
#' @return A list with the following fields (NULL if not specified):
#' \describe{
#'   \item{prefix}{Character: Database prefix for track namespacing}
#'   \item{global}{Logical: TRUE if explicitly marked as global (no prefix)}
#'   \item{description}{Character: Human-readable database description}
#'   \item{author}{Character: Database author or maintainer}
#'
#'   \item{version}{Character: Database version string}
#'   \item{created_at}{Character: Creation date (ISO 8601 format)}
#' }
#' @seealso \code{\link{gdb.prefixes}}, \code{\link{gsetroot}}, \code{\link{gdb.ls}}
#' @keywords ~db ~database ~config ~prefix
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gdb.config()
#'
#' @export gdb.config
gdb.config <- function(db = NULL) {
    .gcheckroot()

    db_path <- .gresolve_db_arg(db)
    if (is.null(db_path)) {
        stop("Database not found", call. = FALSE)
    }

    configs <- get("GDB_CONFIGS", envir = .misha)
    if (!is.null(configs) && db_path %in% names(configs)) {
        return(configs[[db_path]])
    }

    # Read config directly if not cached
    .gdb.read_config(db_path)
}


#' List database prefixes
#'
#' Returns a mapping of database prefixes to their paths.
#'
#' When databases are connected with `.misha` configuration files that specify
#' a `prefix` field, this function returns a named vector mapping each prefix
#' to its database path. Global databases (without prefixes) are not included.
#'
#' @return Named character vector where names are prefixes and values are
#'   database paths. Returns an empty named character vector if no prefixed
#'   databases are connected.
#' @seealso \code{\link{gdb.config}}, \code{\link{gsetroot}}, \code{\link{gdb.ls}}
#' @keywords ~db ~database ~prefix
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gdb.prefixes()
#'
#' @export gdb.prefixes
gdb.prefixes <- function() {
    .gcheckroot()

    prefix_map <- get("GPREFIX_MAP", envir = .misha)
    if (is.null(prefix_map) || length(prefix_map) == 0) {
        return(structure(character(0), names = character(0)))
    }

    # Convert list to named vector
    result <- unlist(prefix_map, use.names = TRUE)
    result
}


# Resolve a db argument to a database path
#
# @param db NULL (current GWD), a path, or a prefix string
# @return Database path or NULL if not found
.gresolve_db_arg <- function(db) {
    if (is.null(db)) {
        # Use current GWD's database
        gwd <- get("GWD", envir = .misha)
        return(.gdb.resolve_db_for_path(gwd))
    }

    # Try as prefix first
    prefix_map <- get("GPREFIX_MAP", envir = .misha)
    if (!is.null(prefix_map) && db %in% names(prefix_map)) {
        return(prefix_map[[db]])
    }

    # Try as path
    groots <- get("GROOTS", envir = .misha)
    if (is.null(groots)) {
        return(NULL)
    }

    # Normalize and check against connected databases
    db_normalized <- tryCatch(
        normalizePath(db, mustWork = FALSE),
        error = function(e) db
    )

    for (groot in groots) {
        if (identical(groot, db_normalized) || identical(groot, db)) {
            return(groot)
        }
    }

    NULL
}
# Note: .gdb.resolve_db_for_path is defined in db-management.R
