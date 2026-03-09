# Interval set attribute management

# ---- Internal helpers ----

#' Resolve path to .iattr file for an interval set
#' @noRd
.gintervals.attr_path <- function(intervals.set) {
    # Look up which database this interval set belongs to
    intervals_db <- get("GINTERVALS_DATASET", envir = .misha)
    if (!is.null(intervals_db) && intervals.set %in% names(intervals_db)) {
        db_root <- intervals_db[intervals.set]
    } else {
        db_root <- get("GROOT", envir = .misha)
    }

    path <- gsub("\\.", "/", intervals.set)
    interv_path <- file.path(db_root, "tracks", paste0(path, ".interv"))

    if (dir.exists(interv_path)) {
        # Big set: store inside the directory
        file.path(interv_path, ".iattr")
    } else {
        # Small set: replace .interv with .iattr
        file.path(db_root, "tracks", paste0(path, ".iattr"))
    }
}

#' Read binary .iattr file into a named character vector
#' @noRd
.gintervals.attr_read <- function(path) {
    if (!file.exists(path)) {
        return(character(0))
    }
    size <- file.info(path)$size
    if (is.na(size) || size == 0) {
        return(character(0))
    }
    raw <- readBin(path, "raw", size)
    if (length(raw) == 0) {
        return(character(0))
    }

    # Find null byte positions
    null_pos <- which(raw == as.raw(0L))
    if (length(null_pos) == 0) {
        return(character(0))
    }

    # Extract strings between null bytes
    starts <- c(1L, null_pos + 1L)
    ends <- null_pos - 1L

    # Only keep complete pairs
    n_strings <- length(null_pos)
    if (n_strings %% 2 != 0) {
        # Odd number of null-terminated strings; drop the last incomplete pair
        n_strings <- n_strings - 1L
    }
    if (n_strings == 0) {
        return(character(0))
    }

    strings <- vapply(seq_len(n_strings), function(i) {
        if (starts[i] > ends[i]) {
            ""
        } else {
            rawToChar(raw[starts[i]:ends[i]])
        }
    }, character(1))

    # Pair up as key/value
    keys <- strings[seq(1L, n_strings, by = 2L)]
    values <- strings[seq(2L, n_strings, by = 2L)]
    stats::setNames(values, keys)
}

#' Write named character vector to binary .iattr file
#' @noRd
.gintervals.attr_write <- function(path, attrs) {
    if (length(attrs) == 0) {
        if (file.exists(path)) {
            unlink(path)
        }
        return(invisible())
    }

    old_umask <- Sys.umask("002")
    on.exit(Sys.umask(old_umask), add = TRUE)

    con <- file(path, "wb")
    on.exit(close(con), add = TRUE)
    for (nm in names(attrs)) {
        writeBin(charToRaw(nm), con)
        writeBin(as.raw(0L), con)
        writeBin(charToRaw(as.character(attrs[nm])), con)
        writeBin(as.raw(0L), con)
    }
    invisible()
}

#' Check that an interval set belongs to the working database (not read-only)
#' @noRd
.gintervals.attr_check_writable <- function(intervals.set) {
    intervals_db <- get("GINTERVALS_DATASET", envir = .misha)
    groot <- get("GROOT", envir = .misha)
    if (!is.null(intervals_db) && intervals.set %in% names(intervals_db)) {
        db <- intervals_db[intervals.set]
        if (!identical(unname(db), groot)) {
            stop(sprintf(
                "Intervals set %s belongs to a read-only dataset and cannot be modified",
                intervals.set
            ), call. = FALSE)
        }
    }
}


# ---- Public functions ----

#' Returns value of an interval set attribute
#'
#' Returns value of an interval set attribute.
#'
#' This function returns the value of an interval set attribute. If the
#' attribute does not exist an empty string is returned.
#'
#' @param intervals.set interval set name
#' @param attr attribute name
#' @return Interval set attribute value (character string).
#' @seealso \code{\link{gintervals.attr.set}},
#' \code{\link{gintervals.attr.export}}, \code{\link{gintervals.attr.import}}
#' @keywords ~attr ~attribute ~intervals
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gintervals.attr.set("annotations1", "test_attr", "value")
#' gintervals.attr.get("annotations1", "test_attr")
#' gintervals.attr.set("annotations1", "test_attr", "")
#'
#' @export gintervals.attr.get
gintervals.attr.get <- function(intervals.set = NULL, attr = NULL) {
    if (is.null(intervals.set) || is.null(attr)) {
        stop("Usage: gintervals.attr.get(intervals.set, attr)", call. = FALSE)
    }
    .gcheckroot()

    res <- gintervals.attr.export(intervals.set, attr)
    res[1, 1]
}


#' Assigns value to an interval set attribute
#'
#' Assigns value to an interval set attribute.
#'
#' This function creates an interval set attribute and assigns 'value' to it.
#' If the attribute already exists its value is overwritten.
#'
#' If 'value' is an empty string the attribute is removed.
#'
#' @param intervals.set interval set name
#' @param attr attribute name
#' @param value value (character string)
#' @return None.
#' @seealso \code{\link{gintervals.attr.get}},
#' \code{\link{gintervals.attr.export}}, \code{\link{gintervals.attr.import}}
#' @keywords ~attr ~attribute ~intervals
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gintervals.attr.set("annotations1", "test_attr", "value")
#' gintervals.attr.get("annotations1", "test_attr")
#' gintervals.attr.set("annotations1", "test_attr", "")
#'
#' @export gintervals.attr.set
gintervals.attr.set <- function(intervals.set = NULL, attr = NULL, value = NULL) {
    if (is.null(intervals.set) || is.null(attr) || is.null(value)) {
        stop("Usage: gintervals.attr.set(intervals.set, attr, value)", call. = FALSE)
    }
    .gcheckroot()

    table <- data.frame(value, stringsAsFactors = FALSE)
    colnames(table) <- attr
    rownames(table) <- intervals.set

    gintervals.attr.import(table, remove.others = FALSE)
    retv <- 0 # suppress return value
}


#' Returns interval set attributes values
#'
#' Returns interval set attributes values.
#'
#' This function returns a data frame that contains interval set attribute
#' values. Column names of the data frame consist of the attribute names, row
#' names contain the interval set names.
#'
#' The list of required interval sets is specified by 'intervals.set'
#' argument. If 'intervals.set' is 'NULL' the attribute values of all existing
#' interval sets are returned.
#'
#' Likewise the list of required attributes is controlled by 'attrs' argument.
#' If 'attrs' is 'NULL' all attribute values of the specified interval sets are
#' returned. The columns are sorted then by "popularity" of an attribute, i.e.
#' the number of interval sets containing this attribute. This sorting is not
#' applied if 'attrs' is not 'NULL'.
#'
#' Empty character string in a table cell marks a non-existing attribute.
#'
#' @param intervals.set a vector of interval set names or 'NULL'
#' @param attrs a vector of attribute names or 'NULL'
#' @return A data frame containing interval set attribute values.
#' @seealso \code{\link{gintervals.attr.import}},
#' \code{\link{gintervals.attr.get}}, \code{\link{gintervals.attr.set}}
#' @keywords ~attr ~attribute ~intervals
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gintervals.attr.set("annotations1", "test_attr", "value1")
#' gintervals.attr.export()
#' gintervals.attr.export(intervals.set = "annotations1")
#' gintervals.attr.export(attrs = "test_attr")
#' gintervals.attr.set("annotations1", "test_attr", "")
#'
#' @export gintervals.attr.export
gintervals.attr.export <- function(intervals.set = NULL, attrs = NULL) {
    .gcheckroot()

    if (!is.null(intervals.set)) {
        intervals.set <- unique(intervals.set)
    }
    if (!is.null(attrs)) {
        attrs <- unique(attrs)
    }

    if (is.null(intervals.set)) {
        intervals.set <- get("GINTERVS", envir = .misha)
    } else {
        idx <- which(!(intervals.set %in% get("GINTERVS", envir = .misha)))[1]
        if (!is.na(idx)) {
            stop(sprintf("Intervals set %s does not exist", intervals.set[idx]), call. = FALSE)
        }
    }

    # Read all attributes for each interval set
    all_attrs <- list()
    for (iset in intervals.set) {
        path <- .gintervals.attr_path(iset)
        all_attrs[[iset]] <- .gintervals.attr_read(path)
    }

    # Collect all attribute names
    if (is.null(attrs)) {
        all_names <- unlist(lapply(all_attrs, names))
        if (length(all_names) == 0) {
            # No attributes at all - return empty data frame with correct row names
            result <- data.frame(row.names = intervals.set)
            return(result)
        }
        # Sort by popularity (most common first)
        name_counts <- table(all_names)
        attrs_sorted <- names(sort(name_counts, decreasing = TRUE))
    } else {
        attrs_sorted <- attrs
    }

    if (length(attrs_sorted) == 0) {
        result <- data.frame(row.names = intervals.set)
        return(result)
    }

    # Build result data frame
    result <- data.frame(
        matrix("", nrow = length(intervals.set), ncol = length(attrs_sorted)),
        stringsAsFactors = FALSE
    )
    colnames(result) <- attrs_sorted
    rownames(result) <- intervals.set

    for (iset in intervals.set) {
        iattrs <- all_attrs[[iset]]
        for (a in attrs_sorted) {
            if (a %in% names(iattrs)) {
                result[iset, a] <- iattrs[a]
            }
        }
    }

    result
}


#' Imports interval set attributes values
#'
#' Imports interval set attributes values.
#'
#' This function imports attribute values contained in a data frame 'table'.
#' The format of a table is similar to the one returned by
#' 'gintervals.attr.export'. The values of the table must be character strings.
#' Column names of the table should specify the attribute names, while row
#' names should contain the interval set names.
#'
#' The specified attributes of the specified interval sets are modified. If an
#' attribute value is an empty string this attribute is removed from the
#' interval set.
#'
#' If 'remove.others' is 'TRUE' all attributes that do not appear in the table
#' are removed, otherwise they are preserved unchanged.
#'
#' @param table a data frame containing attribute values
#' @param remove.others specifies what to do with the attributes that are not
#' in the table
#' @return None.
#' @seealso \code{\link{gintervals.attr.export}},
#' \code{\link{gintervals.attr.set}}, \code{\link{gintervals.attr.get}}
#' @keywords ~attr ~attribute ~intervals
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' t <- data.frame(
#'     myattr = "val1", row.names = "annotations1",
#'     stringsAsFactors = FALSE
#' )
#' gintervals.attr.import(t)
#' gintervals.attr.export(attrs = "myattr")
#'
#' # roll-back the changes
#' t$myattr <- ""
#' gintervals.attr.import(t)
#'
#' @export gintervals.attr.import
gintervals.attr.import <- function(table = NULL, remove.others = FALSE) {
    if (is.null(table)) {
        stop("Usage: gintervals.attr.import(table, remove.others = FALSE)", call. = FALSE)
    }
    .gcheckroot()

    isets <- rownames(table)
    attrs <- colnames(table)

    if (!is.data.frame(table) || any(dim(table) < 1) || !length(isets) || !length(attrs) || any(is.na(isets)) || any(is.na(attrs)) || any(attrs == "")) {
        stop("Invalid format of attributes table", call. = FALSE)
    }

    idx <- which(!(isets %in% get("GINTERVS", envir = .misha)))[1]
    if (!is.na(idx)) {
        stop(sprintf("Intervals set %s does not exist", isets[idx]), call. = FALSE)
    }

    idx <- which(duplicated(isets))[1]
    if (!is.na(idx)) {
        stop(sprintf("Intervals set %s appears more than once", isets[idx]), call. = FALSE)
    }

    idx <- which(duplicated(attrs))[1]
    if (!is.na(idx)) {
        stop(sprintf("Attribute %s appears more than once", attrs[idx]), call. = FALSE)
    }

    # Coerce all columns to character
    table[, 1:ncol(table)] <- sapply(table[, 1:ncol(table)], as.character)

    for (iset in isets) {
        .gintervals.attr_check_writable(iset)

        path <- .gintervals.attr_path(iset)
        existing <- .gintervals.attr_read(path)

        if (remove.others) {
            # Start fresh - only keep what's in the table
            new_attrs <- character(0)
        } else {
            new_attrs <- existing
        }

        # Apply updates from table
        for (a in attrs) {
            val <- table[iset, a]
            if (is.na(val) || val == "") {
                # Remove attribute
                new_attrs <- new_attrs[names(new_attrs) != a]
            } else {
                new_attrs[a] <- val
            }
        }

        .gintervals.attr_write(path, new_attrs)
    }

    retv <- 0 # suppress return value
}
