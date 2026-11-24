# Import of figwig files (.fw, https://github.com/jmschrei/figwig)

# The text a value of 'select' is compared with: figwig attribute values are text as they become
# misha attributes (a string as it is, anything else as its JSON text)
.figwig_filter_text <- function(x) {
    if (is.factor(x)) {
        x <- as.character(x)
    }
    if (is.character(x)) {
        return(x)
    }
    if (is.logical(x)) {
        return(ifelse(x, "true", "false"))
    }
    if (is.numeric(x)) {
        return(ifelse(is.finite(x) & x == round(x) & abs(x) < 1e15, sprintf("%.0f", x), as.character(x)))
    }
    as.character(x)
}

# The misha track name a figwig track name becomes: a trailing "+" or "-" (a strand) is spelled
# "plus" or "minus", each run of characters other than ASCII letters, digits and "_" (dots
# included) becomes one "_", and a name that does not start with a letter gets an "X" in front,
# unless a prefix comes before it. PyMisha's _sanitize_name.
.figwig_sanitize <- function(name, has_prefix) {
    name <- ifelse(grepl("\\+$", name, perl = TRUE), sub("\\+$", "plus", name, perl = TRUE), sub("-$", "minus", name, perl = TRUE))
    name <- gsub("[^A-Za-z0-9_]+", "_", name, perl = TRUE)
    if (!has_prefix) {
        name <- ifelse(grepl("^[A-Za-z]", name, perl = TRUE), name, paste0("X", name))
    }
    name
}

# The indices (1-based) of the figwig tracks to import
.figwig_select <- function(info, tracks, select, file) {
    n <- length(info$name)
    if (is.null(tracks)) {
        idx <- seq_len(n)
    } else {
        if (is.character(tracks)) {
            idx <- match(tracks, info$name)
            if (anyNA(idx)) {
                stop(sprintf(
                    "%s has no track(s) named %s; its tracks: %s%s", file,
                    paste0("\"", tracks[is.na(idx)], "\"", collapse = ", "),
                    paste(utils::head(info$name, 20), collapse = ", "), if (n > 20) " ..." else ""
                ), call. = FALSE)
            }
        } else if (is.numeric(tracks) && all(is.finite(tracks)) && all(tracks == round(tracks))) {
            idx <- as.integer(tracks)
            if (any(idx < 1 | idx > n)) {
                stop(sprintf("Track index %d is out of range: %s holds %d track(s)", idx[idx < 1 | idx > n][1], file, n), call. = FALSE)
            }
        } else {
            stop("tracks must be figwig track names or (1-based) indices", call. = FALSE)
        }
        if (anyDuplicated(idx)) {
            stop("tracks names a figwig track more than once", call. = FALSE)
        }
        if (!length(idx)) {
            stop("No figwig tracks to import", call. = FALSE)
        }
    }

    if (length(select)) {
        if (!is.list(select) && !is.atomic(select) || is.null(names(select)) || any(names(select) == "")) {
            stop("select must be a named list of attribute values", call. = FALSE)
        }
        # the text of each track's attribute 'key' (NA where it has none or it is null)
        text_of <- function(key) {
            if (key == "name") {
                return(info$name)
            }
            if (key == "dtype") {
                return(info$dtype)
            }
            if (key == "group") {
                return(as.character(info$group))
            }
            vapply(info$attrs, function(a) if (key %in% names(a)) a[[key]] else NA_character_, character(1))
        }
        keys <- unique(c("name", "dtype", "group", unlist(lapply(info$attrs, names))))
        for (key in names(select)) {
            if (!key %in% keys) {
                stop(sprintf("No track of %s has the attribute \"%s\"; its attributes: %s", file, key, paste(sort(keys), collapse = ", ")), call. = FALSE)
            }
            values <- text_of(key)
            wanted <- .figwig_filter_text(unlist(select[[key]], use.names = FALSE))
            idx <- idx[!is.na(values[idx]) & values[idx] %in% wanted]
        }
        if (!length(idx)) {
            stop(sprintf("No track of %s matches 'select'", file), call. = FALSE)
        }
    }
    idx
}

# The attributes a track takes from its figwig record (md, a named character vector, NA for null):
# all but those the import sets (created.by, created.date, created.user), type and binsize,
# description unless the description given is empty, and the names in keep (attributes the caller
# set). Null and empty values are left out.
.figwig_restore_attrs <- function(track, md, description, keep = NULL) {
    skip <- c("created.by", "created.date", "created.user", "type", "binsize", if (nzchar(description)) "description", keep)
    md <- md[!is.na(names(md)) & nzchar(names(md)) & !names(md) %in% skip & !is.na(md) & nzchar(md)]
    md <- md[!duplicated(names(md))]
    if (length(md)) {
        table <- as.data.frame(as.list(md), stringsAsFactors = FALSE, check.names = FALSE)
        rownames(table) <- track
        .gtrack.attr.import(table, FALSE, TRUE)
    }
    invisible(NULL)
}


#' Creates tracks from the tracks of a figwig file
#'
#' Creates one misha track per selected track of a figwig file ('.fw', the
#' multi-track format of the figwig package, \url{https://github.com/jmschrei/figwig}).
#' The file is read by misha itself; neither Python nor figwig is needed.
#'
#' With 'binsize' 0 (or NULL) each track is 'Sparse': one interval per maximal
#' run of bases with the same value, leaving out the bases whose value is 0 or
#' NaN. That suits tracks that are mostly zero, such as counts or cut sites; a
#' dense float track makes one interval per base. With a positive 'binsize'
#' each track is 'Dense', its bins computed as \code{\link{gtrack.import}}
#' computes the bins of a bigWig: the mean of the values of the bin's bases,
#' each run of equal values weighted by its length. A figwig file has no "no
#' data": a base that was never written is 0, as is one written as 0 (figwig's
#' 'bw2fw' writes a bigWig's uncovered bases as 0 unless given '--missing nan').
#' NaN bases, bases past the end of the file's chromosome and chromosomes the
#' file does not have have no value: they count as 'defval' unless it is NaN,
#' and a bin with no value is NaN.
#'
#' Tracks are chosen with 'tracks' (figwig track names, or 1-based indices),
#' in the order given, or all of them in the file's order, and then with
#' 'select', a named list of attribute values (a vector of several values
#' matches any of them): a track is imported when every attribute given
#' matches. Values are compared as text, as they become attributes (below),
#' so \code{list(rep = 1)} and \code{list(rep = "1")} both match a
#' \code{"rep": 1} record. 'name', 'dtype' and 'group' can be selected on too.
#'
#' The track made from figwig track 'n' is named '[track.prefix]n2', where 'n2'
#' is the entry of 'names' for it if given (a character vector, one per
#' imported track, or one named by figwig track names) and otherwise 'n'
#' sanitized: a trailing '+' or '-' (a strand) is spelled 'plus' or 'minus',
#' each run of characters other than ASCII letters, digits and '_' becomes one
#' '_', and without a prefix a name that does not start with a letter gets an
#' 'X' in front ('atac.+' becomes 'atac_plus', 'H3K27ac-rep1' 'H3K27ac_rep1').
#' Dots in 'track.prefix' are directories, which must exist. All the names are
#' checked before anything is written, and existing tracks are not overwritten.
#'
#' Each figwig track's attributes become attributes of its track, except
#' 'name', 'dtype' and 'group'; 'created.by', 'created.date' and
#' 'created.user', which the import sets; 'type' and 'binsize'; and
#' 'description', unless the 'description' given is not empty. A string
#' attribute is kept as it is, a null one is left out, and a number, boolean,
#' list or object becomes its JSON text. 'attrs' win over them.
#'
#' Chromosome names are matched through the chromosome alias chain as in
#' \code{\link{gtrack.import}}, with the same reporting of names the database
#' does not have; a chromosome whose length differs from the database's is
#' warned about. The import is all or nothing: if it fails or is interrupted,
#' none of the tracks is created.
#'
#' Every layout of the figwig format version 1 is read: position and track
#' order, stored or zlib-compressed chunks, dense, sparse, runs and matrix
#' chunks and every element type. PyMisha's \code{gtrack_import_figwig} takes
#' the same arguments and makes the same tracks.
#'
#' @param file path of the figwig file
#' @param binsize bin size of 'Dense' tracks, or 0 (or NULL) for 'Sparse' tracks
#' @param track.prefix prefix of the track names
#' @param tracks figwig track names or 1-based indices to import, or NULL for all
#' @param select a named list of attribute values to choose tracks by
#' @param names names of the created tracks (before the prefix), instead of the sanitized figwig names
#' @param defval value of the bases with no value in 'Dense' tracks
#' @param description description of every created track; when empty, a figwig 'description' attribute is kept
#' @param attrs a named vector or list of attributes set on every created track, over those of the file
#' @return Names of the created tracks, invisibly.
#' @seealso \code{\link{gtrack.import}}, \code{\link{gtrack.import_parquet}},
#' \code{\link{gtrack.attr.get}}
#' @keywords ~figwig ~track
#' @examples
#' \dontrun{
#' # every track of the file, as sparse tracks named after the figwig tracks
#' gtrack.import_figwig("tracks.fw")
#'
#' # the ATAC-seq tracks, as dense tracks of 20 bp bins
#' gtrack.import_figwig("tracks.fw",
#'     binsize = 20, track.prefix = "fw_",
#'     select = list(assay = "ATAC-seq")
#' )
#' }
#'
#' @export gtrack.import_figwig
gtrack.import_figwig <- function(file, binsize = 0, track.prefix = NULL, tracks = NULL, select = NULL, names = NULL,
                                 defval = NaN, description = "", attrs = NULL) {
    if (missing(file)) {
        stop("Usage: gtrack.import_figwig(file, binsize = 0, track.prefix = NULL, tracks = NULL, select = NULL, names = NULL, defval = NaN, description = \"\", attrs = NULL)", call. = FALSE)
    }
    .gcheckroot()

    if (!is.character(file) || length(file) != 1 || is.na(file)) {
        stop("file must be a single string", call. = FALSE)
    }
    if (is.null(binsize)) {
        binsize <- 0
    }
    if (!is.numeric(binsize) || length(binsize) != 1 || !is.finite(binsize) || binsize < 0 || binsize != round(binsize) || binsize >= 2^32) {
        stop("binsize must be 0 or a positive integer below 2^32", call. = FALSE)
    }
    if (!is.numeric(defval) || length(defval) != 1) {
        stop("defval must be a number", call. = FALSE)
    }
    if (is.null(track.prefix)) {
        track.prefix <- ""
    }
    if (!is.character(track.prefix) || length(track.prefix) != 1 || is.na(track.prefix)) {
        stop("track.prefix must be a string", call. = FALSE)
    }
    if (!is.character(description) || length(description) != 1 || is.na(description)) {
        stop("description must be a string", call. = FALSE)
    }
    if (!is.null(attrs) && (is.null(base::names(attrs)) || any(base::names(attrs) == ""))) {
        stop("attrs must be a named vector or list", call. = FALSE)
    }
    if (!file.exists(file) || dir.exists(file)) {
        stop(sprintf("File %s does not exist", file), call. = FALSE)
    }

    info <- .gcall("figwig_info", file, .misha_env())
    idx <- .figwig_select(info, tracks, select, file)
    fw_names <- info$name[idx]

    if (is.null(names)) {
        base <- .figwig_sanitize(fw_names, nzchar(track.prefix))
    } else {
        if (!is.character(names) || anyNA(names)) {
            stop("names must be a character vector", call. = FALSE)
        }
        if (!is.null(base::names(names))) {
            unknown <- setdiff(base::names(names), fw_names)
            if (length(unknown)) {
                stop(sprintf("names maps figwig track(s) that are not imported: %s", paste0("\"", unknown, "\"", collapse = ", ")), call. = FALSE)
            }
            base <- ifelse(fw_names %in% base::names(names), names[fw_names], .figwig_sanitize(fw_names, nzchar(track.prefix)))
            base <- unname(base)
        } else {
            if (length(names) != length(idx)) {
                stop(sprintf("names has %d element(s) but %d track(s) are imported", length(names), length(idx)), call. = FALSE)
            }
            base <- names
        }
    }
    new_tracks <- paste0(track.prefix, base)

    bad <- !grepl("^[A-Za-z][A-Za-z0-9_.]*$", new_tracks, perl = TRUE)
    if (any(bad)) {
        stop(sprintf(
            "Invalid track name(s): %s (a letter, then letters, digits, _ or .); choose other names with 'names' or 'track.prefix'",
            paste0("\"", new_tracks[bad], "\"", collapse = ", ")
        ), call. = FALSE)
    }
    dups <- unique(new_tracks[duplicated(new_tracks)])
    if (length(dups)) {
        clash <- vapply(dups, function(t) sprintf("%s from %s", t, paste0("\"", fw_names[new_tracks == t], "\"", collapse = ", ")), character(1))
        stop(sprintf("Figwig tracks would make the same track name: %s. Name them with 'names'.", paste(clash, collapse = "; ")), call. = FALSE)
    }
    for (track in new_tracks) {
        .gconfirmtrackcreate(track)
        if (!dir.exists(dirname(.track_dir(track)))) {
            stop(sprintf("Directory %s does not exist", dirname(gsub(".", "/", track, fixed = TRUE))), call. = FALSE)
        }
    }

    finals <- vapply(new_tracks, .track_dir, character(1), USE.NAMES = FALSE)
    staging <- file.path(dirname(finals), vapply(basename(finals), .gdb.staging_name, character(1), USE.NAMES = FALSE))
    created <- rep(FALSE, length(new_tracks))
    success <- FALSE
    tryCatch(
        {
            for (parent in unique(dirname(finals))) {
                try(.gdb.trash_sweep_old(parent), silent = TRUE)
            }
            for (stage_dir in staging) {
                if (!.gwith_umask(dir.create(stage_dir, mode = "0777"))) {
                    stop(sprintf("Failed to create a directory %s", stage_dir), call. = FALSE)
                }
            }

            .gcall("gtrackimportfigwig", staging, file, file, as.integer(idx - 1L), binsize, defval, .misha_env())

            for (i in seq_along(new_tracks)) {
                .gdb.fsync(staging[i], recursive = TRUE)
                if (dir.exists(finals[i])) {
                    stop(sprintf("Refusing to overwrite existing %s", finals[i]), call. = FALSE)
                }
                created[i] <- file.rename(staging[i], finals[i])
                if (!created[i]) {
                    stop(sprintf("Failed to rename %s -> %s", staging[i], finals[i]), call. = FALSE)
                }
                .gdb.fsync(dirname(finals[i]))
                .gdb.invalidate_dir_cache(finals[i])
                .gdb.add_track(new_tracks[i])
                .gtrack_set_created_attrs(
                    new_tracks[i], description,
                    sprintf(
                        "gtrack.import_figwig(\"%s\", %.0f, \"%s\", \"%s\", %g)",
                        file, binsize, track.prefix, fw_names[i], defval
                    ),
                    attrs
                )
                .figwig_restore_attrs(new_tracks[i], info$attrs[[idx[i]]], description, base::names(attrs))
                if (.gdb.is_indexed()) {
                    gtrack.convert_to_indexed(new_tracks[i])
                }
            }
            success <- TRUE
        },
        finally = {
            if (!success) {
                for (stage_dir in staging[dir.exists(staging)]) {
                    .gdb.trash(stage_dir)
                }
                for (i in which(created)) {
                    if (!.gdb.trash(finals[i]) && dir.exists(finals[i])) {
                        warning(sprintf(
                            "Track %s post-create cleanup left residue at %s; manual cleanup required",
                            new_tracks[i], finals[i]
                        ), call. = FALSE)
                    }
                    try(.gdb.rm_track(new_tracks[i]), silent = TRUE)
                }
            }
        }
    )
    invisible(new_tracks)
}
