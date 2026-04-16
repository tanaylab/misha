#' Extract GLM feature matrix from motif energy tracks
#'
#' Reads multiple motif energy tracks, applies scaling and logistic transforms,
#' and returns the complete feature matrix in a single C++ pass. This replaces
#' the R-side pipeline of gextract -> scale -> reshape -> logistic transform.
#'
#' @param track_names Character vector of motif track names
#'   (e.g., \code{"th_epi.all_db_motifs.JASPAR_Zic2"}).
#' @param intervals Data frame with \code{chrom}, \code{start}, \code{end} columns.
#'   All intervals must have the same width.
#' @param tile_size Tile width in bp (default 200).
#' @param flank_size Flanking region in bp added to each side of the peak (default 350).
#' @param max_cap Named numeric vector of per-motif maximum cap values
#'   (genome-wide quantiles). Names must match \code{track_names}.
#' @param dis_from_cap Distance below cap for normalization (default 10).
#' @param scale_factor Output range scale factor (default 10).
#' @param transforms List of transform configurations. Each element is a list
#'   with fields \code{L}, \code{k}, \code{x_0}, \code{pre_shift}, \code{post_shift}.
#'   Defaults to the 4 standard logistic heads (low-energy, high-energy, sigmoid,
#'   higher-energy).
#' @param gc_track GC content track name (default \code{"seq.G_or_C"}).
#' @param gc_scale_factor Scale factor for GC features (default 10).
#'
#' @return A numeric matrix with one row per interval and columns for:
#'   \itemize{
#'     \item Motif features: \code{n_motifs * n_tiles * n_transforms} columns
#'     \item GC features: \code{n_tiles} columns
#'     \item GC interactions: \code{choose(n_tiles, 2)} columns
#'   }
#'
#' @examples
#' \dontrun{
#' gdb.init("/path/to/trackdb")
#' peaks <- data.frame(chrom = "chr1", start = c(1000, 2000), end = c(1300, 2300))
#' motifs <- c("th_epi.all_db_motifs.JASPAR_Zic2")
#' caps <- c(JASPAR_Zic2 = 10.5)
#' names(caps) <- motifs
#' result <- glm_extract_features(motifs, peaks, max_cap = caps)
#' }
#'
#' @export
glm_extract_features <- function(
  track_names,
  intervals,
  tile_size = 200L,
  flank_size = 350L,
  max_cap,
  dis_from_cap = 10,
  scale_factor = 10,
  transforms = list(
      list(L = 2, k = 0.5, x_0 = 0, pre_shift = 0, post_shift = -1),
      list(L = 2, k = 0.5, x_0 = 10, pre_shift = 0, post_shift = 0),
      list(L = 1, k = 1, x_0 = 0, pre_shift = -5, post_shift = 0),
      list(L = 2, k = 1, x_0 = 10, pre_shift = 0, post_shift = 0)
  ),
  gc_track = "seq.G_or_C",
  gc_scale_factor = 10
) {
    # Validate inputs
    stopifnot(is.character(track_names), length(track_names) > 0)
    stopifnot(is.data.frame(intervals))
    stopifnot(all(c("chrom", "start", "end") %in% names(intervals)))
    stopifnot(length(max_cap) == length(track_names))

    n_peaks <- nrow(intervals)
    n_motifs <- length(track_names)
    n_transforms <- length(transforms)
    if (n_peaks == 0L) {
        stop("'intervals' must contain at least one row")
    }

    # Check uniform peak size
    peak_sizes <- intervals$end - intervals$start
    if (length(unique(peak_sizes)) != 1) {
        stop("All intervals must have the same width")
    }

    # Convert chromosome names to 0-based IDs
    chrom_sizes <- gintervals.chrom_sizes(intervals)
    chrom_ids <- match(as.character(intervals$chrom), chrom_sizes$chrom) - 1L

    if (any(is.na(chrom_ids))) {
        stop("Some chromosome names not found in the genome database")
    }

    # Order max_cap to match track_names
    if (!is.null(names(max_cap))) {
        max_cap_ordered <- max_cap[track_names]
        if (any(is.na(max_cap_ordered))) {
            stop("max_cap names do not match track_names")
        }
    } else {
        max_cap_ordered <- max_cap
    }

    # Build transform matrix (n_transforms x 5, column-major)
    transform_mat <- matrix(0, nrow = n_transforms, ncol = 5)
    for (i in seq_along(transforms)) {
        t <- transforms[[i]]
        transform_mat[i, 1] <- t$L
        transform_mat[i, 2] <- t$k
        transform_mat[i, 3] <- t$x_0
        transform_mat[i, 4] <- t$pre_shift
        transform_mat[i, 5] <- t$post_shift
    }

    # Call C++
    result <- .gcall(
        "C_glm_extract_features",
        as.character(track_names),
        as.integer(chrom_ids),
        as.numeric(intervals$start),
        as.numeric(intervals$end),
        as.integer(tile_size),
        as.integer(flank_size),
        as.numeric(max_cap_ordered),
        as.numeric(dis_from_cap),
        as.numeric(scale_factor),
        transform_mat,
        as.character(gc_track),
        as.numeric(gc_scale_factor),
        .misha_env()
    )

    # Compute column names
    peak_size <- peak_sizes[1]
    extended <- peak_size + 2 * flank_size
    n_tiles <- extended %/% tile_size

    # Tile distances relative to peak center
    half <- peak_size %/% 2
    tile_starts <- (seq_len(n_tiles) - 1L) * tile_size - flank_size - half
    tile_midpoints <- tile_starts + tile_size / 2
    dist_labels <- as.character(as.integer(tile_midpoints))

    # Short motif names (last component of track path)
    motif_short <- vapply(
        strsplit(track_names, "\\."),
        function(x) x[length(x)],
        character(1)
    )

    # Transform suffixes
    transform_names <- c("low-energy", "high-energy", "sigmoid", "higher-energy")
    if (n_transforms != 4 || !identical(
        vapply(transforms, function(t) t$L, numeric(1)),
        c(2, 2, 1, 2)
    )) {
        transform_names <- paste0("t", seq_len(n_transforms))
    }

    # Motif columns: (tile * n_motifs + motif) * n_transforms + transform
    motif_col_names <- character(n_motifs * n_tiles * n_transforms)
    idx <- 1
    for (ti in seq_len(n_tiles)) {
        for (mi in seq_len(n_motifs)) {
            for (tr in seq_len(n_transforms)) {
                motif_col_names[idx] <- paste0(
                    motif_short[mi], "_dist_",
                    gsub("-", "neg_", dist_labels[ti]),
                    "_", transform_names[tr]
                )
                idx <- idx + 1
            }
        }
    }

    # GC columns
    gc_col_names <- paste0(
        "gc_content_tile_dist_",
        gsub("-", "neg_", dist_labels)
    )

    # GC interaction columns (choose(n_tiles, 2) entries in row-major order)
    n_gc_inter <- if (n_tiles >= 2L) n_tiles * (n_tiles - 1L) %/% 2L else 0L
    gc_inter_names <- character(n_gc_inter)
    if (n_gc_inter > 0L) {
        idx <- 1L
        for (a in seq_len(n_tiles - 1L)) {
            for (b in seq.int(a + 1L, n_tiles)) {
                gc_inter_names[idx] <- paste0(
                    "int_", gc_col_names[a], "_x_", gc_col_names[b]
                )
                idx <- idx + 1L
            }
        }
    }

    colnames(result) <- c(motif_col_names, gc_col_names, gc_inter_names)
    result
}


