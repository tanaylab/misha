# ============================================================================
# GLM Predictor Virtual Track
# ============================================================================
# Fused GLM linear prediction: per-position computation of
# bias + Σ(weight × transform(scale(smooth(track)))) + Σ(inter_weight × transform(product))

# Cut() labels for a break vector - same labels glm_pred uses for stratum binning
# (matches BinFinder::val2bin with include_lowest=true, right=true on the C++ side).
.glm_pred_cut_labels <- function(breaks) {
    midpoints <- breaks[-length(breaks)] + diff(breaks) / 2
    levels(cut(midpoints, breaks = breaks, include.lowest = TRUE, right = TRUE))
}

# Validate / auto-derive dimnames for a stratified weights / bias / interaction array.
# x: numeric array. Plain numeric vector accepted only when M=1 and there is no leading axis (bias).
# K_per: integer vector of expected trailing-axis lengths (length M).
# selector_tracks: character vector of expected trailing-axis names (length M).
# selector_breaks: list of break vectors (length M) - labels are derived from these.
# leading_dim: integer expected leading-axis length, or NA when there is no leading axis.
# leading_name: character(1) name to assign to leading axis when present and unnamed.
# param: character(1) - parameter name used in error messages.
.glm_pred_normalize_strata_array <- function(x, K_per, selector_tracks, selector_breaks,
                                             leading_dim, leading_name, param) {
    M <- length(selector_tracks)
    has_leading <- !is.na(leading_dim)
    expected_dim <- if (has_leading) c(leading_dim, K_per) else K_per

    if (!is.numeric(x)) {
        stop(sprintf("'%s' must be numeric", param), call. = FALSE)
    }

    actual_dim <- dim(x)
    if (is.null(actual_dim)) {
        # Plain numeric vector. Allowed only for bias (no leading axis) at M=1.
        if (has_leading) {
            stop(sprintf(
                "'%s' must be an array with dim = c(%s); got a length-%d numeric vector. To migrate from a flat matrix, use array(%s, dim = c(%s)) instead of matrix(%s, nrow = %d).",
                param, paste(expected_dim, collapse = ", "), length(x),
                param, paste(expected_dim, collapse = ", "),
                param, expected_dim[1L]
            ), call. = FALSE)
        }
        if (M != 1L) {
            stop(sprintf(
                "'%s' must be an array with dim = c(%s); got a length-%d numeric vector.",
                param, paste(expected_dim, collapse = ", "), length(x)
            ), call. = FALSE)
        }
        if (length(x) != K_per[1L]) {
            stop(sprintf(
                "'%s' must be numeric of length %d (or scalar); got length %d.",
                param, K_per[1L], length(x)
            ), call. = FALSE)
        }
        dim(x) <- K_per
    } else if (length(actual_dim) != length(expected_dim) || any(actual_dim != expected_dim)) {
        migration_hint <- if (M >= 2L && has_leading && length(actual_dim) == 2L) {
            sprintf(
                " To migrate from a flat matrix, use array(%s, dim = c(%s)) instead of matrix(%s, nrow = %d).",
                param, paste(expected_dim, collapse = ", "), param, expected_dim[1L]
            )
        } else {
            ""
        }
        stop(sprintf(
            "'%s' must be an array with dim = c(%s); got dim = c(%s).%s",
            param, paste(expected_dim, collapse = ", "),
            paste(actual_dim, collapse = ", "), migration_hint
        ), call. = FALSE)
    }

    dn <- dimnames(x)
    if (is.null(dn)) dn <- vector("list", length(expected_dim))
    dn_names <- names(dn)
    if (is.null(dn_names)) dn_names <- character(length(expected_dim))

    for (m in seq_len(M)) {
        ax <- if (has_leading) m + 1L else m
        expected_labs <- .glm_pred_cut_labels(selector_breaks[[m]])
        if (is.null(dn[[ax]])) {
            dn[[ax]] <- expected_labs
        } else if (!identical(as.character(dn[[ax]]), expected_labs)) {
            stop(sprintf(
                "dimnames(%s)[[%d]] does not match cut() labels for selector_breaks[[%d]]; expected c(%s), got c(%s). Either omit dimnames (auto-derived) or match the cut() labels exactly.",
                param, ax, m,
                paste(sprintf('"%s"', expected_labs), collapse = ", "),
                paste(sprintf('"%s"', dn[[ax]]), collapse = ", ")
            ), call. = FALSE)
        }
        if (is.na(dn_names[ax]) || dn_names[ax] == "") {
            dn_names[ax] <- selector_tracks[m]
        } else if (dn_names[ax] != selector_tracks[m]) {
            stop(sprintf(
                "names(dimnames(%s))[%d] must be '%s'; got '%s'. Either omit dimnames (auto-derived) or match the selector_tracks name exactly.",
                param, ax, selector_tracks[m], dn_names[ax]
            ), call. = FALSE)
        }
    }
    if (has_leading) {
        if (is.na(dn_names[1L]) || dn_names[1L] == "") {
            dn_names[1L] <- leading_name
        }
    }

    names(dn) <- dn_names
    dimnames(x) <- dn
    storage.mode(x) <- "double"
    x
}

#' Create a GLM predictor virtual track
#'
#' Creates a virtual track that computes a fused generalized linear model
#' prediction at each genome position. The pipeline for each entry is:
#' smooth → scale (cap + normalize) → transform (logistic) → weight.
#' Interactions are computed post-scaling, pre-transform.
#'
#' @param name character(1) Virtual track name
#' @param tracks character(N) Genomic track names (repeated OK)
#' @param inner_func character(N) \code{"sum"} or \code{"lse"} per entry
#' @param weights LM coefficients per entry. Shape depends on selectors:
#'   \itemize{
#'     \item No selector: \code{numeric(N)}.
#'     \item One selector (M=1): \code{matrix(N, K_1)}.
#'     \item Two or more selectors (M>=2): array with
#'       \code{dim = c(N, K_1, ..., K_M)}. Plain
#'       \code{N x prod(K_m)} matrices are rejected to avoid flatten-order
#'       ambiguity. Build with \code{array(coefs, dim = c(N, K_1, ..., K_M))}.
#'   }
#'   When a selector is specified, trailing-axis dimnames are auto-derived
#'   from \code{selector_breaks} (\code{cut(include.lowest=TRUE, right=TRUE)}
#'   labels) and trailing-axis names from \code{selector_tracks}. If you
#'   supply your own dimnames they must match exactly; mismatches error
#'   with the offending axis named.
#' @param bias Intercept term (default 0). Shape depends on selectors:
#'   \itemize{
#'     \item No selector: \code{numeric(1)}.
#'     \item One selector (M=1): scalar (recycled) or \code{numeric(K_1)}.
#'     \item Two or more selectors (M>=2): scalar (recycled to all strata)
#'       or array with \code{dim = c(K_1, ..., K_M)}.
#'   }
#'   Auto-derived dimnames mirror those of \code{weights} on the trailing
#'   axes.
#' @param kernels list of numeric vectors (N or 1, recycled) Sub-bin kernel
#'   weights, or NULL for direct aggregation
#' @param kernel_bins numeric(B) Offsets relative to shift window center, or NULL
#' @param shifts list of numeric(2) Per-entry \code{c(sshift, eshift)}, or NULL
#' @param trans_family character(N or 1 or NULL) \code{"logist"} or NA per entry
#' @param trans_params list of lists (N or NULL) Each element is a list with
#'   fields \code{L}, \code{k}, \code{x_0}, and optionally \code{pre_shift},
#'   \code{post_shift}
#' @param max_cap numeric(N or NULL) Capping threshold per entry (NA to skip)
#' @param dis_from_cap numeric(N or NULL) Distance from cap per entry (NA to skip)
#' @param simple_cap logical(N or NULL) If TRUE, use simple cap-and-divide
#'   scaling: \code{min(raw, max_cap) / dis_from_cap}. Default NULL (standard
#'   cap-normalize scaling).
#' @param scale_factor numeric(1) Global scale factor (default 10)
#' @param interactions list of integer(2) (M or NULL) Pairs of entry indices (1-based)
#' @param interaction_weights Per-interaction LM coefficients. Shape mirrors
#'   \code{weights} with the leading axis being the interaction index:
#'   \itemize{
#'     \item No selector: \code{numeric(M_int)}.
#'     \item One selector (M=1): \code{matrix(M_int, K_1)}.
#'     \item Two or more selectors (M>=2): array with
#'       \code{dim = c(M_int, K_1, ..., K_M)}.
#'   }
#'   Same dimnames rules as \code{weights}; leading axis name is
#'   \code{"interaction"}.
#' @param interaction_trans_family character(M or 1 or NULL) Transform for interactions
#' @param inter_trans_params list of lists (M or NULL) Logistic params per interaction
#' @param selector_tracks character(M) or NULL Names of fixed-bin (dense) tracks
#'   used for per-position model selection. With M selectors, strata are the
#'   Cartesian product of per-selector bins (column-major: first selector varies
#'   fastest). When NULL (default), K = 1 and a single set of weights is applied
#'   everywhere.
#' @param selector_breaks list of M numeric vectors, or NULL. Each element is
#'   the break vector (length K_m + 1, K_m >= 1) for the corresponding selector
#'   track in \code{selector_tracks}. Break points use
#'   \code{cut(include.lowest = TRUE, right = TRUE)} semantics. Required and
#'   length-matched to \code{selector_tracks} when that is non-NULL.
#'
#' @return Invisibly returns \code{name}.
#' @export
glm_pred.create <- function(name,
                            tracks,
                            inner_func,
                            weights,
                            bias = 0,
                            kernels = NULL,
                            kernel_bins = NULL,
                            shifts = NULL,
                            trans_family = NULL,
                            trans_params = NULL,
                            max_cap = NULL,
                            dis_from_cap = NULL,
                            simple_cap = NULL,
                            scale_factor = 10,
                            interactions = NULL,
                            interaction_weights = NULL,
                            interaction_trans_family = NULL,
                            inter_trans_params = NULL,
                            selector_tracks = NULL,
                            selector_breaks = NULL) {
    .gcheckroot()

    # --- Validate core vectors ---
    if (!is.character(name) || length(name) != 1) {
        stop("'name' must be a single character string", call. = FALSE)
    }
    if (!is.character(tracks) || length(tracks) == 0) {
        stop("'tracks' must be a non-empty character vector", call. = FALSE)
    }
    N <- length(tracks)

    # Validate that all tracks exist
    existing <- gtrack.ls()
    missing_tracks <- setdiff(unique(tracks), existing)
    if (length(missing_tracks) > 0) {
        stop(
            sprintf(
                "Track(s) not found: %s",
                paste0("'", missing_tracks, "'", collapse = ", ")
            ),
            call. = FALSE
        )
    }

    if (!is.character(inner_func) || length(inner_func) != N) {
        stop(sprintf("'inner_func' must be a character vector of length %d", N), call. = FALSE)
    }
    invalid_if <- setdiff(unique(inner_func), c("sum", "lse"))
    if (length(invalid_if) > 0) {
        stop(
            sprintf(
                "'inner_func' values must be 'sum' or 'lse' (got: %s)",
                paste0("'", invalid_if, "'", collapse = ", ")
            ),
            call. = FALSE
        )
    }

    # --- Validate selector tracks and determine K_total = prod K_m (Cartesian product) ---
    K <- 1L
    K_per <- integer(0)
    if (!is.null(selector_tracks)) {
        if (is.null(selector_breaks)) {
            stop("'selector_breaks' required when 'selector_tracks' is specified", call. = FALSE)
        }
        if (!is.character(selector_tracks) || length(selector_tracks) < 1) {
            stop("'selector_tracks' must be a character vector of length >= 1", call. = FALSE)
        }
        if (!is.list(selector_breaks) || length(selector_breaks) != length(selector_tracks)) {
            stop(sprintf(
                "'selector_breaks' must be a list of length %d (one per selector)",
                length(selector_tracks)
            ), call. = FALSE)
        }
        existing_tracks <- gtrack.ls()
        K_per <- integer(length(selector_tracks))
        for (m in seq_along(selector_tracks)) {
            tname <- selector_tracks[m]
            br <- selector_breaks[[m]]
            if (!is.numeric(br) || length(br) < 2) {
                stop(sprintf(
                    "'selector_breaks[[%d]]' must be a numeric vector of length >= 2",
                    m
                ), call. = FALSE)
            }
            if (is.unsorted(br, strictly = TRUE)) {
                stop(sprintf("'selector_breaks[[%d]]' must be strictly increasing", m), call. = FALSE)
            }
            if (!(tname %in% existing_tracks)) {
                stop(sprintf("Selector track '%s' not found", tname), call. = FALSE)
            }
            sel_info <- gtrack.info(tname)
            if (sel_info$type != "dense") {
                stop(sprintf(
                    "Selector track '%s' must be a fixed-bin (dense) track, got '%s'",
                    tname, sel_info$type
                ), call. = FALSE)
            }
            K_per[m] <- length(br) - 1L
        }
        K <- as.integer(prod(K_per))
        if (K < 1L) stop("Product of selector bin counts must be >= 1", call. = FALSE)
    }

    # --- Validate weights and bias ---
    if (is.null(selector_tracks)) {
        # No selector: weights = numeric(N), bias = numeric(1).
        if (!is.numeric(weights) || !is.null(dim(weights)) || length(weights) != N) {
            stop(sprintf(
                "'weights' must be a numeric vector of length %d when no selector is specified", N
            ), call. = FALSE)
        }
        if (!is.numeric(bias) || length(bias) != 1L) {
            stop("'bias' must be a numeric scalar when no selector is specified", call. = FALSE)
        }
        weights <- as.numeric(weights)
        bias <- as.numeric(bias)
    } else {
        weights <- .glm_pred_normalize_strata_array(
            weights,
            K_per = K_per, selector_tracks = selector_tracks, selector_breaks = selector_breaks,
            leading_dim = N, leading_name = "entry", param = "weights"
        )
        if (!is.numeric(bias)) {
            stop("'bias' must be numeric", call. = FALSE)
        }
        if (length(bias) == 1L) {
            bias <- array(rep(as.numeric(bias), prod(K_per)), dim = K_per)
        }
        bias <- .glm_pred_normalize_strata_array(
            bias,
            K_per = K_per, selector_tracks = selector_tracks, selector_breaks = selector_breaks,
            leading_dim = NA_integer_, leading_name = NA_character_, param = "bias"
        )
    }

    if (!is.numeric(scale_factor) || length(scale_factor) != 1 || scale_factor <= 0) {
        stop("'scale_factor' must be a positive numeric scalar", call. = FALSE)
    }

    # --- Validate and process shifts ---
    sshifts <- rep(0, N)
    eshifts <- rep(0, N)
    if (!is.null(shifts)) {
        if (!is.list(shifts) || length(shifts) != N) {
            stop(sprintf("'shifts' must be a list of length %d", N), call. = FALSE)
        }
        for (i in seq_len(N)) {
            s <- shifts[[i]]
            if (!is.numeric(s) || length(s) != 2) {
                stop(sprintf("shifts[[%d]] must be numeric(2)", i), call. = FALSE)
            }
            if (s[1] >= s[2]) {
                stop(sprintf("shifts[[%d]]: sshift (%g) must be < eshift (%g)", i, s[1], s[2]),
                    call. = FALSE
                )
            }
            sshifts[i] <- s[1]
            eshifts[i] <- s[2]
        }
    }

    # --- Validate kernel_bins and kernels ---
    if (!is.null(kernel_bins)) {
        if (!is.numeric(kernel_bins)) {
            stop("'kernel_bins' must be a numeric vector", call. = FALSE)
        }
        if (is.null(kernels)) {
            stop("'kernels' must be provided when 'kernel_bins' is specified", call. = FALSE)
        }
    }
    if (!is.null(kernels)) {
        if (!is.list(kernels)) {
            stop("'kernels' must be a list of numeric vectors", call. = FALSE)
        }
        B <- if (!is.null(kernel_bins)) length(kernel_bins) else NA
        kn_len <- length(kernels)
        if (kn_len != 1 && kn_len != N) {
            stop(sprintf("'kernels' must be length 1 or %d", N), call. = FALSE)
        }
        for (i in seq_along(kernels)) {
            if (!is.numeric(kernels[[i]])) {
                stop(sprintf("kernels[[%d]] must be numeric", i), call. = FALSE)
            }
            if (!is.na(B) && length(kernels[[i]]) != B) {
                stop(sprintf("kernels[[%d]] must have length %d (matching kernel_bins)", i, B),
                    call. = FALSE
                )
            }
        }
        # Recycle length-1 kernels
        if (kn_len == 1 && N > 1) {
            kernels <- rep(kernels, N)
        }
    }

    # --- Validate and recycle trans_family / trans_params ---
    if (!is.null(trans_family)) {
        if (length(trans_family) == 1) {
            trans_family <- rep(trans_family, N)
        }
        if (length(trans_family) != N) {
            stop(sprintf("'trans_family' must be length 1, %d, or NULL", N), call. = FALSE)
        }
        invalid_tf <- setdiff(unique(trans_family[!is.na(trans_family)]), "logist")
        if (length(invalid_tf) > 0) {
            stop(
                sprintf(
                    "'trans_family' must be 'logist' or NA (got: %s)",
                    paste0("'", invalid_tf, "'", collapse = ", ")
                ),
                call. = FALSE
            )
        }
    } else {
        trans_family <- rep(NA_character_, N)
    }

    if (!is.null(trans_params)) {
        if (length(trans_params) == 1) {
            trans_params <- rep(trans_params, N)
        }
        if (length(trans_params) != N) {
            stop(sprintf("'trans_params' must be length 1, %d, or NULL", N), call. = FALSE)
        }
    } else {
        trans_params <- vector("list", N)
    }

    # --- Validate max_cap / dis_from_cap ---
    if (is.null(max_cap)) max_cap <- rep(NA_real_, N)
    if (is.null(dis_from_cap)) dis_from_cap <- rep(NA_real_, N)
    if (length(max_cap) != N) {
        stop(sprintf("'max_cap' must be length %d or NULL", N), call. = FALSE)
    }
    if (length(dis_from_cap) != N) {
        stop(sprintf("'dis_from_cap' must be length %d or NULL", N), call. = FALSE)
    }
    # max_cap and dis_from_cap must be specified together per entry
    cap_mismatch <- xor(is.na(max_cap), is.na(dis_from_cap))
    if (any(cap_mismatch)) {
        bad <- which(cap_mismatch)[1]
        stop(
            sprintf(
                "'max_cap' and 'dis_from_cap' must both be specified or both NA for each entry (mismatch at entry %d)",
                bad
            ),
            call. = FALSE
        )
    }

    # --- Flatten transform params into parallel vectors ---
    trans_L <- vapply(seq_len(N), function(i) {
        p <- trans_params[[i]]
        if (is.null(p) || is.na(trans_family[i])) NA_real_ else p$L %||% 1
    }, numeric(1))
    trans_k <- vapply(seq_len(N), function(i) {
        p <- trans_params[[i]]
        if (is.null(p) || is.na(trans_family[i])) NA_real_ else p$k %||% 1
    }, numeric(1))
    trans_x0 <- vapply(seq_len(N), function(i) {
        p <- trans_params[[i]]
        if (is.null(p) || is.na(trans_family[i])) NA_real_ else p$x_0 %||% 0
    }, numeric(1))
    trans_pre <- vapply(seq_len(N), function(i) {
        p <- trans_params[[i]]
        if (is.null(p) || is.na(trans_family[i])) NA_real_ else p$pre_shift %||% 0
    }, numeric(1))
    trans_post <- vapply(seq_len(N), function(i) {
        p <- trans_params[[i]]
        if (is.null(p) || is.na(trans_family[i])) NA_real_ else p$post_shift %||% 0
    }, numeric(1))

    # --- Build params list for C++ ---
    # weights / bias / interaction_weights keep their array shape (with auto-derived
    # dimnames) when a selector is specified. C++ reads via REAL(), which gives the
    # column-major underlying buffer regardless of dim/dimnames attributes.
    params <- list(
        tracks = tracks,
        inner_func = inner_func,
        weights = weights,
        bias = bias,
        num_bins = as.integer(K),
        scale_factor = as.numeric(scale_factor),
        sshifts = as.numeric(sshifts),
        eshifts = as.numeric(eshifts),
        max_cap = as.numeric(max_cap),
        dis_from_cap = as.numeric(dis_from_cap),
        simple_cap = if (!is.null(simple_cap)) as.logical(simple_cap) else NULL,
        trans_family = trans_family,
        trans_L = as.numeric(trans_L),
        trans_k = as.numeric(trans_k),
        trans_x0 = as.numeric(trans_x0),
        trans_pre = as.numeric(trans_pre),
        trans_post = as.numeric(trans_post)
    )

    # Add selector params if present
    if (!is.null(selector_tracks)) {
        params$selector_tracks <- as.character(selector_tracks)
        params$selector_breaks <- lapply(selector_breaks, as.numeric)
    }

    # Add kernel params if present
    if (!is.null(kernel_bins)) {
        params$kernel_bins <- as.numeric(kernel_bins)
    }
    if (!is.null(kernels)) {
        params$kernels <- kernels
    }

    # --- Validate and flatten interactions ---
    if (!is.null(interactions)) {
        if (!is.list(interactions)) {
            stop("'interactions' must be a list of integer(2) pairs", call. = FALSE)
        }
        M <- length(interactions)
        if (is.null(selector_tracks)) {
            if (!is.numeric(interaction_weights) || !is.null(dim(interaction_weights)) ||
                length(interaction_weights) != M) {
                stop(sprintf(
                    "'interaction_weights' must be a numeric vector of length %d when no selector is specified",
                    M
                ), call. = FALSE)
            }
            interaction_weights <- as.numeric(interaction_weights)
        } else {
            interaction_weights <- .glm_pred_normalize_strata_array(
                interaction_weights,
                K_per = K_per,
                selector_tracks = selector_tracks, selector_breaks = selector_breaks,
                leading_dim = M, leading_name = "interaction",
                param = "interaction_weights"
            )
        }

        inter_i <- vapply(interactions, `[`, integer(1), 1L)
        inter_j <- vapply(interactions, `[`, integer(1), 2L)

        if (any(inter_i < 1 | inter_i > N | inter_j < 1 | inter_j > N)) {
            stop(sprintf("Interaction indices must be in [1, %d]", N), call. = FALSE)
        }

        params$inter_i <- as.integer(inter_i)
        params$inter_j <- as.integer(inter_j)
        params$inter_weights <- interaction_weights

        # Process interaction transforms
        if (is.null(interaction_trans_family)) {
            interaction_trans_family <- rep(NA_character_, M)
        }
        if (length(interaction_trans_family) == 1) {
            interaction_trans_family <- rep(interaction_trans_family, M)
        }

        if (is.null(inter_trans_params)) {
            inter_trans_params <- vector("list", M)
        }
        if (length(inter_trans_params) == 1) {
            inter_trans_params <- rep(inter_trans_params, M)
        }

        params$inter_trans_family <- interaction_trans_family
        params$inter_trans_L <- vapply(seq_len(M), function(m) {
            p <- inter_trans_params[[m]]
            if (is.null(p) || is.na(interaction_trans_family[m])) NA_real_ else p$L %||% 1
        }, numeric(1))
        params$inter_trans_k <- vapply(seq_len(M), function(m) {
            p <- inter_trans_params[[m]]
            if (is.null(p) || is.na(interaction_trans_family[m])) NA_real_ else p$k %||% 1
        }, numeric(1))
        params$inter_trans_x0 <- vapply(seq_len(M), function(m) {
            p <- inter_trans_params[[m]]
            if (is.null(p) || is.na(interaction_trans_family[m])) NA_real_ else p$x_0 %||% 0
        }, numeric(1))
        params$inter_trans_pre <- vapply(seq_len(M), function(m) {
            p <- inter_trans_params[[m]]
            if (is.null(p) || is.na(interaction_trans_family[m])) NA_real_ else p$pre_shift %||% 0
        }, numeric(1))
        params$inter_trans_post <- vapply(seq_len(M), function(m) {
            p <- inter_trans_params[[m]]
            if (is.null(p) || is.na(interaction_trans_family[m])) NA_real_ else p$post_shift %||% 0
        }, numeric(1))
    }

    # --- Register as vtrack ---
    var <- list(
        func = "glm.predict",
        params = params
    )

    .gvtrack.set(name, var)

    invisible(name)
}

#' Remove a GLM predictor virtual track
#'
#' @param name character(1) Virtual track name
#' @return None.
#' @export
glm_pred.rm <- function(name) {
    gvtrack.rm(name)
}

#' List GLM predictor virtual tracks
#'
#' @return Character vector of virtual track names that use glm.predict.
#' @export
glm_pred.ls <- function() {
    all_vt <- gvtrack.ls()
    if (length(all_vt) == 0) {
        return(character(0))
    }
    is_glm <- vapply(all_vt, function(vt) {
        info <- gvtrack.info(vt)
        identical(info$func, "glm.predict")
    }, logical(1))
    all_vt[is_glm]
}

#' Get info for a GLM predictor virtual track
#'
#' Returns the virtual track definition. When the track uses one or more
#' selectors, \code{weights}, \code{bias}, and \code{interaction_weights}
#' are returned as labeled arrays with \code{dim = c(N, K_1, ..., K_M)}
#' (or \code{c(K_1, ..., K_M)} for \code{bias}), with axis names taken
#' from \code{selector_tracks} and axis labels derived from
#' \code{selector_breaks} via \code{cut()} semantics.
#'
#' @param name character(1) Virtual track name
#' @return List with virtual track definition.
#' @export
glm_pred.info <- function(name) {
    gvtrack.info(name)
}
