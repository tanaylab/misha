#' Choose the multitasking strategy for gextract
#'
#' Read \code{getOption("gmultitasking.strategy", "auto")} and resolve "auto"
#' to either "tiles" (the historical misha multitask: each kid handles a tile
#' range across all tracks) or "tracks" (R-side mclapply: each worker handles
#' a track subset across all tiles). Track-parallel is a major win on
#' cold-NFS many-track gextract workloads (~5× measured on Tamar's 30-track
#' × 2.19M-bin bench) because each worker only mmap-faults its own files,
#' avoiding the working-set thrashing that 24 kids × 500 tracks creates on
#' the kernel page cache.
#'
#' Falls back to "tiles" when track-parallel doesn't apply (single track,
#' file/bigset output, 2D band iteration, or workload too small to amortize
#' fork overhead).
#'
#' @param tracks character vector of track expressions
#' @param intervals data.frame of intervals (or NULL for big intervals sets)
#' @param file optional output file path (NULL = data.frame return)
#' @param intervals.set.out optional bigset output name
#' @param band optional 2D band parameter
#' @return one of "tiles" or "tracks"
#' @keywords internal
.gmultitasking_strategy <- function(tracks, intervals, file = NULL,
                                    intervals.set.out = NULL, band = NULL) {
    strategy <- getOption("gmultitasking.strategy", "auto")

    # Explicit override always wins.
    if (identical(strategy, "tiles") || identical(strategy, "tracks")) {
        return(strategy)
    }
    if (!identical(strategy, "auto")) {
        warning(
            sprintf("Unknown gmultitasking.strategy '%s'; falling back to 'auto'.",
                    as.character(strategy)),
            call. = FALSE
        )
    }

    # Track-parallel disqualifiers.
    if (length(tracks) < 2) return("tiles")
    if (!is.null(file) || !is.null(intervals.set.out)) return("tiles")
    if (!is.null(band)) return("tiles")

    # Number of intervals — works for plain data.frames and big sets we can
    # measure without forcing a load.
    n_intervals <- tryCatch(
        if (is.data.frame(intervals)) nrow(intervals) else NA_integer_,
        error = function(e) NA_integer_
    )
    if (is.na(n_intervals)) {
        # Couldn't determine size — fall back to tiles (existing behaviour).
        return("tiles")
    }

    # Workload threshold. Below ~1e6 (track, interval) pairs the per-worker
    # fork + result-merge cost outweighs the I/O parallelism benefit.
    if ((as.numeric(length(tracks)) * as.numeric(n_intervals)) < 1e6) {
        return("tiles")
    }

    "tracks"
}

#' Track-parallel gextract via mclapply
#'
#' Splits \code{tracks} into chunks across at most \code{getOption("gmax.processes")}
#' worker processes. Each worker runs gextract on its track subset with
#' \code{gmultitasking=FALSE} so misha's tile-parallel multitask doesn't nest.
#' Results are merged column-wise: interval/intervalID columns come from the
#' first worker, value columns are cbind'd from each.
#'
#' @keywords internal
.gextract_track_parallel <- function(intervals, tracks, colnames, iterator,
                                     band, file, intervals.set.out, envir) {
    if (.Platform$OS.type != "unix") {
        # mclapply forks; on Windows fall back to the regular tile-parallel path.
        return(.gcall("gextract_multitask", intervals, tracks, colnames,
                      iterator, band, file, intervals.set.out, envir))
    }

    n_workers <- as.integer(.ggetOption("gmax.processes"))
    if (is.na(n_workers) || n_workers < 1) n_workers <- 1L
    n_workers <- min(n_workers, length(tracks))

    # Split track expressions into n_workers contiguous chunks. Round-robin
    # assignment helps balance per-track cost when tracks have similar sizes.
    chunks <- split(tracks, rep(seq_len(n_workers), length.out = length(tracks)))
    if (!is.null(colnames)) {
        chunks_names <- split(colnames,
                              rep(seq_len(n_workers), length.out = length(colnames)))
    } else {
        chunks_names <- replicate(n_workers, NULL, simplify = FALSE)
    }

    # Disable the C++ multitask inside each worker — we already forked here.
    worker <- function(idx) {
        old_mt <- options(gmultitasking = FALSE)
        on.exit(options(old_mt), add = TRUE)
        .gcall("C_gextract", intervals, chunks[[idx]], chunks_names[[idx]],
               iterator, band, NULL, NULL, envir)
    }

    results <- parallel::mclapply(seq_len(n_workers), worker,
                                  mc.cores = n_workers, mc.preschedule = FALSE)

    # Surface any worker errors.
    errs <- vapply(results, inherits, logical(1), what = "try-error")
    if (any(errs)) {
        first_err <- attr(results[[which(errs)[1]]], "condition")
        stop(if (!is.null(first_err)) conditionMessage(first_err)
             else "track-parallel gextract worker failed", call. = FALSE)
    }
    nullish <- vapply(results, is.null, logical(1))
    if (all(nullish)) return(NULL)

    # Pick the first non-null result as the row scaffold; cbind value columns
    # from every worker.
    first_idx <- which(!nullish)[1]
    base <- results[[first_idx]]
    n_track_cols <- length(chunks[[first_idx]])

    # Identify which columns are intervals/intervalID vs value columns. The
    # final column is intervalID for 1D output (per gextract.cpp); intervals
    # cols are the leading 3 (chrom/start/end) or 6 (2D).
    n_total <- ncol(base)
    iv_n <- if ("chrom2" %in% names(base)) 6L else 3L
    has_id <- isTRUE(names(base)[n_total] == "intervalID")
    val_cols_first <- seq.int(iv_n + 1L, n_total - as.integer(has_id))

    out_value_cols <- list()
    out_value_names <- character(0)
    for (i in seq_len(n_workers)) {
        if (nullish[i]) next
        n_chunk <- length(chunks[[i]])
        # value columns occupy positions [iv_n+1 .. iv_n+n_chunk]
        idx <- seq.int(iv_n + 1L, iv_n + n_chunk)
        out_value_cols[[length(out_value_cols) + 1L]] <- results[[i]][, idx, drop = FALSE]
        out_value_names <- c(out_value_names, names(results[[i]])[idx])
    }
    value_df <- do.call(cbind, out_value_cols)
    names(value_df) <- out_value_names

    iv_df <- base[, seq_len(iv_n), drop = FALSE]
    out <- cbind(iv_df, value_df)
    if (has_id) {
        out[["intervalID"]] <- base[["intervalID"]]
    }
    rownames(out) <- NULL
    class(out) <- "data.frame"
    out
}
