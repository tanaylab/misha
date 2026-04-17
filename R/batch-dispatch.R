# Internal helpers for the batched fast-path dispatcher used by gsummary
# (Phase 3), gscreen (Phase 4), and gquantiles (Phase 6).
#
# Contract: detect_fast_path(exprs, iterator, intervals, band) returns
# either NULL (not fast-path eligible, caller should use slow path) or a
# list with fields:
#   $tracks   character vector of underlying track names
#   $func     "lse" | "avg" | "sum" | "max" | "min"
#   $sshift   integer
#   $eshift   integer
#
# Preconditions for eligibility:
#   1. Each expr is a bare track name OR a vtrack wrapping a single
#      source track with func in {avg, sum, max, min, lse}.
#   2. All exprs share the same (func, sshift, eshift) tuple.
#   3. iterator is a fixed integer step (not a track-based iterator).
#   4. band is NULL.
#   5. intervals is either NULL, ALLGENOME, or a 1D intervals data.frame.

.describe_single_expr <- function(e) {
    # Bare track?
    if (gtrack.exists(e)) {
        info <- tryCatch(gtrack.info(e), error = function(err) NULL)
        if (is.null(info)) return(NULL)
        if (!identical(info$type, "dense") && !identical(info$type, "sparse"))
            return(NULL)
        # Bare track → per-bin scan semantics. We return sshift=0 and
        # eshift=bin_size for dense tracks so the window covers exactly one
        # bin at each iterator position (and func=avg collapses to the bin
        # value). For sparse tracks there's no native bin size; use a
        # placeholder 1 and rely on the caller's iterator. The caller may
        # override by wrapping in a vtrack.
        bsz <- info$bin.size
        if (is.null(bsz) || !is.numeric(bsz)) bsz <- 1L
        return(list(track = e, func = "avg", sshift = 0L,
                    eshift = as.integer(bsz)))
    }
    # Virtual track?
    if (exists("GVTRACKS", envir = misha:::.misha)) {
        gwd <- get("GWD", envir = misha:::.misha)
        vts <- get("GVTRACKS", envir = misha:::.misha)[[gwd]]
        if (!is.null(vts) && e %in% names(vts)) {
            v <- vts[[e]]
            if (!is.character(v$src) || length(v$src) != 1) return(NULL)
            if (!gtrack.exists(v$src)) return(NULL)
            if (!is.character(v$func) || length(v$func) != 1) return(NULL)
            if (!(v$func %in% c("avg", "sum", "max", "min", "lse")))
                return(NULL)
            if (is.null(v$itr) || !identical(v$itr$type, "1d")) return(NULL)
            sshift <- as.integer(v$itr$sshift)
            eshift <- as.integer(v$itr$eshift)
            if (is.na(sshift) || is.na(eshift)) return(NULL)
            return(list(track = v$src, func = v$func,
                        sshift = sshift, eshift = eshift))
        }
    }
    NULL
}

.detect_fast_path <- function(exprs, iterator, intervals, band) {
    if (!is.character(exprs) || length(exprs) == 0) return(NULL)
    if (!is.null(band)) return(NULL)

    infos <- lapply(exprs, .describe_single_expr)
    if (any(vapply(infos, is.null, logical(1)))) return(NULL)

    # iterator is optional when every expression is a bare track with a
    # known bin size — default to the common bin size (all must match).
    # Otherwise iterator is required.
    if (is.null(iterator)) {
        eshifts <- vapply(infos, `[[`, integer(1), "eshift")
        sshifts <- vapply(infos, `[[`, integer(1), "sshift")
        is_bare <- sshifts == 0L
        if (!all(is_bare)) return(NULL)
        if (length(unique(eshifts)) != 1) return(NULL)
        it_int <- eshifts[1]
    } else {
        if (!is.numeric(iterator) || length(iterator) != 1) return(NULL)
        it_int <- as.integer(iterator)
        if (is.na(it_int) || it_int <= 0) return(NULL)
    }

    funcs  <- vapply(infos, `[[`, character(1), "func")
    sshift <- vapply(infos, `[[`, integer(1), "sshift")
    eshift <- vapply(infos, `[[`, integer(1), "eshift")

    if (length(unique(funcs)) != 1) return(NULL)
    if (length(unique(sshift)) != 1) return(NULL)
    if (length(unique(eshift)) != 1) return(NULL)

    # Intervals: accept 1D data.frame. ALLGENOME is a list of
    # (1D_df, 2D_df); unwrap to the 1D part. Reject bigset handles
    # (character strings) and anything 2D.
    if (!is.null(intervals)) {
        iv <- intervals
        if (is.list(iv) && !is.data.frame(iv) && length(iv) == 2 &&
            is.data.frame(iv[[1]])) {
            iv <- iv[[1]]
        }
        if (!is.data.frame(iv)) return(NULL)
        if (!all(c("chrom", "start", "end") %in% colnames(iv)))
            return(NULL)
        if ("chrom1" %in% colnames(iv)) return(NULL)
    }

    list(
        tracks = vapply(infos, `[[`, character(1), "track"),
        func = funcs[1],
        sshift = sshift[1],
        eshift = eshift[1],
        iterator = it_int
    )
}


# Screen-specific parser. Each expression must be a single comparison:
#   "<lhs> <op> <const>"  where op in <, <=, ==, >=, >
# and <lhs> must satisfy .describe_single_expr (bare track or simple
# vtrack). Returns either NULL or a list with tracks/func/sshift/eshift/
# iterator/ops/thresholds fields.
.detect_screen_fast_path <- function(exprs, iterator, intervals, band) {
    if (!is.character(exprs) || length(exprs) == 0) return(NULL)
    if (!is.null(band)) return(NULL)
    # Parse "<lhs> <op> <const>" from each expr.
    rx <- "^\\s*(.+?)\\s*(<=|>=|==|<|>)\\s*([-+0-9.eE]+)\\s*$"
    m <- regmatches(exprs, regexec(rx, exprs))
    if (any(vapply(m, function(x) length(x) != 4, logical(1)))) return(NULL)
    lhs <- vapply(m, `[`, character(1), 2)
    ops <- vapply(m, `[`, character(1), 3)
    thr <- suppressWarnings(as.numeric(vapply(m, `[`, character(1), 4)))
    if (any(is.na(thr))) return(NULL)

    infos <- lapply(lhs, .describe_single_expr)
    if (any(vapply(infos, is.null, logical(1)))) return(NULL)

    funcs <- vapply(infos, `[[`, character(1), "func")
    sshifts <- vapply(infos, `[[`, integer(1), "sshift")
    eshifts <- vapply(infos, `[[`, integer(1), "eshift")

    if (is.null(iterator)) {
        is_bare <- sshifts == 0L
        if (!all(is_bare)) return(NULL)
        if (length(unique(eshifts)) != 1) return(NULL)
        it_int <- eshifts[1]
    } else {
        if (!is.numeric(iterator) || length(iterator) != 1) return(NULL)
        it_int <- as.integer(iterator)
        if (is.na(it_int) || it_int <= 0) return(NULL)
    }

    if (length(unique(funcs)) != 1) return(NULL)
    if (length(unique(sshifts)) != 1) return(NULL)
    if (length(unique(eshifts)) != 1) return(NULL)

    if (!is.null(intervals)) {
        iv <- intervals
        if (is.list(iv) && !is.data.frame(iv) && length(iv) == 2 &&
            is.data.frame(iv[[1]])) {
            iv <- iv[[1]]
        }
        if (!is.data.frame(iv)) return(NULL)
        if (!all(c("chrom", "start", "end") %in% colnames(iv))) return(NULL)
        if ("chrom1" %in% colnames(iv)) return(NULL)
    }

    list(
        tracks = vapply(infos, `[[`, character(1), "track"),
        func = funcs[1],
        sshift = sshifts[1],
        eshift = eshifts[1],
        iterator = it_int,
        ops = ops,
        thresholds = thr
    )
}

# Map comparison operator strings to the C-side CmpOp enum integers
# (see ThresholdScreen::CmpOp in BatchScreen.cpp).
.screen_op_to_int <- function(op_chr) {
    tbl <- c("<" = 0L, "<=" = 1L, "==" = 2L, ">=" = 3L, ">" = 4L)
    out <- tbl[op_chr]
    if (any(is.na(out))) stop("unknown comparison operator", call. = FALSE)
    as.integer(out)
}

# Emit a one-time per-session informational message explaining which path
# was taken (or why the fast path was declined). Suppressed via
# options(misha.quiet_dispatch = TRUE).
.fast_dispatch_msg <- function(fn, reason) {
    if (isTRUE(getOption("misha.quiet_dispatch"))) return(invisible())
    key <- paste0("misha.dispatch_msg.", fn)
    if (isTRUE(getOption(key))) return(invisible())
    packageStartupMessage(sprintf("[%s] %s", fn, reason))
    args <- setNames(list(TRUE), key)
    do.call(options, args)
}
