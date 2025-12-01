# Core track helper functions

# tries to locate track name in a path, return the track name
# examples:
#   .gfindtrackinpath("aaa/bbb/ccc.track/ddd/eee")       returns "aaa.bbb.ccc"
#   .gfindtrackinpath("aaa/bbb/ccc.track/ddd.track/eee") returns "aaa.bbb.ccc"
#   .gfindtrackinpath("aaa/bbb/ccc/ddd/eee")             returns NULL
.gfindtrackinpath <- function(path) {
    dirs <- unlist(strsplit(path, split = "/"))
    r <- grep("\\.track$", unlist(dirs))
    if (length(r) > 0) {
        idx <- r[1]
        dirs[idx] <- paste(substr(dirs[idx], 0, nchar(dirs[idx]) - nchar(".tracks") + 1))
        return(paste(dirs[1:idx], collapse = "."))
    }
    NULL
}


.gtrack.prepare.pvals <- function(track) {
    .gcheckroot()

    trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())

    if (is.na(match(trackstr, get("GTRACKS", envir = .misha)))) {
        stop(sprintf("Track %s does not exist", trackstr), call. = FALSE)
    }

    quantile.edge.data.size <- .ggetOption("gquantile.edge.data.size", 100000)
    middle.size <- .ggetOption("gpv.middle.size", 0.98)
    edge.size <- (1 - middle.size) / 2
    middle.precision <- .ggetOption("gpv.middle.precision", 10^(-4))
    edge.precision <- .ggetOption("gpv.edge.precision", 10^(-9))

    # In the middle section the percentiles increase by a constant step = middle.precision.
    # At the edges the step is exponential. For instance, step i at the edge close to 0 follows the following function:
    #   k+i
    #  b     , where b and k are unknown. We also don't know n - the number of total steps required to cover the edge.
    #
    # b, k, n be calculated from the following:
    #
    # 1. b^(k+1) - b^k = edge.precision
    # 2. b^(k+n) - b^(k+n-1) = middle.precision
    # 3. b^(k+n) - b^k = edge.size

    b <- (edge.precision - edge.size) / (middle.precision - edge.size)
    k <- log(edge.precision / (b - 1), b)
    n <- ceiling(log(edge.size + b^k, b) - k)

    percentiles <- 0
    percentiles <- c(percentiles, b^(k + (0:n)) - b^k)
    percentiles <- c(percentiles, 1 - percentiles)
    num.middle.steps <- middle.size / middle.precision
    percentiles <- c(percentiles, edge.size + (0:num.middle.steps) * middle.precision)
    percentiles <- sort(percentiles)

    selected.percentiles <- NULL

    multitasking <- .ggetOption("gmultitasking")
    on.exit(options(gmultitasking = multitasking))
    tryCatch(
        {
            suppressWarnings({ # disable warnings since gquantiles is going to warn about random sampling
                options(gmultitasking = FALSE)
                quantiles <- do.call(gquantiles, list(substitute(track), percentiles = c(0, percentiles)), envir = parent.frame())
                names(quantiles) <- NULL
                minval <- quantiles[1]
                maxval <- quantiles[length(quantiles)]
                quantiles <- quantiles[2:length(quantiles)]

                # for each group of quantiles with identical value choose the maximal one
                selected.percentiles <- sapply(split(percentiles, quantiles), max)
                names(selected.percentiles) <- NULL

                # if all percentiles are equal create an artificial table
                if (length(selected.percentiles) == 1) {
                    selected.percentiles <- c(1, 1)
                    attr(selected.percentiles, "breaks") <- c(minval, maxval + 1)
                } else {
                    indices <- match(selected.percentiles, percentiles)
                    selected.quantiles <- quantiles[indices]
                    attr(selected.percentiles, "breaks") <- selected.quantiles
                }

                attr(selected.percentiles, "minval") <- minval
                attr(selected.percentiles, "maxval") <- maxval
            })
        },
        finally = {
            options(gmultitasking = multitasking)
        }
    )

    # save the percentiles
    .gtrack.var.set(trackstr, "pv.percentiles", selected.percentiles)

    retv <- 0
}
