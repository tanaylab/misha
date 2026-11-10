# The 2D value summarizers (avg, min, max, weighted.sum, and a 2D track expression) skip NaN values,
# as the 1D summarizers do; area stays the area covered by the track rectangles, NaN ones included.
# min and max used to depend on where a NaN cell sat - std::min / std::max let a NaN take over until
# a later cell replaced it - and avg and weighted.sum were NaN. The statistics the quad-tree stores
# per node are written with the NaN in them, so a query that covers such a node got a NaN or a wrong
# min / max.
#
# The track writer is unchanged, so the tracks created here carry the same node statistics, NaN
# included, as tracks written by earlier versions: reading them right is reading old tracks right.

create_isolated_test_db()

nan2d_funcs <- c("avg", "min", "max", "weighted.sum", "area")
nan2d_value_funcs <- c("expr", "avg", "min", "max", "weighted.sum")

nan2d_extract <- function(track, q, band = NULL) {
    vtracks <- paste0("v_", sub(".", "_", nan2d_funcs, fixed = TRUE))
    for (i in seq_along(vtracks)) gvtrack.create(vtracks[i], track, nan2d_funcs[i])
    withr::defer(for (v in vtracks) gvtrack.rm(v))
    # one query at a time: the queries overlap, which a 2D iterator does not allow
    res <- do.call(rbind, lapply(seq_len(nrow(q)), function(i) {
        gextract(c(track, vtracks), q[i, ], iterator = q[i, ], band = band, colnames = c("expr", nan2d_funcs))
    }))
    as.list(res[, c("expr", nan2d_funcs)])
}

# the value summarizers over the cells that have a value, and the area of all the cells
nan2d_expected <- function(cells, vals, q) {
    lapply(seq_len(nrow(q)), function(i) {
        a <- pmax(0, pmin(cells$end1, q$end1[i]) - pmax(cells$start1, q$start1[i])) *
            pmax(0, pmin(cells$end2, q$end2[i]) - pmax(cells$start2, q$start2[i]))
        keep <- a > 0 & !is.nan(vals)
        if (!any(keep)) {
            return(c(expr = NaN, avg = NaN, min = NaN, max = NaN, weighted.sum = NaN, area = sum(a)))
        }
        ws <- sum(vals[keep] * a[keep])
        c(
            expr = ws / sum(a[keep]), avg = ws / sum(a[keep]), min = min(vals[keep]), max = max(vals[keep]),
            weighted.sum = ws, area = sum(a)
        )
    })
}

test_that("a NaN cell is skipped wherever it sits", {
    withr::defer(if (gtrack.exists("temp.nan2d")) gtrack.rm("temp.nan2d", force = TRUE))
    rects <- gintervals.2d(1, c(0, 1000), c(100, 1100), 1, c(0, 1000), c(100, 1100))
    q <- gintervals.2d(1, 0, 2000, 1, 0, 2000)

    for (vals in list(c(NaN, 5), c(5, NaN))) {
        gtrack.2d.create("temp.nan2d", "", rects, vals)
        got <- nan2d_extract("temp.nan2d", q)
        expect_equal(unlist(got), c(expr = 5, avg = 5, min = 5, max = 5, weighted.sum = 5e4, area = 2e4))
        gtrack.rm("temp.nan2d", force = TRUE)
    }

    gtrack.2d.create("temp.nan2d", "", rects, c(NaN, NaN))
    got <- nan2d_extract("temp.nan2d", q)
    expect_equal(unlist(got), c(expr = NaN, avg = NaN, min = NaN, max = NaN, weighted.sum = NaN, area = 2e4))
})

test_that("node statistics written with a NaN in them are not used", {
    kinds <- c("rects", "points")
    tracks <- c(outer(paste0("temp.nan2d_", kinds), c("", "_ref", "_filled"), paste0))
    withr::defer(for (tr in tracks) if (gtrack.exists(tr)) gtrack.rm(tr, force = TRUE))

    # 900 cells on a grid: far more than the 20 objects of a quad-tree leaf, so the queries below
    # cover whole nodes and read their stored statistics
    set.seed(17)
    grid <- expand.grid(i = 0:29, j = 0:29)
    cells <- gintervals.2d(1, grid$i * 200, grid$i * 200 + 100, 1, grid$j * 200, grid$j * 200 + 100)
    vals <- sample(1:50, nrow(cells), replace = TRUE)
    vals[sample(nrow(cells), 120)] <- NaN
    vals[nrow(cells)] <- NaN # the last cell written: the last value every node above it saw
    has_val <- !is.nan(vals)

    # each track three times: as is, without its NaN cells (_ref), and with 1 for NaN (_filled)
    all_cells <- rep(TRUE, nrow(cells))
    versions <- list(
        list(suffix = "", keep = all_cells, vals = vals),
        list(suffix = "_ref", keep = has_val, vals = vals),
        list(suffix = "_filled", keep = all_cells, vals = ifelse(has_val, vals, 1))
    )
    for (v in versions) {
        suffix <- v$suffix
        gtrack.2d.create(paste0("temp.nan2d_rects", suffix), "", cells[v$keep, ], v$vals[v$keep])
        src <- withr::local_tempfile(fileext = ".txt")
        utils::write.table(
            data.frame(
                chrom1 = "chr1", start1 = cells$start1[v$keep], end1 = cells$start1[v$keep] + 1,
                chrom2 = "chr1", start2 = cells$start2[v$keep], end2 = cells$start2[v$keep] + 1, value = v$vals[v$keep]
            ),
            src,
            sep = "\t", row.names = FALSE, quote = FALSE, na = "NaN"
        )
        gtrack.2d.import(paste0("temp.nan2d_points", suffix), "", src)
    }
    expect_equal(gtrack.info("temp.nan2d_points")$type, "points")

    q <- rbind(
        gintervals.2d(1, 0, 6000, 1, 0, 6000),
        gintervals.2d(1, 0, 3000, 1, 3000, 6000),
        gintervals.2d(1, c(150, 1230, 2970), c(5310, 4400, 3020), 1, c(40, 610, 0), c(4890, 5990, 2900))
    )
    expected <- nan2d_expected(cells, vals, q)
    got <- nan2d_extract("temp.nan2d_rects", q)
    for (f in names(got)) expect_equal(got[[f]], sapply(expected, `[[`, f), tolerance = 1e-6, info = f)

    # also with a band and for points: the values are those of the track without its NaN cells, the
    # area that of the track with a value in every cell
    for (band in list(NULL, c(-2500, 0), c(-1e5, -1300), c(700, 1e5))) {
        for (kind in kinds) {
            tr <- paste0("temp.nan2d_", kind)
            got <- nan2d_extract(tr, q, band)
            ref <- nan2d_extract(paste0(tr, "_ref"), q, band)
            filled <- nan2d_extract(paste0(tr, "_filled"), q, band)
            info <- paste(kind, toString(band))
            for (f in nan2d_value_funcs) expect_equal(got[[f]], ref[[f]], tolerance = 1e-6, info = paste(f, info))
            expect_equal(got$area, filled$area, info = info)
        }
    }

    # a 2D track used as intervals still covers its NaN cells
    expect_equal(misha:::.gintervals.big.meta("temp.nan2d_rects")$stats$surface, 100 * 100 * nrow(cells))
})
