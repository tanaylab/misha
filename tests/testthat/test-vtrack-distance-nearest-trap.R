create_isolated_test_db()

# Regression tests for the greedy local-minimum scan bug in the distance
# vtrack family (distance / distance.edge / distance.center). With overlapping
# or nested source intervals the old sequential scan stopped at the first local
# minimum of the distance sequence and returned a non-nearest interval - in the
# extreme, a nonzero distance for a query that overlaps the source.
#
# Trap geometry (single chromosome):
#   C  = [800, 900]    upstream of query
#   C2 = [850, 870]    nested inside C  (sources overlap each other)
#   D  = [950, 1500]   OVERLAPS the query  -> true nearest, distance 0
#   E  = [1250, 1300]  downstream, placed to also trap the end-sorted scan
# query = [1000, 1100], center 1050, sits inside D.
trap_src <- function() {
    gintervals(1, c(800, 850, 950, 1250), c(900, 870, 1500, 1300))
}
trap_query <- gintervals(1, 1000, 1100)

test_that("distance.edge returns 0 for a query overlapping the source (trap geometry)", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    gvtrack.create("d_edge_trap", trap_src(), "distance.edge")
    res <- gextract("d_edge_trap", trap_query, iterator = trap_query)

    # Query overlaps D=[950,1500] -> distance must be 0, matching gintervals.neighbors.
    nb <- gintervals.neighbors(trap_query, trap_src())
    expect_equal(res$d_edge_trap, 0)
    expect_equal(res$d_edge_trap, nb$dist[1])
})

test_that("distance returns 0 when the query center sits inside an overlapping source (trap geometry)", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    gvtrack.create("d_plain_trap", trap_src(), "distance")
    res <- gextract("d_plain_trap", trap_query, iterator = trap_query)

    # center 1050 is inside D -> distance 0 (dist2coord returns 0 inside an interval).
    expect_equal(res$d_plain_trap, 0)
})

test_that("distance.edge matches gintervals.neighbors on a dense overlapping source (stress)", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    set.seed(42)
    n_src <- 2000
    s <- sample.int(1e6, n_src)
    w <- sample(c(10:50, 500:2000), n_src, replace = TRUE) # mix of narrow and very wide
    src <- gintervals(1, s, s + w)
    gvtrack.create("d_edge_stress", src, "distance.edge")

    n_q <- 3000
    qs <- sample.int(1e6, n_q)
    query <- gintervals(1, qs, qs + 30L)
    query <- gintervals.canonic(query)

    res <- gextract("d_edge_stress", query, iterator = query)
    res <- res[order(res$intervalID), ]
    nb <- gintervals.neighbors(query, src)
    key <- match(paste(res$chrom, res$start, res$end), paste(nb$chrom, nb$start, nb$end))

    expect_equal(res$d_edge_stress, nb$dist[key])
})

test_that("distance.center accepts overlapping source intervals (no longer errors)", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Previously this threw "intervals are overlapping and hence incompatible".
    expect_no_error(gvtrack.create("dc_overlap_src", trap_src(), "distance.center"))
})

test_that("distance.center returns distance to the containing interval's center (overlapping source)", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    gvtrack.create("dc_single", trap_src(), "distance.center")
    res <- gextract("dc_single", trap_query, iterator = trap_query)

    # center 1050 is inside D=[950,1500] only; D's center is 1225 -> |1050-1225| = 175.
    expect_equal(res$dc_single, 175)
})

test_that("distance.center picks the nearest center among overlapping containing intervals", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # C=[800,900] center 850, C2=[850,870] center 860 (C2 nested in C).
    src <- gintervals(1, c(800, 850), c(900, 870))
    gvtrack.create("dc_multi", src, "distance.center")

    # query center 860 is inside both C and C2; nearest center is C2 (860) -> 0.
    query <- gintervals(1, 855, 865)
    res <- gextract("dc_multi", query, iterator = query)
    expect_equal(res$dc_multi, 0)
})
