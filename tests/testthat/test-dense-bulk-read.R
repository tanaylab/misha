# A dense track read at its own bin size, with no virtual track in the expressions, is read a run of
# bins at a time (TrackExprScanner::fill_bins_bulk) instead of through read_interval() per iterator
# interval. Each call below is run both ways - MISHA_DISABLE_DENSE_BULK_READ=1 forces the
# per-interval path - and the results must be bit-identical: -0 is read as 0 and a NaN with a payload
# or a sign as R's NaN, as read_interval() reads them. The calls that may not use the bulk path
# (other bin sizes, virtual tracks, other variables) would give other values if they did.

# A database of its own: chromosome ends off the 50 bp grid and a chromosome shorter than a bin,
# per-chromosome and indexed copies of the tracks, NaN runs, a chromosome with no values, and a
# track whose files hold -0, -NaN, a NaN with a payload and infinities.
local_bulk_db <- function(env = parent.frame()) {
    local_db_state(env)
    db <- withr::local_tempdir(.local_envir = env)
    sizes <- data.frame(chrom = c("chr1", "chr2", "chr3", "chrX"), size = c(100007, 55555, 20000, 3))
    create_test_db(db, chrom_sizes = sizes)
    suppressMessages(gsetroot(db))
    withr::local_options(gmax.data.size = 1e9, .local_envir = env)

    set.seed(17)
    all <- gintervals.all()
    create <- function(track, binsize, nan_chrom = NULL) {
        iv <- do.call(rbind, lapply(seq_len(nrow(all)), function(i) {
            s <- seq(0, all$end[i] - 1, by = binsize)
            data.frame(chrom = all$chrom[i], start = s, end = pmin(s + binsize, all$end[i]))
        }))
        v <- round(rnorm(nrow(iv)), 3)
        for (st in sample.int(nrow(iv), 20)) v[st:min(nrow(iv), st + sample.int(100, 1))] <- NaN
        v[c(max(which(iv$chrom == "chr1")), min(which(iv$chrom == "chr2")))] <- NaN
        if (!is.null(nan_chrom)) v[iv$chrom == nan_chrom] <- NaN
        gtrack.create_dense(track, "bulk", iv, v, binsize, NaN)
    }
    create("d50a", 50)
    create("d50b", 50, nan_chrom = "chr3")
    create("d20", 20)
    create("d1", 1)

    gtrack.copy("d50a", "d50z")
    write_bins <- function(chrom, bins, bytes) {
        con <- file(file.path(db, "tracks", "d50z.track", chrom), "r+b")
        on.exit(close(con))
        for (b in bins) {
            seek(con, 4 + 4 * b, rw = "write")
            writeBin(as.raw(bytes), con)
        }
    }
    write_bins("chr1", c(0, 21, 22, 500, 2000), c(0, 0, 0, 0x80)) # -0
    write_bins("chr1", c(23, 501), c(0, 0, 0xc0, 0xff)) # -NaN
    write_bins("chr1", c(24, 502), c(1, 0, 0xa0, 0x7f)) # NaN with a payload
    write_bins("chr1", 25, c(0, 0, 0x80, 0x7f)) # inf
    write_bins("chr2", c(0, 1, 1111), c(0, 0, 0, 0x80))
    write_bins("chr2", 2, c(0, 0, 0x80, 0xff)) # -inf
    write_bins("chrX", 0, c(0, 0, 0, 0x80))
    gdb.reload()

    for (track in c("d50a", "d50b", "d20", "d1", "d50z")) {
        indexed <- sub("^d", "i", track)
        gtrack.copy(track, indexed)
        gtrack.convert_to_indexed(indexed)
    }
    gdb.reload()
    invisible(db)
}

# expr evaluated with the bulk path allowed and with it turned off, bit for bit
expect_same_both_paths <- function(expr) {
    q <- substitute(expr)
    env <- parent.frame()
    per_interval <- withr::with_envvar(c(MISHA_DISABLE_DENSE_BULK_READ = "1"), eval(q, env))
    bulk <- eval(q, env)
    expect_identical(bulk, per_interval)
    expect_identical(serialize(bulk, NULL), serialize(per_interval, NULL))
}

bulk_scope <- function() {
    gintervals(c("chr1", "chr1", "chr2", "chr3", "chrX"), c(1005, 7777, 3, 19990, 0), c(5003, 50001, 55555, 20000, 3))
}

test_that("a dense track at its own bin size reads the same in bulk", {
    local_bulk_db()
    withr::local_options(gmultitasking = FALSE)
    scope <- bulk_scope()

    for (shared in c(FALSE, TRUE)) {
        withr::local_options(gshared.track.backend = shared)
        expect_same_both_paths(gextract("d50a", gintervals.all()))
        expect_same_both_paths(gextract("d50a", scope, iterator = 50))
        expect_same_both_paths(gextract("d50a", scope, iterator = 50L))
        expect_same_both_paths(gextract("i50a", scope, iterator = "i50a"))
        expect_same_both_paths(gextract(c("d50a", "d50b"), gintervals.all(), iterator = "i50b"))
        expect_same_both_paths(gextract(c("d50a", "d50b", "i50a", "i50b"), scope))
        expect_same_both_paths(gextract("d50a * 2 + log(abs(i50b) + 1)", scope, iterator = 50))
        expect_same_both_paths(gextract(c("d50a", "ifelse(is.na(d50b), 0, d50b)", "pmax(d50a, i50a, na.rm = TRUE)"),
            gintervals.all(),
            colnames = c("a", "b", "c")
        ))
        expect_same_both_paths(gextract(c("d20", "i20"), gintervals.all()))
        expect_same_both_paths(gextract("d1 - i1", gintervals(c(1, 2), c(10, 99), c(100007, 5000)), iterator = 1))
        expect_same_both_paths(gextract(c("d50z", "i50z", "1 / d50z", "1 / i50z"), gintervals.all()))
        # intervals inside one bin, and overlapping ones
        expect_same_both_paths(gextract("d50a", gintervals("chr1", c(100, 120, 140, 1049), c(110, 130, 160, 1051)), iterator = 50))
        expect_same_both_paths(gextract("d50a", gintervals(c("chr1", "chr1", "chr2"), c(100, 120, 500), c(300, 251, 777)), iterator = 50))
        expect_same_both_paths(gextract("d50a", scope, iterator = 50, intervals_join = "intervals"))
        expect_same_both_paths(gextract(c("d50b", "i50b"), gintervals("chr3")))
        expect_same_both_paths(gextract("d50a", gintervals("chrX")))
        # the iterator intervals as an expression sees them, and an expression with no track
        itr <- "misha:::.misha$GITERATOR.INTERVALS"
        expect_same_both_paths(gextract(
            c("d50a", sprintf("%s$end - %s$start", itr, itr), sprintf("as.integer(%s$chrom) * 1e6 + %s$start", itr, itr)),
            scope,
            iterator = 50, colnames = c("v", "len", "pos")
        ))
        expect_same_both_paths(gextract("1 + 2", scope, iterator = 50))
        # an expression evaluated one interval at a time, and batches that end inside a scope interval
        expect_same_both_paths(gextract("sum(d50a)", scope, iterator = 50))
        expect_same_both_paths(withr::with_options(list(gbuf.size = 7), gextract(c("d50a", "i50b"), scope)))

        expect_same_both_paths(gsummary("d50z"))
        expect_same_both_paths(gsummary("1 / i50z"))
        expect_same_both_paths(gsummary("d1"))
        expect_same_both_paths(gquantiles("i50a", c(0.1, 0.5, 0.9), scope))
        expect_same_both_paths(gscreen("d50a > 0 & i50b < 0", scope))
        expect_same_both_paths(gdist("d50a", seq(-3, 3, by = 1), "i50b", c(-2, 0, 2), intervals = scope))
        expect_same_both_paths(gcor("d1", "i1", intervals = scope))
        expect_same_both_paths(gintervals.summary("d50a", scope, iterator = 50))
        expect_same_both_paths(gbins.summary("d50a", seq(-2, 2, by = 1), expr = "d50b"))
        expect_same_both_paths(gpartition("d50a", c(-1, 0, 1), gintervals.all()))
        expect_same_both_paths({
            set.seed(3)
            gsample("d50a", 100)
        })
        expect_same_both_paths({
            gtrack.create("bulk_created", "", "d50a * 3 + i50b", iterator = 50)
            on.exit(gtrack.rm("bulk_created", force = TRUE))
            gextract("bulk_created", gintervals.all())
        })
    }
})

test_that("a dense track reads the same in bulk in multitasking", {
    local_bulk_db()
    withr::local_options(gmultitasking = TRUE, gmax.processes = 4, gmin.scope4process = 1000)
    scope <- bulk_scope()

    for (shared in c(FALSE, TRUE)) {
        withr::local_options(gshared.track.backend = shared)
        expect_same_both_paths(gextract(c("d50a", "i50b", "d50z"), gintervals.all()))
        expect_same_both_paths(gextract("d1 + i1", scope, iterator = "d1"))
        expect_same_both_paths(gsummary("d50a + d50b", scope))
        expect_same_both_paths(gquantiles("d1", c(0.01, 0.5, 0.99)))
        expect_same_both_paths(gscreen("d20 > 0.5"))
        expect_same_both_paths(gdist("i50z", seq(-3, 3, by = 0.5)))
    }
})

test_that("intervals that do not lie in one bin of each track are read one by one", {
    local_bulk_db()
    withr::local_options(gmultitasking = FALSE)
    scope <- bulk_scope()
    withr::defer(for (v in c("v_shift", "v_exists", "v_avg", "v_kmer", "v_dist", "v_filter", "v_pv")) {
        try(gvtrack.rm(v), silent = TRUE)
    })
    gvtrack.create("v_shift", "d50a", "avg")
    gvtrack.iterator("v_shift", sshift = -50, eshift = 50)
    gvtrack.create("v_exists", "d50a", "exists")
    gvtrack.create("v_avg", "d50a", "avg")
    gvtrack.create("v_kmer", NULL, "kmer.count", kmer = "AA", strand = 1)
    gvtrack.create("v_dist", gintervals(c(1, 2), c(1000, 7000), c(1500, 9000)), "distance")
    gvtrack.create("v_filter", "d50a", "avg")
    gvtrack.filter("v_filter", gintervals("chr1", 2000, 3000))
    gvtrack.create("v_pv", "d50a", "global.percentile")

    # an iterator of another bin size, tracks of different bin sizes
    expect_same_both_paths(gextract("d50a", scope, iterator = 100))
    expect_same_both_paths(gextract("d50a", scope, iterator = 20))
    expect_same_both_paths(gextract(c("d50a", "d20"), scope, iterator = 50))
    expect_same_both_paths(gsummary("i50a", iterator = 100))
    # virtual tracks, other variables, other iterators
    expect_same_both_paths(gextract("v_shift", scope, iterator = 50))
    expect_same_both_paths(gextract(c("d50a", "v_exists"), scope, iterator = 50))
    expect_same_both_paths(gextract("v_avg", scope, iterator = 50))
    expect_same_both_paths(gextract(c("d50a", "v_kmer"), scope, iterator = 50))
    expect_same_both_paths(gextract(c("d50a", "v_dist"), scope, iterator = 50))
    expect_same_both_paths(gextract(c("d50a", "v_filter"), scope, iterator = 50))
    expect_same_both_paths(gextract(c("d50a", "v_pv"), scope, iterator = 50))
    expect_same_both_paths(gextract("d50a", scope, iterator = gintervals("chr1", c(100, 1049), c(110, 1051))))
    expect_same_both_paths(gscreen("v_shift > 0", scope, iterator = 50))
})
