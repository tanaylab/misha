# gtrack.2d.import and gtrack.2d.import_contacts read the input files and build the
# per-chromosome-pair quad trees in forked kids when gmultitasking is on. Each test runs the
# same import with multitasking off and on and requires byte-identical track files and
# identical gextract output. The parallel runs also assert that kids were really forked
# (their CPU time shows up in proc.time()), so a silent fallback to the serial path cannot
# pass for parallel. These tests compare serial with parallel on this branch only: that the
# serial path still writes the same track files as before the parallel import was added
# (only the intermediate file names changed) was checked by md5 against master on the PR's
# benchmark inputs, not here.

ensure_valid_groot()

import_chroms <- c(chr1 = 1e6, chr2 = 8e5, chr3 = 5e5)

setup_import_db <- function() {
    create_test_db("import_db", chrom_sizes = data.frame(chrom = names(import_chroms), size = import_chroms))
    gdb.init("import_db")
}

# Random points (or non-overlapping rectangles on a 1kb grid when rects = TRUE), unique per
# chromosome pair, plus a few rows on chromosomes the database does not have.
random_intervals <- function(n, seed, rects = FALSE, unknown = 5) {
    set.seed(seed)
    chrom1 <- sample(names(import_chroms), n, replace = TRUE)
    chrom2 <- sample(names(import_chroms), n, replace = TRUE)
    if (rects) {
        start1 <- floor(runif(n) * (import_chroms[chrom1] / 1000 - 1)) * 1000
        start2 <- floor(runif(n) * (import_chroms[chrom2] / 1000 - 1)) * 1000
        width <- 500
    } else {
        start1 <- floor(runif(n) * (import_chroms[chrom1] - 1))
        start2 <- floor(runif(n) * (import_chroms[chrom2] - 1))
        width <- 1
    }
    d <- data.frame(
        chrom1 = chrom1, start1 = as.integer(start1), end1 = as.integer(start1 + width),
        chrom2 = chrom2, start2 = as.integer(start2), end2 = as.integer(start2 + width),
        value = round(runif(n, 0, 100), 3)
    )
    d <- d[!duplicated(d[, c("chrom1", "start1", "chrom2", "start2")]), ]
    unk <- d[seq_len(unknown), ]
    unk$chrom1 <- rep("chrUn_1", unknown)
    rbind(d, unk)
}

# Writes consecutive slices of d into length(fracs) files and returns their paths
write_slices <- function(d, fracs, prefix) {
    ends <- round(cumsum(fracs) / sum(fracs) * nrow(d))
    starts <- c(1, head(ends, -1) + 1)
    vapply(seq_along(fracs), function(i) {
        path <- file.path(getwd(), sprintf("%s_%d.tsv", prefix, i))
        write.table(d[starts[i]:ends[i], ], path, sep = "\t", quote = FALSE, row.names = FALSE)
        path
    }, character(1))
}

track_files_md5 <- function(track) {
    files <- list.files(.track_dir(track), full.names = TRUE)
    setNames(unname(tools::md5sum(files)), basename(files))
}

# misha's lifetime invariant (src/rdbutils.h, C6): each import runs two RdbInitializer scopes
# in turn, and both must be unwound on success and on error alike.
expect_counters_reset <- function() {
    expect_identical(.glifetime_counters(), c(ref_count = 0L, protect_count = 0L))
}

# Runs import(track) with multitasking off and on. Returns, per mode, the md5 of the track
# files, the gextract over the whole genome, the leftover hidden files and the kids' CPU time.
import_both_ways <- function(import, max_data_size = NULL, max_mem_usage = NULL) {
    res <- list()
    for (mt in c(FALSE, TRUE)) {
        track <- if (mt) "imp_parallel" else "imp_serial"
        opts <- list(gmultitasking = mt, gmax.processes = 4)
        if (!is.null(max_data_size)) {
            opts$gmax.data.size <- max_data_size
        }
        if (!is.null(max_mem_usage)) {
            opts$gmax.mem.usage <- max_mem_usage
        }
        if (gtrack.exists(track)) {
            gtrack.rm(track, force = TRUE)
        }
        t0 <- proc.time()
        withr::with_options(opts, import(track))
        dt <- proc.time() - t0
        expect_counters_reset()
        hidden <- list.files(.track_dir(track), all.files = TRUE, no.. = TRUE, pattern = "^\\.")
        res[[if (mt) "parallel" else "serial"]] <- list(
            md5 = track_files_md5(track),
            data = gextract(track, gintervals.2d.all(), colnames = "v"),
            hidden = setdiff(hidden, ".attributes"),
            kids_cpu = dt[["user.child"]] + dt[["sys.child"]]
        )
    }
    res
}

# import_both_ways, plus what stage 3 of the parallel run reports through MISHA_2D_IMPORT_STATS:
# the budget and the largest estimate of a pair in bytes, and the most memory reserved and pairs
# built at once
import_with_stats <- function(...) {
    stats_file <- withr::local_tempfile(fileext = ".tsv")
    withr::local_envvar(MISHA_2D_IMPORT_STATS = stats_file)
    res <- import_both_ways(...)
    res$stats <- read.delim(stats_file)
    res
}

# The pairs built at once stayed within the budget, and were at most max_running
expect_budget_kept <- function(stats, budget_kb, max_running) {
    expect_equal(stats$budget, budget_kb * 1000)
    expect_lte(stats$max_reserved, stats$budget)
    expect_lte(stats$max_running, max_running)
}

expect_same_both_ways <- function(res) {
    expect_gt(length(res$serial$md5), 0)
    expect_identical(res$parallel$md5, res$serial$md5)
    expect_identical(res$parallel$data, res$serial$data)
    expect_length(res$serial$hidden, 0)
    expect_length(res$parallel$hidden, 0)
    expect_gt(res$parallel$kids_cpu, 0)
}

sort_2d <- function(d) {
    d <- d[order(d$chrom1, d$start1, d$chrom2, d$start2), ]
    rownames(d) <- NULL
    d
}

# What gtrack.2d.import_contacts should produce: contacts at the interval centers, duplicates
# summed, every contact also stored mirrored.
expected_contacts <- function(d) {
    d <- d[d$chrom1 %in% names(import_chroms) & d$chrom2 %in% names(import_chroms), ]
    x1 <- (d$start1 + d$end1) %/% 2
    x2 <- (d$start2 + d$end2) %/% 2
    fwd <- data.frame(chrom1 = d$chrom1, start1 = x1, chrom2 = d$chrom2, start2 = x2, v = d$value)
    self <- d$chrom1 == d$chrom2 & x1 == x2
    rev <- data.frame(chrom1 = d$chrom2, start1 = x2, chrom2 = d$chrom1, start2 = x1, v = d$value)[!self, ]
    agg <- aggregate(v ~ chrom1 + start1 + chrom2 + start2, rbind(fwd, rev), sum)
    sort_2d(agg)
}

test_that("gtrack.2d.import of points from several files is identical in parallel", {
    local_db_state()
    withr::with_tempdir({
        setup_import_db()
        d <- random_intervals(60000, seed = 1)
        files <- write_slices(d, c(5, 1, 3, 2), "points")
        res <- import_both_ways(function(track) gtrack.2d.import(track, "test", files))
        expect_same_both_ways(res)

        expect_equal(gtrack.info("imp_parallel")$type, "points")
        known <- d[d$chrom1 %in% names(import_chroms), ]
        got <- sort_2d(res$parallel$data)
        want <- sort_2d(known)
        expect_equal(nrow(got), nrow(want))
        expect_equal(got$start1, want$start1)
        expect_equal(got$start2, want$start2)
        expect_equal(got$v, want$value, tolerance = 1e-6)
    })
})

test_that("gtrack.2d.import of rectangles is identical in parallel, also when one file holds only points", {
    local_db_state()
    withr::with_tempdir({
        setup_import_db()
        d <- random_intervals(30000, seed = 2, rects = TRUE)
        files <- write_slices(d, c(2, 3, 1), "rects")
        # one file of points only: its kid sees nothing but points, the others do not
        pts <- random_intervals(200, seed = 3, unknown = 0)
        pts$start1 <- pts$start1 %/% 1000L * 1000L + 700L
        pts$end1 <- pts$start1 + 1L
        pts <- pts[!duplicated(pts[, c("chrom1", "start1", "chrom2", "start2")]), ]
        points_file <- file.path(getwd(), "only_points.tsv")
        write.table(pts, points_file, sep = "\t", quote = FALSE, row.names = FALSE)
        files <- c(points_file, files)

        res <- import_both_ways(function(track) gtrack.2d.import(track, "test", files))
        expect_same_both_ways(res)
        expect_equal(gtrack.info("imp_parallel")$type, "rectangles")
        expect_equal(nrow(res$parallel$data), sum(d$chrom1 %in% names(import_chroms)) + nrow(pts))
    })
})

test_that("gtrack.2d.import with several subtrees per pair is identical in parallel", {
    local_db_state()
    withr::with_tempdir({
        setup_import_db()
        d <- random_intervals(40000, seed = 4)
        files <- write_slices(d, c(1, 1), "subtrees")
        # well below the records of the cis pairs, so they are split into subtrees
        res <- import_both_ways(function(track) gtrack.2d.import(track, "test", files), max_data_size = 1000)
        expect_same_both_ways(res)
        expect_equal(nrow(res$parallel$data), sum(d$chrom1 %in% names(import_chroms)))
    })
})

test_that("gtrack.2d.import_contacts sums duplicates across files and mirrors cis contacts in parallel", {
    local_db_state()
    withr::with_tempdir({
        setup_import_db()
        d <- random_intervals(50000, seed = 5)
        # the same contacts again, in both orientations, in other files: summed across kids
        dup <- d[sample(nrow(d), 3000), ]
        dup_rev <- dup[1:1500, c("chrom2", "start2", "end2", "chrom1", "start1", "end1", "value")]
        names(dup_rev) <- names(d)
        all <- rbind(d, dup, dup_rev)
        files <- c(write_slices(d, c(3, 1, 2), "contacts"), write_slices(rbind(dup, dup_rev), c(1, 1), "dups"))

        res <- import_both_ways(function(track) gtrack.2d.import_contacts(track, "test", files))
        expect_same_both_ways(res)

        got <- sort_2d(res$parallel$data[, c("chrom1", "start1", "chrom2", "start2", "v")])
        got$chrom1 <- as.character(got$chrom1)
        got$chrom2 <- as.character(got$chrom2)
        want <- expected_contacts(all)
        expect_equal(nrow(got), nrow(want))
        expect_equal(got[, 1:4], want[, 1:4], ignore_attr = TRUE)
        expect_equal(got$v, want$v, tolerance = 1e-5)
    })
})

test_that("gtrack.2d.import_contacts with several subtrees per pair is identical in parallel", {
    local_db_state()
    withr::with_tempdir({
        setup_import_db()
        d <- random_intervals(40000, seed = 6)
        d <- rbind(d, d[1:2000, ])
        files <- write_slices(d, c(1, 2, 1), "contacts")
        res <- import_both_ways(function(track) gtrack.2d.import_contacts(track, "test", files), max_data_size = 1000)
        expect_same_both_ways(res)
        expect_equal(nrow(res$parallel$data), nrow(expected_contacts(d)))
    })
})

# gmax.mem.usage (KB) caps the estimated memory of the pairs built at once: a fixed 4 MiB per pair
# plus its records. A 1 MB budget is below every pair, so the pairs are built one at a time. The
# larger budgets fit two of these pairs (contacts pairs are larger), and the other kids wait.
# The 4 kids start together, so a budget that were not kept would show up in the stats.
test_that("gtrack.2d.import under a memory budget is identical in parallel", {
    local_db_state()
    withr::with_tempdir({
        setup_import_db()
        points <- write_slices(random_intervals(40000, seed = 11), c(2, 1, 1), "points")
        rects <- write_slices(random_intervals(20000, seed = 12, rects = TRUE), c(1, 1), "rects")
        for (budget in c(1000, 10000)) {
            for (files in list(points, rects)) {
                res <- import_with_stats(function(track) gtrack.2d.import(track, "test", files), max_mem_usage = budget)
                expect_same_both_ways(res)
                expect_budget_kept(res$stats, budget, if (budget == 1000) 1 else 2)
            }
        }
    })
})

test_that("gtrack.2d.import_contacts under a memory budget is identical in parallel", {
    local_db_state()
    withr::with_tempdir({
        setup_import_db()
        d <- random_intervals(40000, seed = 13)
        files <- write_slices(rbind(d, d[1:3000, ]), c(1, 2, 1), "contacts")
        for (budget in c(1000, 13000)) {
            res <- import_with_stats(function(track) gtrack.2d.import_contacts(track, "test", files), max_mem_usage = budget)
            expect_same_both_ways(res)
            expect_budget_kept(res$stats, budget, if (budget == 1000) 1 else 2)
            expect_equal(nrow(res$parallel$data), nrow(expected_contacts(rbind(d, d[1:3000, ]))))
        }
        # with subtrees too
        res <- import_with_stats(function(track) gtrack.2d.import_contacts(track, "test", files), max_data_size = 1000, max_mem_usage = 1000)
        expect_same_both_ways(res)
        expect_budget_kept(res$stats, 1000, 1)
    })
})

# A pair split into subtrees is built one subtree at a time, so its estimate counts the records of
# its largest subtree, not of the whole pair: lowering gmax.data.size lets more pairs fit the budget.
test_that("the memory estimate of a pair split into subtrees is that of its largest subtree", {
    local_db_state()
    withr::with_tempdir({
        setup_import_db()
        points <- write_slices(random_intervals(40000, seed = 14), c(1, 1), "points")
        rects <- write_slices(random_intervals(20000, seed = 15, rects = TRUE), c(1, 1), "rects")
        for (import in list(
            function(track) gtrack.2d.import(track, "test", points),
            function(track) gtrack.2d.import(track, "test", rects),
            function(track) gtrack.2d.import_contacts(track, "test", points)
        )) {
            whole <- import_with_stats(import)
            parts <- import_with_stats(import, max_data_size = 500)
            expect_same_both_ways(parts)
            expect_lt(parts$stats$max_mem, whole$stats$max_mem)
        }
    })
})

test_that("gmax.mem.usage given as an integer is not truncated", {
    local_db_state()
    withr::with_tempdir({
        setup_import_db()
        points <- write_slices(random_intervals(2000, seed = 16), c(1, 1), "points")
        # 4 GB as KB: overflows a 32-bit int multiplied by 1000
        res <- import_with_stats(function(track) gtrack.2d.import(track, "test", points), max_mem_usage = 4000000L)
        expect_equal(res$stats$budget, 4e9)
    })
})

test_that("gtrack.2d.import_contacts with allow.duplicates = FALSE errors the same way in parallel", {
    local_db_state()
    withr::with_tempdir({
        setup_import_db()
        d <- random_intervals(20000, seed = 7, unknown = 0)
        files <- write_slices(d, c(1, 1, 1), "nodup")

        # no duplicates: imports
        res <- import_both_ways(function(track) gtrack.2d.import_contacts(track, "test", files, allow.duplicates = FALSE))
        expect_same_both_ways(res)

        # one contact repeated in another file, the second time with its ends swapped
        dup <- d[100, c("chrom2", "start2", "end2", "chrom1", "start1", "end1", "value")]
        names(dup) <- names(d)
        dup_file <- file.path(getwd(), "dup.tsv")
        write.table(dup, dup_file, sep = "\t", quote = FALSE, row.names = FALSE)

        errs <- vapply(c(FALSE, TRUE), function(mt) {
            track <- paste0("dup_", mt)
            err <- withr::with_options(list(gmultitasking = mt, gmax.processes = 4), tryCatch(
                gtrack.2d.import_contacts(track, "test", c(files, dup_file), allow.duplicates = FALSE),
                error = conditionMessage
            ))
            expect_counters_reset()
            expect_false(gtrack.exists(track))
            err
        }, character(1))
        expect_match(errs[1], "Duplicated contact")
        expect_identical(errs[2], errs[1])
    })
})

test_that("gtrack.2d.import_contacts with fends is identical in parallel", {
    local_db_state()
    withr::with_tempdir({
        setup_import_db()
        set.seed(8)
        nfends <- 5000
        chr <- sample(names(import_chroms), nfends, replace = TRUE)
        fends <- data.frame(
            fend = seq_len(nfends), chr = chr,
            coord = as.integer(floor(runif(nfends) * (import_chroms[chr] - 1)))
        )
        fends_file <- file.path(getwd(), "redb.fends")
        write.table(fends, fends_file, sep = "\t", quote = FALSE, row.names = FALSE)

        n <- 40000
        contacts <- data.frame(
            fend1 = sample(nfends + 50, n, replace = TRUE), # fends above nfends are undefined: skipped
            fend2 = sample(nfends + 50, n, replace = TRUE),
            count = sample(1:5, n, replace = TRUE)
        )
        files <- write_slices(contacts, c(2, 1, 1), "fends_contacts")

        res <- import_both_ways(function(track) gtrack.2d.import_contacts(track, "test", files, fends_file))
        expect_same_both_ways(res)

        ok <- contacts$fend1 <= nfends & contacts$fend2 <= nfends
        f1 <- fends[contacts$fend1[ok], ]
        f2 <- fends[contacts$fend2[ok], ]
        as_intervals <- data.frame(
            chrom1 = f1$chr, start1 = f1$coord, end1 = f1$coord + 1L,
            chrom2 = f2$chr, start2 = f2$coord, end2 = f2$coord + 1L, value = contacts$count[ok]
        )
        expect_equal(sum(res$parallel$data$v), sum(expected_contacts(as_intervals)$v), tolerance = 1e-6)
    })
})

test_that("malformed input gives the serial error message in parallel", {
    local_db_state()
    withr::with_tempdir({
        setup_import_db()
        d <- random_intervals(9000, seed = 9, unknown = 0)
        files <- write_slices(d, c(1, 1, 1), "bad")
        lines <- readLines(files[2])
        lines[500] <- sub("\t[^\t]*$", "\tnot_a_number", lines[500])
        writeLines(lines, files[2])

        for (import in list(gtrack.2d.import, gtrack.2d.import_contacts)) {
            errs <- vapply(c(FALSE, TRUE), function(mt) {
                err <- withr::with_options(list(gmultitasking = mt, gmax.processes = 4), tryCatch(
                    import("bad_track", "test", files),
                    error = conditionMessage
                ))
                expect_counters_reset()
                err
            }, character(1))
            expect_match(errs[1], sprintf("File %s, line 500: invalid value", files[2]), fixed = TRUE)
            expect_identical(errs[2], errs[1])
            expect_false(gtrack.exists("bad_track"))
        }
    })
})

test_that("parallel gtrack.2d.import_contacts is converted in an indexed database", {
    local_db_state()
    tmp_root <- withr::local_tempdir()
    fasta <- file.path(tmp_root, "genome.fasta")
    cat(unlist(lapply(names(import_chroms), function(chr) {
        c(">", chr, "\n", strrep("A", import_chroms[[chr]]), "\n")
    })), sep = "", file = fasta)
    db_path <- file.path(tmp_root, "testdb")

    withr::with_options(list(gmulticontig.indexed_format = TRUE), {
        gdb.create(groot = db_path, fasta = fasta, verbose = FALSE)
        gdb.init(db_path)
        expect_true(.gdb.is_indexed())

        withr::with_dir(tmp_root, {
            d <- random_intervals(20000, seed = 10)
            files <- write_slices(d, c(1, 1), "contacts")
            res <- import_both_ways(function(track) gtrack.2d.import_contacts(track, "test", files))
        })
        expect_true(file.exists(file.path(.track_dir("imp_parallel"), "track.idx")))
        expect_same_both_ways(res)
    })
})
