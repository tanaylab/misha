# Third adversarial pass: targets fix CODE the first two passes never exercised -
# H2's 2D resize line (only 1D was tested), H3's gextract_multitask block + 2D/gzip
# branches (the second rewrite), H4's compute_n_mutations read site (only the
# edit_distance solvers were tested), H7's explicit-colnames split path, plus M6.

create_isolated_test_db()

test_that("ADV H2-2D: 2D quantiles bigset == in-memory on the non-multitask streaming path", {
    # exercises the 2D save (line 770) - pass 1 only covered the 1D save (line 755).
    # Uses a small subset here; the >gbuf.size truncation that this once worked
    # around (issue #149) is now fixed and covered by test-issue-149-*.R.
    withr::local_options(list(gmultitasking = FALSE, gbig.intervals.size = 5))
    full <- gscreen("test.rects > 40", gintervals.2d(c(1, 2), 0, -1))
    scope <- full[seq_len(min(80, nrow(full))), ]
    for (pct in list(0.5, c(0.2, 0.5, 0.9), c(0.1, 0.25, 0.75, 0.95))) {
        tmp <- paste0("test.aq2d_", sample(1:1e9, 1))
        gintervals.rm(tmp, force = TRUE)
        gintervals.quantiles("test.rects", pct, scope, iterator = c(1, 1), intervals.set.out = tmp)
        from_disk <- gintervals.load(tmp)
        in_mem <- gintervals.quantiles("test.rects", pct, scope, iterator = c(1, 1))
        ord <- function(d) d[order(d$chrom1, d$start1, d$chrom2, d$start2), ]
        expect_equal(ord(from_disk), ord(in_mem), ignore_attr = TRUE, info = paste(pct, collapse = ","))
        gintervals.rm(tmp, force = TRUE)
    }
})

test_that("ADV H3-multipath: gextract(file=) matches in-memory for 2D, gzip, and the multitask block", {
    parse_tsv <- function(path, con = path) {
        read.table(con, header = TRUE, sep = "\t", na.strings = "nan", check.names = FALSE, stringsAsFactors = FALSE)
    }
    # --- 2D file branch (use a scope that actually has rects) ---
    scope2d <- head(gscreen("test.rects > 40", gintervals.2d(1, 0, -1)), 50)
    inmem2d <- gextract("test.rects", scope2d, iterator = "test.rects")
    expect_gt(nrow(inmem2d), 0)
    f2 <- tempfile()
    withr::defer(unlink(f2))
    gextract("test.rects", scope2d, iterator = "test.rects", file = f2)
    ff2 <- parse_tsv(f2)
    o2a <- order(ff2$start1, ff2$start2)
    o2b <- order(inmem2d$start1, inmem2d$start2)
    expect_equal(nrow(ff2), nrow(inmem2d))
    expect_equal(ff2[["test.rects"]][o2a], inmem2d[["test.rects"]][o2b], tolerance = 1e-10)

    # --- gzip output ---
    scope <- gintervals(1, 0, 50000)
    inmem <- gextract("test.fixedbin", scope, iterator = 50)
    fz <- tempfile(fileext = ".gz")
    withr::defer(unlink(fz))
    gextract("test.fixedbin", scope, iterator = 50, file = fz)
    ffz <- parse_tsv(fz, gzfile(fz))
    expect_equal(ffz[["test.fixedbin"]], inmem[["test.fixedbin"]], tolerance = 1e-10)

    # --- C++ multitask block (gmultitasking=TRUE + file= -> tiles -> gextract_multitask) ---
    big <- gintervals(1, 0, 1e6)
    inmemM <- gextract("test.fixedbin", big, iterator = 1000)
    fm <- tempfile()
    withr::defer(unlink(fm))
    withr::with_options(
        list(gmultitasking = TRUE, gmax.processes = 3),
        gextract("test.fixedbin", big, iterator = 1000, file = fm)
    )
    ffm <- parse_tsv(fm)
    expect_equal(nrow(ffm), nrow(inmemM))
    o1 <- order(ffm$start)
    o2 <- order(inmemM$start)
    expect_equal(ffm[["test.fixedbin"]][o1], inmemM[["test.fixedbin"]][o2], tolerance = 1e-10)
})

test_that("ADV H4-modes: pwm.n_mutations and .pos work at L=70 (heap-table modes)", {
    remove_all_vtracks()
    consensus_pssm <- function(seqstr) {
        L <- nchar(seqstr)
        m <- matrix(0.02, L, 4)
        colnames(m) <- c("A", "C", "G", "T")
        for (i in seq_len(L)) m[i, substr(seqstr, i, i)] <- 0.94
        m / rowSums(m)
    }
    base <- 1000000
    sq <- toupper(gseq.extract(gintervals(1, base, base + 70)))
    skip_if(grepl("N", sq), "region has N bases")
    pssm <- consensus_pssm(sq)
    pos1 <- gintervals(1, base, base + 1)

    # n_mutations has its own read site (compute_n_mutations) over the [64] tables;
    # a perfect consensus already passes -> 0 mutations.
    gvtrack.create("nm70", NULL, "pwm.n_mutations", pssm = pssm, score.thresh = -40)
    expect_equal(gextract("nm70", pos1, iterator = pos1)$nm70, 0)

    # .pos goes through the same heap tables and must return a value (no OOB/crash).
    gvtrack.create("edp70", NULL, "pwm.edit_distance.pos", pssm = pssm, score.thresh = -40, max_edits = 5L)
    v <- gextract("edp70", pos1, iterator = pos1)$edp70
    expect_true(is.na(v) || is.finite(v))
})

test_that("ADV H7-colnames: track-parallel preserves explicit colnames in requested order", {
    exprs <- c(
        "test.fixedbin", "test.fixedbin+1", "test.fixedbin+2",
        "test.fixedbin+3", "test.fixedbin+4", "test.fixedbin+5"
    )
    cn <- c("aa", "bb", "cc", "dd", "ee", "ff")
    scope <- gintervals(1, 0, 100000)
    res <- withr::with_options(
        list(gmultitasking = TRUE, gmultitasking.strategy = "tracks", gmax.processes = 3),
        gextract(exprs, scope, colnames = cn)
    )
    val_names <- setdiff(names(res), c("chrom", "start", "end", "intervalID"))
    expect_identical(val_names, cn) # custom names, in requested order
    # names stay attached to the right data
    base <- res[["aa"]]
    for (k in 2:6) expect_equal(res[[cn[k]]], base + (k - 1))
})

test_that("ADV M6: gtrack.import_set names tracks from the file basename, not a dotted parent dir", {
    d <- file.path(tempdir(), paste0("GSE777.RAW_", sample(1:1e6, 1)))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    withr::defer(unlink(d, recursive = TRUE))
    for (nm in c("sampleA.bed", "sampleB.bed")) {
        writeLines(c("chr1\t0\t100\tx\t5", "chr1\t100\t200\ty\t7"), file.path(d, nm))
    }
    withr::defer({
        gtrack.rm("sampleA", force = TRUE)
        gtrack.rm("sampleB", force = TRUE)
    })
    failed <- gtrack.import_set("desc", file.path(d, "*.bed"), binsize = 50, track.prefix = "")
    expect_true(gtrack.exists("sampleA")) # not collapsed to "GSE777"
    expect_true(gtrack.exists("sampleB"))
})

test_that("ADV H8-rm: gtrack.rm from a database subdirectory marks the cache dirty", {
    iv <- gintervals(1, 0, 10000)
    gdir.create("rmsub", showWarnings = FALSE)
    gtrack.create_sparse("rmsub.tA", "a", iv, 1)
    gtrack.create_sparse("rmsub.tB", "b", iv, 1)
    gdb.reload(rescan = TRUE)
    withr::defer(suppressMessages(gdb.init(get("GROOT", envir = .misha))))
    gdir.cd("rmsub")
    gtrack.rm("tB", force = TRUE) # remove while in the subdir
    expect_true(misha:::.gdb.cache_is_dirty(get("GROOT", envir = .misha)))
})
