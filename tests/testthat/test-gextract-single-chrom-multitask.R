create_isolated_test_db()

test_that("single-chrom fixed-bin multitask split matches serial output", {
    intervals <- gintervals("chr1", 0, 5000000)

    withr::local_options(list(
        gmax.processes = 8L,
        gmax.processes2core = 4L,
        gmax.data.size = 1e9,
        gmin.scope4process = 1L
    ))

    withr::local_options(list(gmultitasking = FALSE))
    serial <- gextract("1", intervals = intervals, iterator = 20, colnames = "value")

    withr::local_options(list(gmultitasking = TRUE))
    parallel <- gextract("1", intervals = intervals, iterator = 20, colnames = "value")

    ord_serial <- order(serial$chrom, serial$start, serial$end, serial$intervalID)
    ord_parallel <- order(parallel$chrom, parallel$start, parallel$end, parallel$intervalID)

    serial <- serial[ord_serial, ]
    parallel <- parallel[ord_parallel, ]

    expect_equal(parallel, serial)
})

test_that("single-chrom shifted vtrack extraction matches serial output", {
    track_name <- "test.mt.single.chrom.sum"
    if (track_name %in% gvtrack.ls()) {
        gvtrack.rm(track_name)
    }
    on.exit(gvtrack.rm(track_name), add = TRUE)

    gvtrack.create(track_name, src = "test.fixedbin", func = "sum")
    gvtrack.iterator(track_name, sshift = -100, eshift = 100)

    intervals <- gintervals("chr1", 0, 5000000)

    withr::local_options(list(
        gmax.processes = 8L,
        gmax.processes2core = 4L,
        gmax.data.size = 1e9,
        gmin.scope4process = 1L
    ))

    withr::local_options(list(gmultitasking = FALSE))
    serial <- gextract(track_name, intervals = intervals, iterator = 20)

    withr::local_options(list(gmultitasking = TRUE))
    parallel <- gextract(track_name, intervals = intervals, iterator = 20)

    ord_serial <- order(serial$chrom, serial$start, serial$end, serial$intervalID)
    ord_parallel <- order(parallel$chrom, parallel$start, parallel$end, parallel$intervalID)

    serial <- serial[ord_serial, ]
    parallel <- parallel[ord_parallel, ]

    expect_equal(parallel[, c("chrom", "start", "end", "intervalID")], serial[, c("chrom", "start", "end", "intervalID")])
    expect_equal(parallel[[track_name]], serial[[track_name]], tolerance = 1e-6)
})
