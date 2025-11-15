test_that("gintervals.load_chain handles source overlaps with 'error' policy", {
    local_db_state()

    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg38")
    gdb.reload()

    set.seed(60427)
    random_intervals <- gintervals.random(size = 500, n = 1e4)

    chain_file <- "/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19ToHg38.over.chain"
    chain <- gintervals.load_chain(chain_file)

    expect_regression(chain, "liftover-ucsc-random-chain")
})
