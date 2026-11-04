create_isolated_test_db()

test_that("gtrack.create_pwm_energy works", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    allgenome <- .misha$ALLGENOME
    # The snapshot predates the rule that a window holding N (or any non-ACGT base) is not scored,
    # so it is compared only on bins none of whose windows (motif 0 is 33 bp) reach such a base
    chroms <- unique(as.character(allgenome[[1]]$chrom))
    seqs <- vapply(chroms, function(cc) toupper(gseq.extract(gintervals(cc, 0, 1000040))), character(1))
    .misha$ALLGENOME[[1]]$end <- 1000000
    .misha$ALLGENOME[[2]]$end1 <- 1000000
    .misha$ALLGENOME[[2]]$end2 <- 1000000
    gtrack.rm(tmptrack, force = TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.create_pwm_energy(tmptrack, "", "misha_motifs", 0, 0.02, iterator = 50)
    r <- gextract(tmptrack, .misha$ALLGENOME, colnames = "test.tmptrack")
    .misha$ALLGENOME <- allgenome
    old <- readr::read_rds("/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/gtrack.create_pwm_energy.rds")
    r <- r[order(as.character(r$chrom), r$start), ]
    old <- old[order(as.character(old$chrom), old$start), ]
    expect_equal(r$start, old$start)
    acgt_only <- !mapply(
        function(cc, s, e) grepl("[^ACGT]", substr(seqs[[cc]], s + 1, e + 33)),
        as.character(r$chrom), r$start, r$end
    )
    expect_gt(sum(acgt_only), 300000)
    expect_equal(r$test.tmptrack[acgt_only], old$test.tmptrack[acgt_only], tolerance = 1e-5)
})

test_that("gtrack.create_pwm_energy multitasking matches sequential on an indexed database", {
    # gtrack.create_pwm_energy used to open the genome sequence before forking, and forked
    # processes share the file offset of an inherited descriptor -- on an indexed database, where
    # set_seqdir() opens genome.seq immediately, the workers then handed each other's bytes back.
    # Reproducing the race needs slow (networked) reads, so this test is a structural guard on the
    # multitasking path rather than a reliable detector; see the warning in GenomeSeqFetch.h.
    test_fasta <- tempfile(fileext = ".fasta")
    set.seed(17)
    seqs <- vapply(1:8, function(i) paste(sample(c("A", "C", "G", "T"), 2e5, replace = TRUE), collapse = ""), character(1))
    cat(paste0(">chr", seq_along(seqs), "\n", seqs, "\n", collapse = ""), file = test_fasta)

    test_db <- tempfile()
    original_groot <- .misha$GROOT
    withr::defer({
        # don't leave GROOT pointing at the temporary DB we are about to delete
        if (!is.null(original_groot) && dir.exists(original_groot)) {
            suppressMessages(gdb.init(original_groot))
        }
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    withr::with_options(list(gmulticontig.indexed_format = TRUE), {
        gdb.create(groot = test_db, fasta = test_fasta, verbose = FALSE)
        gdb.init(test_db)

        dir.create(file.path(test_db, "pssms"), showWarnings = FALSE)
        file.copy(
            list.files(file.path(shared_test_db_path(), "pssms"), full.names = TRUE),
            file.path(test_db, "pssms")
        )

        energies <- function(multitasking) {
            track <- paste0("pwmE", if (multitasking) "par" else "ser")
            withr::with_options(list(gmultitasking = multitasking), {
                gtrack.create_pwm_energy(track, "", "misha_motifs", 0, 0.02, iterator = 50)
            })
            on.exit(gtrack.rm(track, force = TRUE))
            gextract(track, gintervals.all(), colnames = "v")$v
        }

        expect_equal(energies(TRUE), energies(FALSE))
    })
})
