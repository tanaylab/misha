load_test_db()

# Tests for parallel gdb.convert_to_indexed (threads argument). Each test
# builds a small isolated per-chromosome database, populates it with several
# tracks (1D + 2D), and converts the whole DB. mclapply requires fork so the
# parallel-path tests skip on non-Unix; the threads=1 path runs everywhere.

# Helper: create a per-chromosome test DB with several convertible tracks.
# Returns a list with the db path and the names of the created tracks.
.make_per_chrom_db_with_tracks <- function(n_1d = 3L, n_2d = 2L) {
    tmp_root <- tempfile(pattern = "misha_par_convert_")
    dir.create(tmp_root, recursive = TRUE)

    # Per-chrom layout: one FASTA per chromosome. gdb.create with
    # gmulticontig.indexed_format=FALSE uses the FASTA basename (sans
    # extension) as the chromosome name, so we name files chr1.fa, chr2.fa.
    chr1_fa <- file.path(tmp_root, "chr1.fa")
    chr2_fa <- file.path(tmp_root, "chr2.fa")
    cat(">chr1\n", paste0(rep("ACGT", 200), collapse = ""), "\n",
        sep = "", file = chr1_fa
    )
    cat(">chr2\n", paste0(rep("ACGT", 200), collapse = ""), "\n",
        sep = "", file = chr2_fa
    )

    db_path <- file.path(tmp_root, "testdb")

    # Force per-chromosome layout (the format gdb.convert_to_indexed migrates
    # away from). We then convert to indexed without remove_old_files so the
    # seq directory keeps both representations.
    withr::with_options(list(gmulticontig.indexed_format = FALSE), {
        gdb.create(groot = db_path, fasta = c(chr1_fa, chr2_fa), verbose = FALSE)
        gdb.init(db_path)

        # Sparse 1D tracks
        sparse_names <- paste0("sparse_", seq_len(n_1d))
        for (nm in sparse_names) {
            intervs <- data.frame(
                chrom = c("chr1", "chr1", "chr2"),
                start = c(0L, 100L, 50L),
                end = c(50L, 200L, 150L),
                stringsAsFactors = FALSE
            )
            gtrack.create_sparse(nm, "sparse test track", intervs, values = c(1.0, 2.0, 3.0))
        }

        # 2D rectangles tracks
        rect_names <- paste0("rect_", seq_len(n_2d))
        for (nm in rect_names) {
            intervs2d <- gintervals.2d(
                chroms1 = c("chr1", "chr2"),
                starts1 = c(10L, 20L),
                ends1 = c(40L, 60L),
                chroms2 = c("chr1", "chr2"),
                starts2 = c(100L, 200L),
                ends2 = c(130L, 240L)
            )
            gtrack.2d.create(nm, "2d test track", intervs2d, c(1.0, 2.0))
        }

        list(
            tmp_root = tmp_root,
            db_path = db_path,
            sparse_tracks = sparse_names,
            rect_tracks = rect_names
        )
    })
}

# Helper: was a track converted to indexed format?
.is_track_indexed <- function(db_path, track_name) {
    trackdir <- file.path(db_path, "tracks", paste0(gsub("\\.", "/", track_name), ".track"))
    file.exists(file.path(trackdir, "track.idx"))
}

test_that("gdb.convert_to_indexed with threads > 1 converts all tracks", {
    skip_if(.Platform$OS.type != "unix")

    local_db_state()
    handles <- .make_per_chrom_db_with_tracks(n_1d = 3L, n_2d = 2L)
    withr::defer(unlink(handles$tmp_root, recursive = TRUE))

    gdb.convert_to_indexed(
        groot = handles$db_path,
        convert_tracks = TRUE,
        remove_old_files = FALSE,
        force = TRUE,
        verbose = FALSE,
        threads = 2L
    )

    gdb.init(handles$db_path)
    for (nm in c(handles$sparse_tracks, handles$rect_tracks)) {
        expect_true(.is_track_indexed(handles$db_path, nm),
            info = sprintf("track %s should be indexed", nm)
        )
    }
})

test_that("threads = 1 uses the serial path and still converts all tracks", {
    local_db_state()
    handles <- .make_per_chrom_db_with_tracks(n_1d = 2L, n_2d = 1L)
    withr::defer(unlink(handles$tmp_root, recursive = TRUE))

    gdb.convert_to_indexed(
        groot = handles$db_path,
        convert_tracks = TRUE,
        remove_old_files = FALSE,
        force = TRUE,
        verbose = FALSE,
        threads = 1L
    )

    gdb.init(handles$db_path)
    for (nm in c(handles$sparse_tracks, handles$rect_tracks)) {
        expect_true(.is_track_indexed(handles$db_path, nm),
            info = sprintf("track %s should be indexed", nm)
        )
    }
})

test_that("parallel conversion reports per-track failures without aborting batch", {
    skip_if(.Platform$OS.type != "unix")

    local_db_state()
    handles <- .make_per_chrom_db_with_tracks(n_1d = 3L, n_2d = 0L)
    withr::defer({
        # Restore permissions before unlink so cleanup succeeds.
        victim_dir <- file.path(
            handles$db_path, "tracks",
            paste0(handles$sparse_tracks[1], ".track")
        )
        if (dir.exists(victim_dir)) Sys.chmod(victim_dir, "0700")
        unlink(handles$tmp_root, recursive = TRUE)
    })

    # Make one track unconvertible by removing write permission on its
    # directory. The converter writes track.dat.tmp inside the track dir,
    # so this triggers a "Permission denied" error from the C++ side
    # without affecting the other tracks.
    victim <- handles$sparse_tracks[1]
    victim_dir <- file.path(handles$db_path, "tracks", paste0(victim, ".track"))
    expect_true(dir.exists(victim_dir))
    Sys.chmod(victim_dir, "0500")

    expect_warning(
        gdb.convert_to_indexed(
            groot = handles$db_path,
            convert_tracks = TRUE,
            remove_old_files = FALSE,
            force = TRUE,
            verbose = FALSE,
            threads = 2L
        ),
        regexp = victim
    )

    # Restore so subsequent operations work.
    Sys.chmod(victim_dir, "0700")
    gdb.init(handles$db_path)

    # Surviving tracks must have been converted.
    survivors <- setdiff(handles$sparse_tracks, victim)
    for (nm in survivors) {
        expect_true(.is_track_indexed(handles$db_path, nm),
            info = sprintf("survivor track %s should still be indexed", nm)
        )
    }

    # Victim track itself must NOT be indexed (the failure prevented it).
    expect_false(.is_track_indexed(handles$db_path, victim),
        info = "victim track must remain unconverted on failure"
    )
})

test_that("threads default resolves to a sensible value", {
    # NULL -> min(detectCores, 8), or 1 on non-Unix.
    n <- misha:::.gdb.convert_to_indexed.resolve_threads(NULL)
    expect_type(n, "integer")
    expect_gte(n, 1L)
    expect_lte(n, 8L)

    # Explicit numeric is coerced to integer.
    expect_identical(misha:::.gdb.convert_to_indexed.resolve_threads(2), 2L)
    expect_identical(misha:::.gdb.convert_to_indexed.resolve_threads(1L), 1L)

    # Invalid -> error.
    expect_error(misha:::.gdb.convert_to_indexed.resolve_threads(0))
    expect_error(misha:::.gdb.convert_to_indexed.resolve_threads(-1))
    expect_error(misha:::.gdb.convert_to_indexed.resolve_threads("hello"))
})

# Restore the test database after this file's tests.
suppressMessages(gdb.init("/net/mraid20/export/tgdata/db/tgdb/misha_test_db/"))
