# Tests for 2D indexed track support (RECTS and POINTS)
# These tests verify that converting 2D tracks to indexed format (track.dat + track.idx)
# produces identical query results compared to the original per-pair format.

create_isolated_test_db()
gdir.create("temp", showWarnings = FALSE)

# Skip if already in indexed format (these tests convert tracks themselves)
skip_if(
    getOption("gmulticontig.indexed_format", FALSE) || gdb.info(.misha$GROOT)$format == "indexed",
    "Indexed format enabled, set gmulticontig.indexed_format = FALSE to run this test"
)

# ============================================================================
# Category 1: Conversion round-trip for RECTS tracks
# ============================================================================

test_that("RECTS track round-trip: gextract identical before and after conversion", {
    withr::defer(gtrack.rm("temp.rects_roundtrip", force = TRUE))

    # Create a RECTS track with known intervals and values
    intervs <- gintervals.2d(
        chroms1 = c(1, 1, 2, 2),
        starts1 = c(200, 400, 100, 300),
        ends1 = c(300, 500, 200, 400),
        chroms2 = c(1, 2, 2, 2),
        starts2 = c(300, 100, 200, 500),
        ends2 = c(400, 200, 300, 600)
    )
    values <- c(1.5, 2.7, 3.3, 4.1)
    gtrack.2d.create("temp.rects_roundtrip", "Test RECTS roundtrip", intervs, values)

    # Extract BEFORE conversion
    query_intervs <- gintervals.2d(chroms1 = c(1, 2), chroms2 = c(1, 2))
    before <- gextract("temp.rects_roundtrip", query_intervs)

    # Convert to indexed format
    gtrack.2d.convert_to_indexed("temp.rects_roundtrip", remove.old = TRUE)

    # Verify physical structure: track.idx and track.dat should exist
    trackdir <- .track_dir("temp.rects_roundtrip")
    expect_true(file.exists(file.path(trackdir, "track.idx")))
    expect_true(file.exists(file.path(trackdir, "track.dat")))

    # Extract AFTER conversion
    after <- gextract("temp.rects_roundtrip", query_intervs)

    # Verify identical results
    rownames(before) <- NULL
    rownames(after) <- NULL
    expect_equal(before, after, info = "RECTS round-trip: gextract should be identical after conversion")
})

test_that("RECTS track round-trip with multiple chromosome pairs", {
    withr::defer(gtrack.rm("temp.rects_multi", force = TRUE))

    # Create a RECTS track spanning many chromosome pairs
    intervs <- gintervals.2d(
        chroms1 = c(rep(1, 3), rep(1, 3), rep(2, 3)),
        starts1 = rep(c(100, 200, 300), 3),
        ends1 = rep(c(150, 250, 350), 3),
        chroms2 = c(rep(1, 3), rep(2, 3), rep(2, 3)),
        starts2 = rep(c(400, 500, 600), 3),
        ends2 = rep(c(450, 550, 650), 3)
    )
    values <- seq(1.0, 9.0, by = 1.0)
    gtrack.2d.create("temp.rects_multi", "Multi-pair RECTS", intervs, values)

    # Extract before
    before <- gextract("temp.rects_multi", gintervals.2d(chroms1 = c(1, 1, 2), chroms2 = c(1, 2, 2)))

    # Convert
    gtrack.2d.convert_to_indexed("temp.rects_multi", remove.old = TRUE)

    # Extract after
    after <- gextract("temp.rects_multi", gintervals.2d(chroms1 = c(1, 1, 2), chroms2 = c(1, 2, 2)))

    rownames(before) <- NULL
    rownames(after) <- NULL
    expect_equal(before, after, info = "Multi-pair RECTS round-trip should be identical")
})

# ============================================================================
# Category 2: Conversion round-trip for POINTS tracks
# ============================================================================

test_that("POINTS track round-trip: gextract identical before and after conversion", {
    withr::defer(gtrack.rm("temp.points_roundtrip", force = TRUE))

    # Create point data (end = start + 1 for both dimensions)
    intervs <- gintervals.2d(
        chroms1 = c(1, 1, 2, 2),
        starts1 = c(100, 200, 150, 350),
        ends1 = c(101, 201, 151, 351),
        chroms2 = c(1, 2, 2, 2),
        starts2 = c(500, 300, 400, 700),
        ends2 = c(501, 301, 401, 701)
    )
    values <- c(10.0, 20.0, 30.0, 40.0)

    # Create via gtrack.2d.import using a temp file (points tracks are created by import)
    tmp_file <- tempfile(fileext = ".txt")
    withr::defer(unlink(tmp_file))
    df <- data.frame(
        chrom1 = intervs$chrom1,
        start1 = intervs$start1,
        end1 = intervs$end1,
        chrom2 = intervs$chrom2,
        start2 = intervs$start2,
        end2 = intervs$end2,
        value = values
    )
    write.table(df, tmp_file, sep = "\t", row.names = FALSE, quote = FALSE)
    gtrack.2d.import("temp.points_roundtrip", "Test POINTS roundtrip", tmp_file)

    # Verify it is indeed a Points track
    info <- gtrack.info("temp.points_roundtrip")
    expect_equal(info$type, "points")

    # Extract BEFORE conversion
    query_intervs <- gintervals.2d(chroms1 = c(1, 2), chroms2 = c(1, 2))
    before <- gextract("temp.points_roundtrip", query_intervs)

    # Convert to indexed format
    gtrack.2d.convert_to_indexed("temp.points_roundtrip", remove.old = TRUE)

    # Verify physical structure
    trackdir <- .track_dir("temp.points_roundtrip")
    expect_true(file.exists(file.path(trackdir, "track.idx")))
    expect_true(file.exists(file.path(trackdir, "track.dat")))

    # Extract AFTER conversion
    after <- gextract("temp.points_roundtrip", query_intervs)

    rownames(before) <- NULL
    rownames(after) <- NULL
    expect_equal(before, after, info = "POINTS round-trip: gextract should be identical after conversion")
})

test_that("POINTS track round-trip with multi-pair data", {
    withr::defer(gtrack.rm("temp.points_multi", force = TRUE))

    # Create point data spanning multiple chromosome pairs
    tmp_file <- tempfile(fileext = ".txt")
    withr::defer(unlink(tmp_file))

    df <- data.frame(
        chrom1 = c("chr1", "chr1", "chr1", "chr2", "chr2"),
        start1 = c(100, 200, 300, 100, 200),
        end1 = c(101, 201, 301, 101, 201),
        chrom2 = c("chr1", "chr1", "chr2", "chr2", "chr2"),
        start2 = c(500, 600, 400, 300, 700),
        end2 = c(501, 601, 401, 301, 701),
        value = c(1.0, 2.0, 3.0, 4.0, 5.0)
    )
    write.table(df, tmp_file, sep = "\t", row.names = FALSE, quote = FALSE)
    gtrack.2d.import("temp.points_multi", "Multi-pair POINTS", tmp_file)

    # Extract before
    before <- gextract("temp.points_multi", gintervals.2d(chroms1 = c(1, 1, 2), chroms2 = c(1, 2, 2)))

    # Convert
    gtrack.2d.convert_to_indexed("temp.points_multi", remove.old = TRUE)

    # Extract after
    after <- gextract("temp.points_multi", gintervals.2d(chroms1 = c(1, 1, 2), chroms2 = c(1, 2, 2)))

    rownames(before) <- NULL
    rownames(after) <- NULL
    expect_equal(before, after, info = "Multi-pair POINTS round-trip should be identical")
})

# ============================================================================
# Category 3: Query parity
# ============================================================================

test_that("gextract with 2D intervals on indexed RECTS track", {
    withr::defer(gtrack.rm("temp.rects_query", force = TRUE))

    intervs <- gintervals.2d(
        chroms1 = c(rep(1, 4), rep(2, 4)),
        starts1 = rep(c(100, 200, 300, 400), 2),
        ends1 = rep(c(150, 250, 350, 450), 2),
        chroms2 = c(rep(1, 4), rep(2, 4)),
        starts2 = rep(c(500, 600, 700, 800), 2),
        ends2 = rep(c(550, 650, 750, 850), 2)
    )
    values <- seq(1.0, 8.0, by = 1.0)
    gtrack.2d.create("temp.rects_query", "RECTS query test", intervs, values)

    # Extract before
    query <- gintervals.2d(chroms1 = c(1, 2), chroms2 = c(1, 2))
    before <- gextract("temp.rects_query", query)

    # Convert
    gtrack.2d.convert_to_indexed("temp.rects_query", remove.old = TRUE)

    # gextract with same 2D intervals
    after <- gextract("temp.rects_query", query)
    rownames(before) <- NULL
    rownames(after) <- NULL
    expect_equal(before, after, info = "gextract with 2D intervals on indexed RECTS track")
})

test_that("gscreen with indexed 2D track", {
    withr::defer(gtrack.rm("temp.rects_screen", force = TRUE))

    intervs <- gintervals.2d(
        chroms1 = c(rep(1, 4), rep(2, 4)),
        starts1 = rep(c(100, 200, 300, 400), 2),
        ends1 = rep(c(150, 250, 350, 450), 2),
        chroms2 = c(rep(1, 4), rep(2, 4)),
        starts2 = rep(c(500, 600, 700, 800), 2),
        ends2 = rep(c(550, 650, 750, 850), 2)
    )
    values <- c(1.0, 5.0, 2.0, 8.0, 3.0, 7.0, 4.0, 6.0)
    gtrack.2d.create("temp.rects_screen", "RECTS screen test", intervs, values)

    query <- gintervals.2d(chroms1 = c(1, 2), chroms2 = c(1, 2))
    before <- gscreen("temp.rects_screen > 4", intervals = query)

    # Convert
    gtrack.2d.convert_to_indexed("temp.rects_screen", remove.old = TRUE)

    after <- gscreen("temp.rects_screen > 4", intervals = query)
    rownames(before) <- NULL
    rownames(after) <- NULL
    expect_equal(before, after, info = "gscreen with indexed 2D track should produce identical results")
    expect_true(nrow(after) > 0, info = "gscreen should return some results")
    expect_true(nrow(after) < 8, info = "gscreen should filter some intervals")
})

test_that("giterator.intervals with 2D track as data source", {
    withr::defer(gtrack.rm("temp.rects_giter", force = TRUE))

    intervs <- gintervals.2d(
        chroms1 = c(rep(1, 3), rep(2, 3)),
        starts1 = rep(c(100, 200, 300), 2),
        ends1 = rep(c(150, 250, 350), 2),
        chroms2 = c(rep(1, 3), rep(2, 3)),
        starts2 = rep(c(400, 500, 600), 2),
        ends2 = rep(c(450, 550, 650), 2)
    )
    values <- seq(1.0, 6.0, by = 1.0)
    gtrack.2d.create("temp.rects_giter", "RECTS giter test", intervs, values)

    # Use giterator.intervals before conversion
    before <- giterator.intervals(iterator = "temp.rects_giter")

    # Convert
    gtrack.2d.convert_to_indexed("temp.rects_giter", remove.old = TRUE)

    # Use giterator.intervals after conversion
    after <- giterator.intervals(iterator = "temp.rects_giter")

    rownames(before) <- NULL
    rownames(after) <- NULL
    expect_equal(before, after, info = "giterator.intervals with indexed 2D track")
})

test_that("band filtering on indexed 2D track", {
    withr::defer(gtrack.rm("temp.rects_band", force = TRUE))

    # Create track with intervals at various distances
    intervs <- gintervals.2d(
        chroms1 = rep(1, 6),
        starts1 = c(1000, 2000, 3000, 4000, 5000, 6000),
        ends1 = c(1100, 2100, 3100, 4100, 5100, 6100),
        chroms2 = rep(1, 6),
        starts2 = c(1500, 3000, 5000, 8000, 12000, 20000),
        ends2 = c(1600, 3100, 5100, 8100, 12100, 20100)
    )
    values <- c(10.0, 20.0, 30.0, 40.0, 50.0, 60.0)
    gtrack.2d.create("temp.rects_band", "RECTS band test", intervs, values)

    query <- gintervals.2d(chroms1 = 1, chroms2 = 1)
    band <- c(500, 5000)

    before <- gextract("temp.rects_band", query, band = band)

    # Convert
    gtrack.2d.convert_to_indexed("temp.rects_band", remove.old = TRUE)

    after <- gextract("temp.rects_band", query, band = band)
    rownames(before) <- NULL
    rownames(after) <- NULL
    expect_equal(before, after, info = "Band filtering on indexed 2D track should be identical")
})

test_that("track expression with arithmetic on indexed 2D track", {
    withr::defer(gtrack.rm("temp.rects_expr", force = TRUE))

    intervs <- gintervals.2d(
        chroms1 = c(1, 1, 2),
        starts1 = c(100, 200, 300),
        ends1 = c(150, 250, 350),
        chroms2 = c(1, 2, 2),
        starts2 = c(400, 500, 600),
        ends2 = c(450, 550, 650)
    )
    values <- c(10.0, 20.0, 30.0)
    gtrack.2d.create("temp.rects_expr", "RECTS expr test", intervs, values)

    query <- gintervals.2d(chroms1 = c(1, 2), chroms2 = c(1, 2))
    expr <- "temp.rects_expr * 2 + 5"

    before <- gextract(expr, query)

    # Convert
    gtrack.2d.convert_to_indexed("temp.rects_expr", remove.old = TRUE)

    after <- gextract(expr, query)
    rownames(before) <- NULL
    rownames(after) <- NULL
    expect_equal(before, after, info = "Arithmetic expression on indexed 2D track")
})

test_that("multiple 2D tracks in expression after conversion", {
    withr::defer({
        gtrack.rm("temp.rects_multi_a", force = TRUE)
        gtrack.rm("temp.rects_multi_b", force = TRUE)
    })

    intervs <- gintervals.2d(
        chroms1 = c(1, 1, 2),
        starts1 = c(100, 200, 300),
        ends1 = c(150, 250, 350),
        chroms2 = c(1, 2, 2),
        starts2 = c(400, 500, 600),
        ends2 = c(450, 550, 650)
    )
    gtrack.2d.create("temp.rects_multi_a", "Track A", intervs, c(1.0, 2.0, 3.0))
    gtrack.2d.create("temp.rects_multi_b", "Track B", intervs, c(10.0, 20.0, 30.0))

    query <- gintervals.2d(chroms1 = c(1, 2), chroms2 = c(1, 2))
    before <- gextract("temp.rects_multi_a", "temp.rects_multi_b",
        query,
        iterator = "temp.rects_multi_a"
    )

    # Convert both tracks
    gtrack.2d.convert_to_indexed("temp.rects_multi_a", remove.old = TRUE)
    gtrack.2d.convert_to_indexed("temp.rects_multi_b", remove.old = TRUE)

    after <- gextract("temp.rects_multi_a", "temp.rects_multi_b",
        query,
        iterator = "temp.rects_multi_a"
    )
    rownames(before) <- NULL
    rownames(after) <- NULL
    expect_equal(before, after, info = "Multiple indexed 2D tracks in expression")
})

# ============================================================================
# Category 4: gtrack.convert_to_indexed dispatches 2D
# ============================================================================

test_that("gtrack.convert_to_indexed dispatches to gtrack.2d.convert_to_indexed for 2D tracks", {
    withr::defer(gtrack.rm("temp.rects_dispatch", force = TRUE))

    intervs <- gintervals.2d(
        chroms1 = c(1, 2),
        starts1 = c(100, 200),
        ends1 = c(150, 250),
        chroms2 = c(1, 2),
        starts2 = c(300, 400),
        ends2 = c(350, 450)
    )
    gtrack.2d.create("temp.rects_dispatch", "Dispatch test", intervs, c(1.0, 2.0))

    query <- gintervals.2d(chroms1 = c(1, 2), chroms2 = c(1, 2))
    before <- gextract("temp.rects_dispatch", query)

    # Call the unified gtrack.convert_to_indexed (not the 2D-specific one)
    gtrack.convert_to_indexed("temp.rects_dispatch")

    # Verify it was converted (track.idx should exist)
    trackdir <- .track_dir("temp.rects_dispatch")
    expect_true(file.exists(file.path(trackdir, "track.idx")),
        info = "gtrack.convert_to_indexed should create track.idx for 2D tracks"
    )
    expect_true(file.exists(file.path(trackdir, "track.dat")),
        info = "gtrack.convert_to_indexed should create track.dat for 2D tracks"
    )

    # Verify data integrity
    after <- gextract("temp.rects_dispatch", query)
    rownames(before) <- NULL
    rownames(after) <- NULL
    expect_equal(before, after, info = "gtrack.convert_to_indexed dispatches correctly for 2D")
})

# ============================================================================
# Category 5: Auto-conversion in indexed databases
# ============================================================================

test_that("2D track is auto-converted when created in indexed database", {
    # Save current root to restore later
    local_db_state()

    # Create a minimal indexed database
    tmp_root <- withr::local_tempdir()
    test_fasta <- file.path(tmp_root, "genome.fasta")
    cat(">chr1\n", paste0(rep("A", 10000), collapse = ""), "\n",
        ">chr2\n", paste0(rep("C", 10000), collapse = ""), "\n",
        sep = "", file = test_fasta
    )

    db_path <- file.path(tmp_root, "testdb")
    withr::with_options(list(gmulticontig.indexed_format = TRUE), {
        gdb.create(groot = db_path, fasta = test_fasta, verbose = FALSE)
        gdb.init(db_path)

        # Verify the DB is indexed
        expect_true(.gdb.is_indexed())

        # Create a 2D track -- should be auto-converted
        intervs <- gintervals.2d(
            chroms1 = c(1, 2),
            starts1 = c(100, 200),
            ends1 = c(150, 250),
            chroms2 = c(1, 2),
            starts2 = c(300, 400),
            ends2 = c(350, 450)
        )
        gtrack.2d.create("test_auto_2d", "Auto-convert test", intervs, c(1.0, 2.0))
        withr::defer(gtrack.rm("test_auto_2d", force = TRUE))

        # Verify track.idx file exists (auto-converted)
        trackdir <- .track_dir("test_auto_2d")
        expect_true(file.exists(file.path(trackdir, "track.idx")),
            info = "2D track should be auto-converted in indexed database"
        )
        expect_true(file.exists(file.path(trackdir, "track.dat")),
            info = "track.dat should exist for auto-converted 2D track"
        )

        # Verify data is accessible
        result <- gextract("test_auto_2d", gintervals.2d(chroms1 = c(1, 2), chroms2 = c(1, 2)))
        expect_equal(nrow(result), 2)
    })
})

test_that("2D track via gtrack.2d.import is auto-converted in indexed database", {
    local_db_state()

    tmp_root <- withr::local_tempdir()
    test_fasta <- file.path(tmp_root, "genome.fasta")
    cat(">chr1\n", paste0(rep("A", 10000), collapse = ""), "\n",
        ">chr2\n", paste0(rep("C", 10000), collapse = ""), "\n",
        sep = "", file = test_fasta
    )

    db_path <- file.path(tmp_root, "testdb")
    withr::with_options(list(gmulticontig.indexed_format = TRUE), {
        gdb.create(groot = db_path, fasta = test_fasta, verbose = FALSE)
        gdb.init(db_path)
        expect_true(.gdb.is_indexed())

        # Create a points file for import
        tmp_file <- tempfile(fileext = ".txt")
        withr::defer(unlink(tmp_file))
        df <- data.frame(
            chrom1 = c("chr1", "chr2"),
            start1 = c(100, 200),
            end1 = c(101, 201),
            chrom2 = c("chr1", "chr2"),
            start2 = c(500, 600),
            end2 = c(501, 601),
            value = c(1.0, 2.0)
        )
        write.table(df, tmp_file, sep = "\t", row.names = FALSE, quote = FALSE)

        gtrack.2d.import("test_auto_import_2d", "Auto-import test", tmp_file)
        withr::defer(gtrack.rm("test_auto_import_2d", force = TRUE))

        # Verify auto-conversion
        trackdir <- .track_dir("test_auto_import_2d")
        expect_true(file.exists(file.path(trackdir, "track.idx")),
            info = "Imported 2D track should be auto-converted"
        )

        # Verify data
        result <- gextract("test_auto_import_2d", gintervals.2d(chroms1 = c(1, 2), chroms2 = c(1, 2)))
        expect_equal(nrow(result), 2)
    })
})

# ============================================================================
# Category 6: gdb.convert_to_indexed includes 2D tracks
# ============================================================================

test_that("gdb.convert_to_indexed converts 2D tracks when convert_tracks=TRUE", {
    local_db_state()

    tmp_root <- withr::local_tempdir()
    test_fasta <- file.path(tmp_root, "genome.fasta")
    cat(">genome\n", paste0(rep("A", 20000), collapse = ""), "\n",
        sep = "", file = test_fasta
    )

    db_path <- file.path(tmp_root, "testdb")
    withr::with_options(list(gmulticontig.indexed_format = FALSE), {
        gdb.create(groot = db_path, fasta = test_fasta, verbose = FALSE)
        gdb.init(db_path)

        # Create a 2D track in per-pair format (single chrom pair: genome-genome)
        intervs <- gintervals.2d(
            chroms1 = c("genome", "genome"),
            starts1 = c(100, 200),
            ends1 = c(150, 250),
            chroms2 = c("genome", "genome"),
            starts2 = c(300, 400),
            ends2 = c(350, 450)
        )
        gtrack.2d.create("test_batch_2d", "Batch convert test", intervs, c(1.0, 2.0))

        # Verify it is NOT indexed yet
        trackdir <- .track_dir("test_batch_2d")
        expect_false(file.exists(file.path(trackdir, "track.idx")))

        # Extract before conversion
        before <- gextract("test_batch_2d", gintervals.2d("genome"))

        # Run full DB conversion including tracks (keep old files so seq still works)
        gdb.convert_to_indexed(
            groot = db_path,
            convert_tracks = TRUE,
            remove_old_files = FALSE,
            force = TRUE,
            verbose = FALSE
        )

        # Reload after conversion
        gdb.init(db_path)

        # Verify 2D track was converted
        trackdir <- .track_dir("test_batch_2d")
        expect_true(file.exists(file.path(trackdir, "track.idx")),
            info = "gdb.convert_to_indexed should convert 2D tracks"
        )

        # Verify data integrity
        after <- gextract("test_batch_2d", gintervals.2d("genome"))
        rownames(before) <- NULL
        rownames(after) <- NULL
        expect_equal(before, after, info = "Data should be identical after batch conversion")
    })
})

# ============================================================================
# Category 7: Edge cases
# ============================================================================

test_that("convert already-indexed track with force=FALSE is a no-op", {
    withr::defer(gtrack.rm("temp.rects_noop", force = TRUE))

    intervs <- gintervals.2d(
        chroms1 = c(1),
        starts1 = c(100),
        ends1 = c(200),
        chroms2 = c(1),
        starts2 = c(300),
        ends2 = c(400)
    )
    gtrack.2d.create("temp.rects_noop", "No-op test", intervs, c(1.0))

    # Convert once
    gtrack.2d.convert_to_indexed("temp.rects_noop", remove.old = TRUE)
    trackdir <- .track_dir("temp.rects_noop")
    expect_true(file.exists(file.path(trackdir, "track.idx")))

    # Get modification time of track.idx
    idx_mtime <- file.mtime(file.path(trackdir, "track.idx"))

    # Brief pause to detect file modification
    Sys.sleep(1)

    # Convert again with force=FALSE -- should be no-op with a message
    expect_message(
        gtrack.2d.convert_to_indexed("temp.rects_noop", force = FALSE),
        "already in indexed format"
    )

    # Modification time should be unchanged
    idx_mtime_after <- file.mtime(file.path(trackdir, "track.idx"))
    expect_equal(idx_mtime, idx_mtime_after, info = "force=FALSE should not modify existing index")
})

test_that("convert already-indexed track with force=TRUE re-converts", {
    withr::defer(gtrack.rm("temp.rects_force", force = TRUE))

    intervs <- gintervals.2d(
        chroms1 = c(1),
        starts1 = c(100),
        ends1 = c(200),
        chroms2 = c(1),
        starts2 = c(300),
        ends2 = c(400)
    )
    gtrack.2d.create("temp.rects_force", "Force re-convert test", intervs, c(42.0))

    # Convert once
    gtrack.2d.convert_to_indexed("temp.rects_force", remove.old = TRUE)
    trackdir <- .track_dir("temp.rects_force")
    idx_mtime <- file.mtime(file.path(trackdir, "track.idx"))

    # Brief pause
    Sys.sleep(1)

    # Convert with force=TRUE -- should re-convert
    gtrack.2d.convert_to_indexed("temp.rects_force", force = TRUE)

    idx_mtime_after <- file.mtime(file.path(trackdir, "track.idx"))
    expect_true(idx_mtime_after > idx_mtime, info = "force=TRUE should re-create the index")

    # Verify data is still intact
    result <- gextract("temp.rects_force", gintervals.2d(chroms1 = 1, chroms2 = 1))
    expect_equal(nrow(result), 1)
    expect_equal(result$temp.rects_force, 42.0)
})

test_that("convert non-existent track errors", {
    expect_error(
        gtrack.2d.convert_to_indexed("nonexistent_track_xyz_2d"),
        "does not exist"
    )
})

test_that("convert 1D track via gtrack.2d.convert_to_indexed errors", {
    # test.fixedbin is a 1D dense track in the test database
    expect_error(
        gtrack.2d.convert_to_indexed("test.fixedbin"),
        "only 2D tracks"
    )
})

test_that("gtrack.info reports correct format after 2D conversion", {
    withr::defer(gtrack.rm("temp.rects_info", force = TRUE))

    intervs <- gintervals.2d(
        chroms1 = c(1, 2),
        starts1 = c(100, 200),
        ends1 = c(150, 250),
        chroms2 = c(1, 2),
        starts2 = c(300, 400),
        ends2 = c(350, 450)
    )
    gtrack.2d.create("temp.rects_info", "Info test", intervs, c(1.0, 2.0))

    # Before conversion
    info_before <- gtrack.info("temp.rects_info")
    expect_equal(info_before$type, "rectangles")

    # Convert
    gtrack.2d.convert_to_indexed("temp.rects_info", remove.old = TRUE)

    # After conversion: type should still be rectangles, format should be indexed
    info_after <- gtrack.info("temp.rects_info")
    expect_equal(info_after$type, "rectangles", info = "Track type should remain 'rectangles'")
    expect_equal(info_after$format, "indexed", info = "Format should report 'indexed' after conversion")
})

test_that("single chromosome pair track converts correctly", {
    withr::defer(gtrack.rm("temp.rects_single_pair", force = TRUE))

    # Track with data on only one chromosome pair
    intervs <- gintervals.2d(
        chroms1 = c(1, 1, 1),
        starts1 = c(100, 200, 300),
        ends1 = c(150, 250, 350),
        chroms2 = c(1, 1, 1),
        starts2 = c(500, 600, 700),
        ends2 = c(550, 650, 750)
    )
    values <- c(1.0, 2.0, 3.0)
    gtrack.2d.create("temp.rects_single_pair", "Single pair test", intervs, values)

    before <- gextract("temp.rects_single_pair", gintervals.2d(chroms1 = 1, chroms2 = 1))

    gtrack.2d.convert_to_indexed("temp.rects_single_pair", remove.old = TRUE)

    after <- gextract("temp.rects_single_pair", gintervals.2d(chroms1 = 1, chroms2 = 1))
    rownames(before) <- NULL
    rownames(after) <- NULL
    expect_equal(before, after, info = "Single pair track converts correctly")
})

test_that("querying empty chromosome pair on indexed 2D track returns no data", {
    withr::defer(gtrack.rm("temp.rects_empty_pair", force = TRUE))

    # Track with data only on chr1-chr1
    intervs <- gintervals.2d(
        chroms1 = c(1, 1),
        starts1 = c(100, 200),
        ends1 = c(150, 250),
        chroms2 = c(1, 1),
        starts2 = c(400, 500),
        ends2 = c(450, 550)
    )
    gtrack.2d.create("temp.rects_empty_pair", "Empty pair test", intervs, c(1.0, 2.0))

    # Get result for chr2-chr2 before conversion (should be empty/NaN)
    before_empty <- gextract("temp.rects_empty_pair", gintervals.2d(chroms1 = 2, chroms2 = 2))

    # Convert
    gtrack.2d.convert_to_indexed("temp.rects_empty_pair", remove.old = TRUE)

    # After conversion, querying an empty pair should give same (empty) result
    after_empty <- gextract("temp.rects_empty_pair", gintervals.2d(chroms1 = 2, chroms2 = 2))
    expect_equal(nrow(before_empty), nrow(after_empty),
        info = "Empty pair query should return same number of rows"
    )
})

# ============================================================================
# Category 8: Re-read parity with existing tracks
# ============================================================================

test_that("existing test.rects track can be converted and queried identically", {
    # The test database should have test.rects as a per-pair 2D track
    # Create a copy so we do not modify the original
    withr::defer(gtrack.rm("temp.rects_copy", force = TRUE))

    # Extract data from original track
    query <- gintervals.2d(chroms1 = c(2, 3), chroms2 = c(2, 4))
    original_data <- gextract("test.rects", query)

    # Create a copy by extracting and re-creating
    all_data <- gextract("test.rects", gintervals.2d(chroms1 = c(2, 3), chroms2 = c(2, 4)),
        iterator = "test.rects"
    )
    if (nrow(all_data) > 0) {
        track_col <- "test.rects"
        intervs_copy <- all_data[, c("chrom1", "start1", "end1", "chrom2", "start2", "end2")]
        values_copy <- all_data[[track_col]]
        gtrack.2d.create("temp.rects_copy", "Copy of test.rects", intervs_copy, values_copy)

        # Extract before conversion
        before <- gextract("temp.rects_copy", query)

        # Convert
        gtrack.2d.convert_to_indexed("temp.rects_copy", remove.old = TRUE)

        # Extract after conversion
        after <- gextract("temp.rects_copy", query)

        rownames(before) <- NULL
        rownames(after) <- NULL
        expect_equal(before, after,
            info = "Copied test.rects track should be identical after conversion"
        )
    }
})

test_that("2D track with iterator parity after conversion", {
    withr::defer(gtrack.rm("temp.rects_iter_parity", force = TRUE))

    intervs <- gintervals.2d(
        chroms1 = c(rep(1, 4), rep(2, 4)),
        starts1 = rep(c(100, 200, 300, 400), 2),
        ends1 = rep(c(150, 250, 350, 450), 2),
        chroms2 = c(rep(1, 4), rep(2, 4)),
        starts2 = rep(c(500, 600, 700, 800), 2),
        ends2 = rep(c(550, 650, 750, 850), 2)
    )
    values <- seq(1.0, 8.0, by = 1.0)
    gtrack.2d.create("temp.rects_iter_parity", "Iterator parity test", intervs, values)

    query <- gintervals.2d(chroms1 = c(1, 2), chroms2 = c(1, 2))

    # Use track as iterator before conversion
    before <- gextract("temp.rects_iter_parity", query, iterator = "temp.rects_iter_parity")

    # Convert
    gtrack.2d.convert_to_indexed("temp.rects_iter_parity", remove.old = TRUE)

    # Use track as iterator after conversion
    after <- gextract("temp.rects_iter_parity", query, iterator = "temp.rects_iter_parity")

    rownames(before) <- NULL
    rownames(after) <- NULL
    expect_equal(before, after, info = "Using indexed track as iterator should give identical results")
})

# ============================================================================
# Additional: Larger-scale and stress tests
# ============================================================================

test_that("RECTS track with many intervals per pair converts correctly", {
    withr::defer(gtrack.rm("temp.rects_large", force = TRUE))

    # Create a larger track: 50 intervals on each of 2 pairs
    n <- 50
    intervs <- gintervals.2d(
        chroms1 = c(rep(1, n), rep(2, n)),
        starts1 = c(seq(100, by = 200, length.out = n), seq(100, by = 200, length.out = n)),
        ends1 = c(seq(150, by = 200, length.out = n), seq(150, by = 200, length.out = n)),
        chroms2 = c(rep(1, n), rep(2, n)),
        starts2 = c(seq(10000, by = 200, length.out = n), seq(10000, by = 200, length.out = n)),
        ends2 = c(seq(10050, by = 200, length.out = n), seq(10050, by = 200, length.out = n))
    )
    values <- seq_len(2 * n) * 1.5
    gtrack.2d.create("temp.rects_large", "Large RECTS test", intervs, values)

    query <- gintervals.2d(chroms1 = c(1, 2), chroms2 = c(1, 2))
    before <- gextract("temp.rects_large", query)

    # Convert
    gtrack.2d.convert_to_indexed("temp.rects_large", remove.old = TRUE)

    after <- gextract("temp.rects_large", query)
    rownames(before) <- NULL
    rownames(after) <- NULL
    expect_equal(before, after, info = "Large RECTS track converts correctly")
    expect_equal(nrow(after), 2 * n)
})

test_that("RECTS track with all chromosome pair combinations", {
    withr::defer(gtrack.rm("temp.rects_allpairs", force = TRUE))

    # Create intervals for all 6 chromosome pair combinations (chr1, chr2, chrX)
    # chr1-chr1, chr1-chr2, chr1-chrX, chr2-chr2, chr2-chrX, chrX-chrX
    intervs <- gintervals.2d(
        chroms1 = c(1, 1, 1, 2, 2, "X"),
        starts1 = c(100, 200, 300, 400, 500, 600),
        ends1 = c(150, 250, 350, 450, 550, 650),
        chroms2 = c(1, 2, "X", 2, "X", "X"),
        starts2 = c(700, 800, 900, 1000, 1100, 1200),
        ends2 = c(750, 850, 950, 1050, 1150, 1250)
    )
    values <- c(1.1, 2.2, 3.3, 4.4, 5.5, 6.6)
    gtrack.2d.create("temp.rects_allpairs", "All pairs RECTS", intervs, values)

    # Extract from each pair individually before conversion
    pairs <- list(
        c(1, 1), c(1, 2), c(1, "X"),
        c(2, 2), c(2, "X"),
        c("X", "X")
    )

    before_results <- list()
    for (i in seq_along(pairs)) {
        p <- pairs[[i]]
        before_results[[i]] <- gextract(
            "temp.rects_allpairs",
            gintervals.2d(chroms1 = p[1], chroms2 = p[2])
        )
    }

    # Convert
    gtrack.2d.convert_to_indexed("temp.rects_allpairs", remove.old = TRUE)

    # Extract from each pair after conversion and compare
    for (i in seq_along(pairs)) {
        p <- pairs[[i]]
        after <- gextract(
            "temp.rects_allpairs",
            gintervals.2d(chroms1 = p[1], chroms2 = p[2])
        )

        rownames(before_results[[i]]) <- NULL
        rownames(after) <- NULL
        expect_equal(before_results[[i]], after,
            info = sprintf("Pair chr%s-chr%s should be identical after conversion", p[1], p[2])
        )
    }
})

test_that("gsummary works with indexed 2D track", {
    withr::defer(gtrack.rm("temp.rects_summary", force = TRUE))

    intervs <- gintervals.2d(
        chroms1 = c(rep(1, 4), rep(2, 4)),
        starts1 = rep(c(100, 200, 300, 400), 2),
        ends1 = rep(c(150, 250, 350, 450), 2),
        chroms2 = c(rep(1, 4), rep(2, 4)),
        starts2 = rep(c(500, 600, 700, 800), 2),
        ends2 = rep(c(550, 650, 750, 850), 2)
    )
    values <- c(1, 2, 3, 4, 5, 6, 7, 8)
    gtrack.2d.create("temp.rects_summary", "Summary test", intervs, values)

    query <- gintervals.2d(chroms1 = c(1, 2), chroms2 = c(1, 2))
    before <- gsummary("temp.rects_summary", query)

    # Convert
    gtrack.2d.convert_to_indexed("temp.rects_summary", remove.old = TRUE)

    after <- gsummary("temp.rects_summary", query)
    expect_equal(before, after, info = "gsummary should be identical after conversion")
})

test_that("mixed indexed and non-indexed 2D tracks work in same expression", {
    withr::defer({
        gtrack.rm("temp.rects_indexed_mix", force = TRUE)
        gtrack.rm("temp.rects_plain_mix", force = TRUE)
    })

    intervs <- gintervals.2d(
        chroms1 = c(1, 1, 2),
        starts1 = c(100, 200, 300),
        ends1 = c(150, 250, 350),
        chroms2 = c(1, 2, 2),
        starts2 = c(400, 500, 600),
        ends2 = c(450, 550, 650)
    )
    gtrack.2d.create("temp.rects_indexed_mix", "Indexed track", intervs, c(1.0, 2.0, 3.0))
    gtrack.2d.create("temp.rects_plain_mix", "Plain track", intervs, c(10.0, 20.0, 30.0))

    query <- gintervals.2d(chroms1 = c(1, 2), chroms2 = c(1, 2))

    # Extract both before any conversion
    before <- gextract("temp.rects_indexed_mix", "temp.rects_plain_mix",
        query,
        iterator = "temp.rects_indexed_mix"
    )

    # Convert only one
    gtrack.2d.convert_to_indexed("temp.rects_indexed_mix", remove.old = TRUE)

    # Mix indexed and non-indexed in same expression
    after <- gextract("temp.rects_indexed_mix", "temp.rects_plain_mix",
        query,
        iterator = "temp.rects_indexed_mix"
    )

    rownames(before) <- NULL
    rownames(after) <- NULL
    expect_equal(before, after,
        info = "Mixed indexed and non-indexed 2D tracks should produce identical results"
    )
})

test_that("virtual track on indexed 2D track works correctly", {
    withr::defer({
        gtrack.rm("temp.rects_vt", force = TRUE)
        try(gvtrack.rm("temp_vt_2d"), silent = TRUE)
    })

    intervs <- gintervals.2d(
        chroms1 = c(1, 1, 2),
        starts1 = c(100, 200, 300),
        ends1 = c(150, 250, 350),
        chroms2 = c(1, 2, 2),
        starts2 = c(400, 500, 600),
        ends2 = c(450, 550, 650)
    )
    gtrack.2d.create("temp.rects_vt", "VTrack test", intervs, c(10.0, 20.0, 30.0))

    gvtrack.create("temp_vt_2d", "temp.rects_vt")
    query <- gintervals.2d(chroms1 = c(1, 2), chroms2 = c(1, 2))

    before <- gextract("temp_vt_2d", query, iterator = "temp.rects_vt")

    # Convert
    gtrack.2d.convert_to_indexed("temp.rects_vt", remove.old = TRUE)

    after <- gextract("temp_vt_2d", query, iterator = "temp.rects_vt")

    rownames(before) <- NULL
    rownames(after) <- NULL
    expect_equal(before, after, info = "Virtual track on indexed 2D track")
})

test_that("gextract with colnames parameter works on indexed 2D track", {
    withr::defer(gtrack.rm("temp.rects_colnames", force = TRUE))

    intervs <- gintervals.2d(
        chroms1 = c(1, 2),
        starts1 = c(100, 200),
        ends1 = c(150, 250),
        chroms2 = c(1, 2),
        starts2 = c(300, 400),
        ends2 = c(350, 450)
    )
    gtrack.2d.create("temp.rects_colnames", "Colnames test", intervs, c(1.0, 2.0))

    query <- gintervals.2d(chroms1 = c(1, 2), chroms2 = c(1, 2))
    before <- gextract("temp.rects_colnames", query, colnames = "my_val")

    # Convert
    gtrack.2d.convert_to_indexed("temp.rects_colnames", remove.old = TRUE)

    after <- gextract("temp.rects_colnames", query, colnames = "my_val")
    rownames(before) <- NULL
    rownames(after) <- NULL
    expect_equal(before, after, info = "colnames parameter with indexed 2D track")
    expect_true("my_val" %in% colnames(after))
})
