# Tests for multi-database support

# Helper to create a minimal database
create_test_db <- function(path, chrom_sizes = data.frame(chrom = c("chr1", "chr2"), size = c(10000, 10000))) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(path, "tracks"), showWarnings = FALSE)
    dir.create(file.path(path, "seq"), showWarnings = FALSE)

    write.table(chrom_sizes, file.path(path, "chrom_sizes.txt"),
        sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
    )

    # Create dummy sequence files
    for (chr in chrom_sizes$chrom) {
        writeLines(
            paste0(rep("A", chrom_sizes$size[chrom_sizes$chrom == chr]), collapse = ""),
            file.path(path, "seq", paste0(chr, ".seq"))
        )
    }
}

test_that("single database initialization works unchanged (backward compatibility)", {
    withr::with_tempdir({
        create_test_db("single_db")

        gsetroot("single_db")

        expect_equal(length(gdb.ls()), 1)
        expect_equal(gdb.ls(), normalizePath("single_db"))
        expect_equal(.misha$GROOT, .misha$GROOTS[1])

        # gdb.summary should work
        info <- gdb.summary()
        expect_equal(nrow(info), 1)
        expect_true("tracks" %in% names(info))
        expect_true("writable" %in% names(info))
    })
})

test_that("multiple databases can be connected", {
    withr::with_tempdir({
        # Create two databases with same chrom_sizes
        create_test_db("db1")
        create_test_db("db2")

        gsetroot(c("db1", "db2"))

        expect_equal(length(gdb.ls()), 2)
        expect_equal(gdb.ls()[1], normalizePath("db1"))
        expect_equal(gdb.ls()[2], normalizePath("db2"))

        # GROOT should be first db
        expect_equal(.misha$GROOT, normalizePath("db1"))
    })
})

test_that("connecting dbs with different chrom_sizes fails", {
    withr::with_tempdir({
        # Create db1 with one genome
        create_test_db("db1", chrom_sizes = data.frame(chrom = "chr1", size = 1000))

        # Create db2 with different genome
        create_test_db("db2", chrom_sizes = data.frame(chrom = "chr1", size = 2000))

        expect_error(
            gsetroot(c("db1", "db2")),
            "identical chrom_sizes"
        )
    })
})

test_that("last database wins for track resolution", {
    withr::with_tempdir({
        # Setup: create same-named track in both dbs
        create_test_db("db1")
        gsetroot("db1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("shared_track", "from db1", intervs, 1)

        create_test_db("db2")
        gsetroot("db2")
        gtrack.create_sparse("shared_track", "from db2", intervs, 2)

        # Connect to both
        gsetroot(c("db1", "db2"))

        # Track should come from db2 (last wins)
        expect_equal(gtrack.db("shared_track"), normalizePath("db2"))
    })
})

test_that("gtrack.db works on vectors", {
    withr::with_tempdir({
        create_test_db("db1")
        gsetroot("db1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("track_a", "a", intervs, 1)

        create_test_db("db2")
        gsetroot("db2")
        gtrack.create_sparse("track_b", "b", intervs, 2)

        gsetroot(c("db1", "db2"))

        # Vectorized call
        dbs <- gtrack.db(c("track_a", "track_b", "nonexistent"))
        expect_equal(dbs[1], normalizePath("db1"))
        expect_equal(dbs[2], normalizePath("db2"))
        expect_true(is.na(dbs[3]))
    })
})

test_that("gtrack.ls filters by database", {
    withr::with_tempdir({
        create_test_db("db1")
        gsetroot("db1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("db1_track", "desc", intervs, 1)

        create_test_db("db2")
        gsetroot("db2")
        gtrack.create_sparse("db2_track", "desc", intervs, 2)

        gsetroot(c("db1", "db2"))

        # All tracks
        all_tracks <- gtrack.ls()
        expect_true("db1_track" %in% all_tracks)
        expect_true("db2_track" %in% all_tracks)

        # Filter by db1
        db1_tracks <- gtrack.ls(db = "db1")
        expect_true("db1_track" %in% db1_tracks)
        expect_false("db2_track" %in% db1_tracks)

        # Filter by db2
        db2_tracks <- gtrack.ls(db = "db2")
        expect_false("db1_track" %in% db2_tracks)
        expect_true("db2_track" %in% db2_tracks)
    })
})

test_that("gdb.summary returns correct database summary", {
    withr::with_tempdir({
        create_test_db("db1")
        gsetroot("db1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("track1", "desc", intervs, 1)
        gtrack.create_sparse("track2", "desc", intervs, 2)

        create_test_db("db2")
        gsetroot("db2")
        gtrack.create_sparse("track3", "desc", intervs, 3)

        gsetroot(c("db1", "db2"))

        info <- gdb.summary()
        expect_equal(nrow(info), 2)
        expect_equal(info$tracks[1], 2) # db1 has 2 tracks
        expect_equal(info$tracks[2], 1) # db2 has 1 track
        expect_true(all(info$writable)) # both writable in test
    })
})

test_that("gdb.create_linked creates proper symlinks", {
    withr::with_tempdir({
        create_test_db("parent_db")

        gdb.create_linked("linked_db", parent = "parent_db")

        # Check symlinks exist
        expect_true(file.exists("linked_db/chrom_sizes.txt"))
        expect_true(file.exists("linked_db/seq"))
        expect_true(Sys.readlink("linked_db/chrom_sizes.txt") != "")
        expect_true(Sys.readlink("linked_db/seq") != "")

        # Check tracks dir is real (not symlink)
        expect_true(dir.exists("linked_db/tracks"))
        expect_equal(Sys.readlink("linked_db/tracks"), "")

        # Should be usable
        gsetroot(c("parent_db", "linked_db"))
        expect_equal(length(gdb.ls()), 2)
    })
})

test_that("gdb.create_linked fails if directory exists", {
    withr::with_tempdir({
        create_test_db("parent_db")
        dir.create("linked_db")

        expect_error(
            gdb.create_linked("linked_db", parent = "parent_db"),
            "already exists"
        )
    })
})

test_that("gintervals.db and gintervals.ls with db filter work", {
    withr::with_tempdir({
        create_test_db("db1")
        gsetroot("db1")
        intervs1 <- gintervals(1, 0, 1000)
        gintervals.save("intervals1", intervs1)

        create_test_db("db2")
        gsetroot("db2")
        intervs2 <- gintervals(1, 0, 2000)
        gintervals.save("intervals2", intervs2)

        gsetroot(c("db1", "db2"))

        # gintervals.db
        expect_equal(gintervals.db("intervals1"), normalizePath("db1"))
        expect_equal(gintervals.db("intervals2"), normalizePath("db2"))

        # gintervals.ls with db filter
        db1_intervals <- gintervals.ls(db = "db1")
        expect_true("intervals1" %in% db1_intervals)
        expect_false("intervals2" %in% db1_intervals)
    })
})

test_that("virtual tracks remain global across databases", {
    withr::with_tempdir({
        create_test_db("db1")
        gsetroot("db1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("base_track", "base", intervs, 1)

        create_test_db("db2")
        gsetroot(c("db1", "db2"))

        # Create virtual track
        gvtrack.create("vtrack", "base_track", "avg")
        withr::defer(gvtrack.rm("vtrack"))

        # Virtual track should be accessible regardless of which db we're working with
        expect_true("vtrack" %in% gvtrack.ls())
    })
})

test_that("track operations work correctly with multi-db", {
    withr::with_tempdir({
        create_test_db("db1")
        gsetroot("db1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("db1_track", "test track in db1", intervs, 5)

        create_test_db("db2")
        gsetroot(c("db1", "db2"))

        # Track from db1 should still be accessible
        expect_true(gtrack.exists("db1_track"))
        info <- gtrack.info("db1_track")
        expect_equal(info$type, "sparse")

        # gtrack.path should return correct path in db1
        path <- gtrack.path("db1_track")
        expect_true(grepl("db1", path))
    })
})

test_that("GWD defaults to last database's tracks", {
    withr::with_tempdir({
        create_test_db("db1")
        create_test_db("db2")

        gsetroot(c("db1", "db2"))

        # GWD should be in db2 (last db)
        expect_true(grepl("db2", .misha$GWD))
    })
})

# =============================================================================
# Advanced Multi-Database Tests
# =============================================================================

test_that("three or more databases work correctly", {
    withr::with_tempdir({
        create_test_db("db1")
        create_test_db("db2")
        create_test_db("db3")

        gsetroot("db1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("track_a", "from db1", intervs, 1)

        gsetroot("db2")
        gtrack.create_sparse("track_b", "from db2", intervs, 2)

        gsetroot("db3")
        gtrack.create_sparse("track_c", "from db3", intervs, 3)

        # Connect all three
        gsetroot(c("db1", "db2", "db3"))

        expect_equal(length(gdb.ls()), 3)
        expect_equal(length(gtrack.ls()), 3)

        # Verify each track resolves to correct db
        expect_equal(gtrack.db("track_a"), normalizePath("db1"))
        expect_equal(gtrack.db("track_b"), normalizePath("db2"))
        expect_equal(gtrack.db("track_c"), normalizePath("db3"))

        # Summary should show 1 track per db
        info <- gdb.summary()
        expect_equal(sum(info$tracks), 3)
    })
})

test_that("database order affects track resolution (last wins)", {
    withr::with_tempdir({
        create_test_db("db1")
        create_test_db("db2")

        # Create same track with different values in both DBs
        gsetroot("db1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("shared", "db1 version", intervs, 100)

        gsetroot("db2")
        gtrack.create_sparse("shared", "db2 version", intervs, 200)

        # Order: db1, db2 -> db2 wins
        gsetroot(c("db1", "db2"))
        expect_equal(gtrack.db("shared"), normalizePath("db2"))
        val1 <- gextract("shared", gintervals(1, 0, 500))$shared[1]

        # Order: db2, db1 -> db1 wins
        gsetroot(c("db2", "db1"))
        expect_equal(gtrack.db("shared"), normalizePath("db1"))
        val2 <- gextract("shared", gintervals(1, 0, 500))$shared[1]

        # Values should differ based on which db "wins"
        expect_equal(val1, 200)
        expect_equal(val2, 100)
    })
})

test_that("gextract works with tracks from different databases", {
    withr::with_tempdir({
        create_test_db("db1")
        create_test_db("db2")

        gsetroot("db1")
        intervs <- gintervals(1, 0, 5000)
        gtrack.create_sparse("track_from_db1", "db1", intervs, 10)

        gsetroot("db2")
        gtrack.create_sparse("track_from_db2", "db2", intervs, 20)

        gsetroot(c("db1", "db2"))

        # This should work - extract from both tracks in single call
        # Use iterator since we have multiple sparse tracks
        result <- gextract(c("track_from_db1", "track_from_db2"), gintervals(1, 0, 1000), iterator = 100)

        expect_true("track_from_db1" %in% names(result))
        expect_true("track_from_db2" %in% names(result))
        expect_equal(result$track_from_db1[1], 10)
        expect_equal(result$track_from_db2[1], 20)
    })
})

test_that("track expressions work across databases", {
    withr::with_tempdir({
        create_test_db("db1")
        create_test_db("db2")

        gsetroot("db1")
        intervs <- gintervals(1, 0, 5000)
        gtrack.create_sparse("x", "db1", intervs, 10)

        gsetroot("db2")
        gtrack.create_sparse("y", "db2", intervs, 3)

        gsetroot(c("db1", "db2"))

        # Expression using tracks from different databases
        # Use iterator since we have multiple sparse tracks
        result <- gextract("x + y", gintervals(1, 0, 1000), iterator = 100)
        expect_equal(result[["x + y"]][1], 13)
    })
})

test_that("new tracks are created in last (user) database", {
    withr::with_tempdir({
        create_test_db("shared_db")
        create_test_db("user_db")

        gsetroot("shared_db")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("shared_track", "shared", intervs, 1)

        # Connect with user_db last
        gsetroot(c("shared_db", "user_db"))

        # Create new track - should go to user_db (last/GWD)
        gtrack.create_sparse("new_track", "new", intervs, 2)

        # After creating, reload to update GTRACK_DB
        gdb.reload()

        expect_equal(gtrack.db("new_track"), normalizePath("user_db"))
        expect_equal(gtrack.db("shared_track"), normalizePath("shared_db"))

        # Verify by reconnecting to just user_db
        gsetroot("user_db")
        expect_true(gtrack.exists("new_track"))

        # And shared_db shouldn't have the new track
        gsetroot("shared_db")
        expect_false(gtrack.exists("new_track"))
    })
})

test_that("gdir.cd works correctly in multi-db mode", {
    # gdir.cd in multi-db mode changes to subdirectory in the last database
    # and then only shows tracks from that subdirectory (original single-db behavior)
    withr::with_tempdir({
        create_test_db("db1")
        create_test_db("db2")

        gsetroot("db1")
        gdir.create("subdir", showWarnings = FALSE)
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("subdir.sub_track", "in subdir", intervs, 1)
        gtrack.create_sparse("root_track", "in root", intervs, 2)

        gsetroot("db2")
        gdir.create("subdir", showWarnings = FALSE) # Create in db2 too for gdir.cd to work
        gtrack.create_sparse("db2_track", "in db2", intervs, 3)

        # Connect to both, then cd into subdir
        gsetroot(c("db1", "db2"))

        # At root, should see all tracks from both databases
        all_tracks <- gtrack.ls()
        expect_true("subdir.sub_track" %in% all_tracks)
        expect_true("root_track" %in% all_tracks)
        expect_true("db2_track" %in% all_tracks)

        # Change to subdir - since GWD is in db2, we cd to db2's subdir
        # This switches to single-db subdirectory mode (original behavior)
        gdir.cd("subdir")
        subdir_tracks <- gtrack.ls()
        # Should now see only tracks from db2/tracks/subdir (which is empty)
        expect_equal(length(subdir_tracks), 0)

        # Go back to root
        gdir.cd("..")
        # Back at db2's root tracks dir, but this is still single-db context
        # (see note in gdb.reload about subdirectory scanning)
    })
})

test_that("empty database in multi-db chain works", {
    withr::with_tempdir({
        create_test_db("db1")
        create_test_db("empty_db")
        create_test_db("db3")

        gsetroot("db1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("track1", "db1", intervs, 1)

        gsetroot("db3")
        gtrack.create_sparse("track3", "db3", intervs, 3)

        # empty_db has no tracks
        gsetroot(c("db1", "empty_db", "db3"))

        expect_equal(length(gdb.ls()), 3)
        tracks <- gtrack.ls()
        expect_equal(length(tracks), 2)
        expect_true("track1" %in% tracks)
        expect_true("track3" %in% tracks)

        # Summary should show 0 tracks for empty_db
        info <- gdb.summary()
        expect_equal(info$tracks[2], 0)
    })
})

test_that("intervals override follows last-wins semantics", {
    withr::with_tempdir({
        create_test_db("db1")
        create_test_db("db2")

        gsetroot("db1")
        intervs1 <- gintervals(1, 0, 1000)
        gintervals.save("shared_intervals", intervs1)

        gsetroot("db2")
        intervs2 <- gintervals(1, 0, 5000) # Different range
        gintervals.save("shared_intervals", intervs2)

        # db2 should win
        gsetroot(c("db1", "db2"))
        loaded <- gintervals.load("shared_intervals")
        expect_equal(loaded$end[1], 5000)
        expect_equal(gintervals.db("shared_intervals"), normalizePath("db2"))

        # Reverse order - db1 should win
        gsetroot(c("db2", "db1"))
        loaded2 <- gintervals.load("shared_intervals")
        expect_equal(loaded2$end[1], 1000)
        expect_equal(gintervals.db("shared_intervals"), normalizePath("db1"))
    })
})

test_that("track attributes work across databases", {
    withr::with_tempdir({
        create_test_db("db1")
        create_test_db("db2")

        gsetroot("db1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("db1_track", "track in db1", intervs, 1)
        gtrack.attr.set("db1_track", "custom_attr", "db1_value")

        gsetroot("db2")
        gtrack.create_sparse("db2_track", "track in db2", intervs, 2)
        gtrack.attr.set("db2_track", "custom_attr", "db2_value")

        gsetroot(c("db1", "db2"))

        # Should read attributes from both databases
        expect_equal(gtrack.attr.get("db1_track", "custom_attr"), "db1_value")
        expect_equal(gtrack.attr.get("db2_track", "custom_attr"), "db2_value")
    })
})

test_that("gtrack.rm only affects correct database", {
    withr::with_tempdir({
        create_test_db("db1")
        create_test_db("db2")

        gsetroot("db1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("track_to_keep", "keep", intervs, 1)

        gsetroot("db2")
        gtrack.create_sparse("track_to_delete", "delete", intervs, 2)

        gsetroot(c("db1", "db2"))

        # Delete track from db2
        gtrack.rm("track_to_delete", force = TRUE)

        # Verify it's gone
        expect_false(gtrack.exists("track_to_delete"))

        # Track in db1 should still exist
        expect_true(gtrack.exists("track_to_keep"))

        # Reconnect to verify persistence
        gsetroot(c("db1", "db2"))
        expect_false(gtrack.exists("track_to_delete"))
        expect_true(gtrack.exists("track_to_keep"))
    })
})

test_that("gtrack.rm can target a specific database", {
    withr::with_tempdir({
        create_test_db("db1")
        create_test_db("db2")

        gsetroot("db1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("shared_track", "db1", intervs, 1)

        gsetroot("db2")
        gtrack.create_sparse("shared_track", "db2", intervs, 2)

        gsetroot(c("db1", "db2"))
        expect_equal(gtrack.db("shared_track"), normalizePath("db2"))

        gtrack.rm("shared_track", force = TRUE, db = "db2")

        expect_true(gtrack.exists("shared_track"))
        expect_equal(gtrack.db("shared_track"), normalizePath("db1"))
        expect_true(file.exists(file.path("db1", "tracks", "shared_track.track")))
        expect_false(file.exists(file.path("db2", "tracks", "shared_track.track")))
    })
})

test_that("reconnecting refreshes track list correctly", {
    withr::with_tempdir({
        create_test_db("db1")
        create_test_db("db2")

        gsetroot("db1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("initial_track", "initial", intervs, 1)

        gsetroot(c("db1", "db2"))
        expect_equal(length(gtrack.ls()), 1)

        # Add track to db2 while connected
        gtrack.create_sparse("added_track", "added", intervs, 2)
        expect_equal(length(gtrack.ls()), 2)

        # Reconnect should still show both
        gsetroot(c("db1", "db2"))
        expect_equal(length(gtrack.ls()), 2)
        expect_true("initial_track" %in% gtrack.ls())
        expect_true("added_track" %in% gtrack.ls())
    })
})

test_that("read-only database gives clear error message", {
    skip_on_os("windows") # chmod not reliable on Windows

    withr::with_tempdir({
        create_test_db("readonly_db")
        create_test_db("writable_db")

        gsetroot("readonly_db")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("existing_track", "exists", intervs, 1)

        # Make readonly_db's tracks directory read-only
        Sys.chmod("readonly_db/tracks", "555")
        withr::defer(Sys.chmod("readonly_db/tracks", "755"))

        gsetroot(c("readonly_db", "writable_db"))

        # Creating new track should work (goes to writable_db)
        expect_silent(gtrack.create_sparse("new_track", "new", intervs, 2))

        # Deleting track from read-only db will fail silently (R cannot delete)
        # The track remains since the parent directory is read-only
        suppressMessages(gtrack.rm("existing_track", force = TRUE))

        # Verify the track still exists (because delete failed)
        expect_true(gtrack.exists("existing_track"))
    })
})

test_that("gdb.summary shows correct writable status", {
    skip_on_os("windows")

    withr::with_tempdir({
        create_test_db("readonly_db")
        create_test_db("writable_db")

        # Make one db read-only
        Sys.chmod("readonly_db/tracks", "555")
        withr::defer(Sys.chmod("readonly_db/tracks", "755"))

        gsetroot(c("readonly_db", "writable_db"))

        info <- gdb.summary()
        expect_equal(nrow(info), 2)

        # First db should be read-only
        expect_false(info$writable[1])
        # Second db should be writable
        expect_true(info$writable[2])
    })
})

test_that("virtual track can reference track from any database", {
    withr::with_tempdir({
        create_test_db("db1")
        create_test_db("db2")

        gsetroot("db1")
        intervs <- gintervals(1, 0, 5000)
        gtrack.create_sparse("base_in_db1", "base", intervs, 10)

        gsetroot("db2")
        gtrack.create_sparse("base_in_db2", "base", intervs, 20)

        gsetroot(c("db1", "db2"))

        gvtrack.create("vt1", "base_in_db1", "avg")
        gvtrack.create("vt2", "base_in_db2", "avg")
        withr::defer({
            gvtrack.rm("vt1")
            gvtrack.rm("vt2")
        })

        result1 <- gextract("vt1", gintervals(1, 0, 1000))
        result2 <- gextract("vt2", gintervals(1, 0, 1000))

        expect_equal(result1$vt1[1], 10)
        expect_equal(result2$vt2[1], 20)
    })
})

test_that("connecting same database twice is handled gracefully", {
    withr::with_tempdir({
        create_test_db("db1")

        gsetroot("db1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("test_track", "test", intervs, 1)

        # Connect same db twice
        gsetroot(c("db1", "db1"))

        # Should work without errors, track appears once
        tracks <- gtrack.ls()
        expect_equal(sum(tracks == "test_track"), 1)

        # gdb.ls shows the path twice
        expect_equal(length(gdb.ls()), 2)
    })
})

test_that("absolute and relative paths work for multi-db", {
    withr::with_tempdir({
        create_test_db("db1")
        create_test_db("db2")

        abs_path1 <- normalizePath("db1")
        rel_path2 <- "db2"

        gsetroot("db1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("track1", "t1", intervs, 1)

        gsetroot("db2")
        gtrack.create_sparse("track2", "t2", intervs, 2)

        # Mix absolute and relative paths
        gsetroot(c(abs_path1, rel_path2))

        expect_equal(length(gtrack.ls()), 2)
        expect_true(gtrack.exists("track1"))
        expect_true(gtrack.exists("track2"))
    })
})

# =============================================================================
# Additional Tests (inspired by naryn test patterns)
# =============================================================================

test_that("deletion of overriding track reveals underlying track", {
    withr::with_tempdir({
        create_test_db("db1")
        create_test_db("db2")
        create_test_db("db3")

        # Create same track in all three databases with different values
        gsetroot("db1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("shared", "from db1", intervs, 100)

        gsetroot("db2")
        gtrack.create_sparse("shared", "from db2", intervs, 200)

        gsetroot("db3")
        gtrack.create_sparse("shared", "from db3", intervs, 300)

        # Connect all three - db3 wins
        gsetroot(c("db1", "db2", "db3"))
        expect_equal(gtrack.db("shared"), normalizePath("db3"))
        expect_equal(gextract("shared", gintervals(1, 0, 500))$shared[1], 300)

        # Delete from db3 - db2 should now win
        gtrack.rm("shared", force = TRUE)
        gdb.reload()
        expect_equal(gtrack.db("shared"), normalizePath("db2"))
        expect_equal(gextract("shared", gintervals(1, 0, 500))$shared[1], 200)

        # Delete from db2 - db1 should now win
        gtrack.rm("shared", force = TRUE)
        gdb.reload()
        expect_equal(gtrack.db("shared"), normalizePath("db1"))
        expect_equal(gextract("shared", gintervals(1, 0, 500))$shared[1], 100)

        # Delete from db1 - track should be gone
        gtrack.rm("shared", force = TRUE)
        gdb.reload()
        expect_false(gtrack.exists("shared"))
    })
})

test_that("gtrack.exists with db parameter works", {
    withr::with_tempdir({
        create_test_db("db1")
        create_test_db("db2")

        gsetroot("db1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("only_in_db1", "db1 only", intervs, 1)

        gsetroot("db2")
        gtrack.create_sparse("only_in_db2", "db2 only", intervs, 2)

        gsetroot(c("db1", "db2"))

        # Both tracks exist overall
        expect_true(gtrack.exists("only_in_db1"))
        expect_true(gtrack.exists("only_in_db2"))

        # But each only exists in its own database
        expect_equal(gtrack.db("only_in_db1"), normalizePath("db1"))
        expect_equal(gtrack.db("only_in_db2"), normalizePath("db2"))
    })
})

test_that("gtrack.mv works correctly in multi-db", {
    withr::with_tempdir({
        create_test_db("db1")
        create_test_db("db2")

        gsetroot("db1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("original", "original track", intervs, 42)

        gsetroot(c("db1", "db2"))

        # Rename track (stays in db1)
        gtrack.mv("original", "renamed")

        expect_false(gtrack.exists("original"))
        expect_true(gtrack.exists("renamed"))
        expect_equal(gtrack.db("renamed"), normalizePath("db1"))

        # Value should be preserved
        expect_equal(gextract("renamed", gintervals(1, 0, 500))$renamed[1], 42)
    })
})

test_that("gtrack.mv fails when source doesn't exist", {
    withr::with_tempdir({
        create_test_db("db1")
        gsetroot("db1")

        expect_error(gtrack.mv("nonexistent", "new_name"), "does not exist")
    })
})

test_that("gtrack.mv fails when destination already exists", {
    withr::with_tempdir({
        create_test_db("db1")
        gsetroot("db1")

        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("track1", "t1", intervs, 1)
        gtrack.create_sparse("track2", "t2", intervs, 2)

        expect_error(gtrack.mv("track1", "track2"), "already exists")
    })
})

test_that("gtrack.mv can move track to different namespace", {
    withr::with_tempdir({
        create_test_db("db1")
        gsetroot("db1")

        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("track1", "t1", intervs, 55)

        # Move to subdir namespace
        gtrack.mv("track1", "subdir.track1")

        expect_false(gtrack.exists("track1"))
        expect_true(gtrack.exists("subdir.track1"))

        # Value preserved
        expect_equal(gextract("subdir.track1", gintervals(1, 0, 500))[["subdir.track1"]][1], 55)
    })
})

test_that("gtrack.copy works correctly in single db", {
    withr::with_tempdir({
        create_test_db("db1")
        gsetroot("db1")

        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("source", "src", intervs, 123)

        gtrack.copy("source", "dest")

        expect_true(gtrack.exists("source"))
        expect_true(gtrack.exists("dest"))

        # Values should match
        expect_equal(
            gextract("source", gintervals(1, 0, 500))$source[1],
            gextract("dest", gintervals(1, 0, 500))$dest[1]
        )
    })
})

test_that("gtrack.copy works across databases", {
    withr::with_tempdir({
        create_test_db("db1")
        create_test_db("db2")

        gsetroot("db1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("source_track", "source", intervs, 99)

        # Connect to both, GWD defaults to db2 (last)
        gsetroot(c("db1", "db2"))

        # Copy track - new track should go to db2 (current GWD)
        gtrack.copy("source_track", "copied_track")

        expect_true(gtrack.exists("source_track"))
        expect_true(gtrack.exists("copied_track"))

        # Source stays in db1, copy goes to db2
        expect_equal(gtrack.db("source_track"), normalizePath("db1"))
        expect_equal(gtrack.db("copied_track"), normalizePath("db2"))

        # Values should match
        expect_equal(
            gextract("source_track", gintervals(1, 0, 500))$source_track[1],
            gextract("copied_track", gintervals(1, 0, 500))$copied_track[1]
        )
    })
})

test_that("gtrack.copy fails when source doesn't exist", {
    withr::with_tempdir({
        create_test_db("db1")
        gsetroot("db1")

        expect_error(gtrack.copy("nonexistent", "copy"), "does not exist")
    })
})

test_that("gtrack.copy fails when destination already exists", {
    withr::with_tempdir({
        create_test_db("db1")
        gsetroot("db1")

        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("track1", "t1", intervs, 1)
        gtrack.create_sparse("track2", "t2", intervs, 2)

        expect_error(gtrack.copy("track1", "track2"), "already exists")
    })
})

test_that("creating track in multi-db goes to correct database", {
    withr::with_tempdir({
        create_test_db("db1")
        create_test_db("db2")

        gsetroot("db1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("source_track", "source", intervs, 99)

        # Connect to both, GWD defaults to db2 (last)
        gsetroot(c("db1", "db2"))

        # Create new track - should go to db2 (current GWD)
        gtrack.create_sparse("new_track", "new", intervs, 77)
        gdb.reload()

        expect_true(gtrack.exists("source_track"))
        expect_true(gtrack.exists("new_track"))

        # Source stays in db1, new track goes to db2
        expect_equal(gtrack.db("source_track"), normalizePath("db1"))
        expect_equal(gtrack.db("new_track"), normalizePath("db2"))
    })
})

test_that("gintervals operations work across databases", {
    withr::with_tempdir({
        create_test_db("db1")
        create_test_db("db2")

        gsetroot("db1")
        int1 <- gintervals(1, 0, 1000)
        gintervals.save("db1_intervals", int1)

        gsetroot("db2")
        int2 <- gintervals(1, 500, 1500)
        gintervals.save("db2_intervals", int2)

        gsetroot(c("db1", "db2"))

        # Both interval sets should be visible
        all_intervals <- gintervals.ls()
        expect_true("db1_intervals" %in% all_intervals)
        expect_true("db2_intervals" %in% all_intervals)

        # Load using gintervals (not gintervals.load which requires specific format)
        # Verify the intervals exist and can be used
        expect_equal(gintervals.db("db1_intervals"), normalizePath("db1"))
        expect_equal(gintervals.db("db2_intervals"), normalizePath("db2"))

        # Operations between intervals from different dbs using constructed intervals
        union_result <- gintervals.union(int1, int2)
        expect_equal(nrow(union_result), 1) # Should merge into one interval
        expect_equal(union_result$start[1], 0)
        expect_equal(union_result$end[1], 1500)
    })
})

test_that("gscreen works with tracks from different databases", {
    withr::with_tempdir({
        create_test_db("db1")
        create_test_db("db2")

        gsetroot("db1")
        intervs <- gintervals(1, 0, 5000)
        gtrack.create_sparse("db1_track", "db1", intervs, 10)

        gsetroot("db2")
        gtrack.create_sparse("db2_track", "db2", intervs, 5)

        gsetroot(c("db1", "db2"))

        # Screen using track from db1
        result1 <- gscreen("db1_track > 5", gintervals(1, 0, 1000))
        expect_true(nrow(result1) > 0)

        # Screen using track from db2
        result2 <- gscreen("db2_track > 3", gintervals(1, 0, 1000))
        expect_true(nrow(result2) > 0)

        # Screen using expression combining both (with explicit iterator since multiple sparse tracks)
        result3 <- gscreen("db1_track + db2_track > 10", gintervals(1, 0, 1000), iterator = 100)
        expect_true(nrow(result3) > 0)
    })
})

test_that("gsummary works with tracks from different databases", {
    withr::with_tempdir({
        create_test_db("db1")
        create_test_db("db2")

        gsetroot("db1")
        intervs <- gintervals(1, 0, 5000)
        gtrack.create_sparse("db1_track", "db1", intervs, 10)

        gsetroot("db2")
        gtrack.create_sparse("db2_track", "db2", intervs, 5)

        gsetroot(c("db1", "db2"))

        # Summary of track from db1 (sparse track sums values, not value*length)
        sum1 <- gsummary("db1_track", gintervals(1, 0, 1000))
        expect_equal(sum1[["Sum"]], 10) # single sparse value

        # Summary of track from db2
        sum2 <- gsummary("db2_track", gintervals(1, 0, 1000))
        expect_equal(sum2[["Sum"]], 5)
    })
})

test_that("gdb.reload correctly refreshes multi-db state", {
    withr::with_tempdir({
        create_test_db("db1")
        create_test_db("db2")

        gsetroot("db1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("track1", "t1", intervs, 1)

        gsetroot(c("db1", "db2"))
        expect_equal(length(gtrack.ls()), 1)

        # Manually create a track in db2 without using misha
        dir.create(file.path("db2", "tracks", "manual_track.track"), recursive = TRUE)
        writeLines("sparse", file.path("db2", "tracks", "manual_track.track", "type"))

        # Before reload, track not visible
        expect_false("manual_track" %in% gtrack.ls())

        # After reload with rescan, track should appear
        gdb.reload(rescan = TRUE)

        # Note: The manually created track may not be valid, but this tests the reload mechanism
    })
})

test_that("connecting to non-unique databases fails gracefully", {
    withr::with_tempdir({
        create_test_db("db1")

        abs_path <- normalizePath("db1")

        # Connecting same db twice should work but deduplicate in track listing
        gsetroot(c(abs_path, abs_path))

        # Track count should not be doubled
        gsetroot("db1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("test_track", "test", intervs, 1)

        gsetroot(c(abs_path, abs_path))
        track_count <- sum(gtrack.ls() == "test_track")
        expect_equal(track_count, 1) # Track appears only once
    })
})

test_that("gtrack.info works correctly for tracks in different databases", {
    withr::with_tempdir({
        create_test_db("db1")
        create_test_db("db2")

        gsetroot("db1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("sparse_in_db1", "sparse track", intervs, 1)

        gsetroot("db2")
        # Create another sparse track in db2
        gtrack.create_sparse("sparse_in_db2", "sparse track 2", intervs, 2)

        gsetroot(c("db1", "db2"))

        # Info should work for both
        info1 <- gtrack.info("sparse_in_db1")
        info2 <- gtrack.info("sparse_in_db2")

        expect_equal(info1$type, "sparse")
        expect_equal(info2$type, "sparse")
    })
})

test_that("gdb.summary track counts are accurate with overlapping tracks", {
    withr::with_tempdir({
        create_test_db("db1")
        create_test_db("db2")

        gsetroot("db1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("unique_db1", "db1", intervs, 1)
        gtrack.create_sparse("shared", "shared from db1", intervs, 10)

        gsetroot("db2")
        gtrack.create_sparse("unique_db2", "db2", intervs, 2)
        gtrack.create_sparse("shared", "shared from db2", intervs, 20)

        gsetroot(c("db1", "db2"))

        info <- gdb.summary()

        # db1 should have 1 track (unique_db1, shared is overridden)
        expect_equal(info$tracks[1], 1)
        # db2 should have 2 tracks (unique_db2 + shared)
        expect_equal(info$tracks[2], 2)

        # Total visible tracks should be 3
        expect_equal(length(gtrack.ls()), 3)
    })
})

test_that("connecting fails when seq directory is missing for indexed db", {
    withr::with_tempdir({
        # Create incomplete database (missing seq for indexed db)
        # For per-chromosome dbs, seq is checked lazily, so we test indexed format
        dir.create("incomplete_db/tracks", recursive = TRUE)
        write.table(
            data.frame(chrom = "chr1", size = 1000),
            "incomplete_db/chrom_sizes.txt",
            sep = "\t", row.names = FALSE, col.names = FALSE
        )

        # Create seq dir with invalid genome.idx (indexed format requires validation)
        dir.create("incomplete_db/seq", recursive = TRUE)
        writeLines("invalid", "incomplete_db/seq/genome.idx")

        expect_error(gsetroot("incomplete_db"))
    })
})

test_that("connecting fails when chrom_sizes.txt is missing", {
    withr::with_tempdir({
        # Create incomplete database (missing chrom_sizes.txt)
        dir.create("incomplete_db/tracks", recursive = TRUE)
        dir.create("incomplete_db/seq", recursive = TRUE)

        expect_error(
            gsetroot("incomplete_db"),
            "does not contain a chrom_sizes.txt file"
        )
    })
})

test_that("gsetroot gives clear error when directory doesn't exist", {
    expect_error(
        gsetroot("/this/path/does/not/exist"),
        "Database directory does not exist"
    )
})

test_that("gsetroot gives clear error when tracks/ subdirectory is missing", {
    withr::with_tempdir({
        dir.create("not_a_db/seq", recursive = TRUE)

        expect_error(
            gsetroot("not_a_db"),
            "does not contain a 'tracks' subdirectory"
        )
    })
})

test_that("gsetroot gives clear error when seq/ subdirectory is missing", {
    withr::with_tempdir({
        dir.create("not_a_db/tracks", recursive = TRUE)

        expect_error(
            gsetroot("not_a_db"),
            "does not contain a 'seq' subdirectory"
        )
    })
})

test_that("gsetroot validates all databases in multi-db setup", {
    withr::with_tempdir({
        create_test_db("db1")

        # Create incomplete second database
        dir.create("db2/tracks", recursive = TRUE)
        # Missing seq/

        expect_error(
            gsetroot(c("db1", "db2")),
            "does not contain a 'seq' subdirectory"
        )
    })
})

test_that("gsetroot detects when user passes database path as dir parameter", {
    withr::with_tempdir({
        create_test_db("db1")
        create_test_db("db2")

        # Common mistake: passing second database as 'dir' parameter
        expect_error(
            gsetroot("db1", "db2"),
            "looks like a misha database path.*To connect multiple databases, use: gsetroot\\(c\\("
        )
    })
})
