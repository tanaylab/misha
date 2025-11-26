create_isolated_test_db()
gdir.create("temp", showWarnings = FALSE)

test_that("get readonly attributes of gdb", {
    gdb.set_readonly_attrs(c("created.by", "created.date", "created.user"))
    r <- gdb.get_readonly_attrs()
    expect_equal(r, c("created.by", "created.date", "created.user"))
})

test_that("set and revert readonly attributes of gdb", {
    gdb.set_readonly_attrs(c("created.by", "created.date", "created.user"))
    old_ro_attrs <- gdb.get_readonly_attrs()
    expect_error(gdb.set_readonly_attrs(c(old_ro_attrs, "")))
    gdb.set_readonly_attrs(old_ro_attrs)
    expect_equal(gdb.get_readonly_attrs(), old_ro_attrs)
})

test_that("set NULL and revert readonly attributes", {
    gdb.set_readonly_attrs(c("created.by", "created.date", "created.user"))
    old_ro_attrs <- gdb.get_readonly_attrs()
    gdb.set_readonly_attrs(NULL)
    expect_null(gdb.get_readonly_attrs())
    gdb.set_readonly_attrs(old_ro_attrs)
    expect_equal(gdb.get_readonly_attrs(), old_ro_attrs)
})

test_that("append and revert readonly attributes", {
    gdb.set_readonly_attrs(c("created.by", "created.date", "created.user"))
    old_ro_attrs <- gdb.get_readonly_attrs()
    gdb.set_readonly_attrs(c(old_ro_attrs, "testattr1", "testattr2"))
    expect_equal(gdb.get_readonly_attrs(), union(old_ro_attrs, c("testattr1", "testattr2")))
    gdb.set_readonly_attrs(old_ro_attrs)
    expect_equal(gdb.get_readonly_attrs(), old_ro_attrs)
})

test_that("get attributes of track", {
    gdb.set_readonly_attrs(c("created.by", "created.date", "created.user"))
    tmptrack <- paste0("temp.tmptrack_", sample(1:1e9, 1))
    gtrack.rm(tmptrack, force = TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.create_sparse(tmptrack, "Test track description", gintervals(c(1, 2)), 1:2)

    # Check that created.by attribute exists and contains expected pattern
    created_by <- gtrack.attr.get(tmptrack, "created.by")
    expect_true(nchar(created_by) > 0)
    expect_match(created_by, "gtrack.create_sparse")

    # Check that non-existent attribute returns empty string
    expect_equal(gtrack.attr.get(tmptrack, "blablabla"), "")
})

test_that("set and get attribute after removing from readonly", {
    gdb.set_readonly_attrs(c("created.by", "created.date", "created.user"))
    attrs <- gdb.get_readonly_attrs()
    gdb.set_readonly_attrs(c(attrs, "testattr1"))
    attrs <- attrs[attrs != "testattr1"]
    gdb.set_readonly_attrs(attrs)

    tmptrack <- paste0("temp.tmptrack_", sample(1:1e9, 1))
    gtrack.rm(tmptrack, force = TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.create_sparse(tmptrack, "Test track", gintervals(c(1, 2)), 1:2)

    gtrack.attr.set(tmptrack, "testattr1", "value")
    r <- gtrack.attr.get(tmptrack, "testattr1")
    expect_equal(r, "value")
})

test_that("set, reset and export attribute after removing from readonly", {
    gdb.set_readonly_attrs(c("created.by", "created.date", "created.user"))
    attrs <- gdb.get_readonly_attrs()
    gdb.set_readonly_attrs(c(attrs, "testattr1"))
    attrs <- attrs[attrs != "testattr1"]
    gdb.set_readonly_attrs(attrs)

    tmptrack <- paste0("temp.tmptrack_", sample(1:1e9, 1))
    gtrack.rm(tmptrack, force = TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.create_sparse(tmptrack, "Test track", gintervals(c(1, 2)), 1:2)

    gtrack.attr.set(tmptrack, "testattr1", "value")
    # Verify it was set
    expect_equal(gtrack.attr.get(tmptrack, "testattr1"), "value")

    gtrack.attr.set(tmptrack, "testattr1", "")
    # When set to empty string, attribute should be removed
    # gtrack.attr.get returns "" for non-existent attributes
    expect_equal(gtrack.attr.get(tmptrack, "testattr1"), "")

    r <- gtrack.attr.export(tmptrack)
    # Verify testattr1 is either missing from export or empty string
    if ("testattr1" %in% names(r)) {
        expect_equal(r[tmptrack, "testattr1"], "")
    }
    # If not in export, that's also valid (attribute was removed)
})

test_that("export track attributes without parameters", {
    gdb.set_readonly_attrs(c("created.by", "created.date", "created.user"))

    # Create test tracks with known attributes
    tmptrack1 <- paste0("temp.tmptrack1_", sample(1:1e9, 1))
    tmptrack2 <- paste0("temp.tmptrack2_", sample(1:1e9, 1))
    gtrack.rm(tmptrack1, force = TRUE)
    gtrack.rm(tmptrack2, force = TRUE)
    withr::defer(gtrack.rm(tmptrack1, force = TRUE))
    withr::defer(gtrack.rm(tmptrack2, force = TRUE))

    gtrack.create_sparse(tmptrack1, "Test track 1 description", gintervals(c(1, 2)), 1:2)
    gtrack.create_sparse(tmptrack2, "Test track 2 description", gintervals(c(1, 2)), 3:4)

    # Export all attributes
    r <- gtrack.attr.export()

    # Verify our test tracks are in the export
    expect_true(tmptrack1 %in% rownames(r))
    expect_true(tmptrack2 %in% rownames(r))

    # Verify expected attributes exist
    expect_true("created.by" %in% names(r))
    expect_true("created.date" %in% names(r))
    expect_true("description" %in% names(r))

    # Verify our track descriptions are correct
    expect_equal(r[tmptrack1, "description"], "Test track 1 description")
    expect_equal(r[tmptrack2, "description"], "Test track 2 description")

    # Verify created.by contains expected pattern
    expect_match(r[tmptrack1, "created.by"], "gtrack.create_sparse")
    expect_match(r[tmptrack2, "created.by"], "gtrack.create_sparse")
})

test_that("export attributes of non-existing track", {
    gdb.set_readonly_attrs(c("created.by", "created.date", "created.user"))
    expect_error(gtrack.attr.export("blablablatrack"))
})

test_that("export specific attributes of tracks", {
    gdb.set_readonly_attrs(c("created.by", "created.date", "created.user"))

    # Create two test tracks with known attributes
    tmptrack1 <- paste0("temp.tmptrack1_", sample(1:1e9, 1))
    tmptrack2 <- paste0("temp.tmptrack2_", sample(1:1e9, 1))
    gtrack.rm(tmptrack1, force = TRUE)
    gtrack.rm(tmptrack2, force = TRUE)
    withr::defer(gtrack.rm(tmptrack1, force = TRUE))
    withr::defer(gtrack.rm(tmptrack2, force = TRUE))

    gtrack.create_sparse(tmptrack1, "Test track 1", gintervals(c(1, 2)), 1:2)
    gtrack.create_sparse(tmptrack2, "Test track 2", gintervals(c(1, 2)), 3:4)

    # Test exporting single attribute
    r1 <- gtrack.attr.export(tracks = c(tmptrack1, tmptrack2), attrs = c("created.by"))
    expect_true(tmptrack1 %in% rownames(r1))
    expect_true(tmptrack2 %in% rownames(r1))
    expect_true("created.by" %in% names(r1))
    expect_true(nchar(r1[tmptrack1, "created.by"]) > 0)
    expect_true(nchar(r1[tmptrack2, "created.by"]) > 0)

    # Test exporting multiple attributes
    r2 <- gtrack.attr.export(tracks = c(tmptrack1, tmptrack2), attrs = c("created.by", "created.date"))
    expect_true("created.by" %in% names(r2))
    expect_true("created.date" %in% names(r2))
    expect_true(nchar(r2[tmptrack1, "created.date"]) > 0)

    # Test exporting all attributes of specific tracks
    r3 <- gtrack.attr.export(c(tmptrack1, tmptrack2))
    expect_true(tmptrack1 %in% rownames(r3))
    expect_true(tmptrack2 %in% rownames(r3))
    expect_true("created.by" %in% names(r3))
    expect_true("description" %in% names(r3))
    expect_equal(r3[tmptrack1, "description"], "Test track 1")
    expect_equal(r3[tmptrack2, "description"], "Test track 2")

    # Test exporting specific attributes of specific tracks
    r4 <- gtrack.attr.export(c(tmptrack1, tmptrack2), attrs = c("created.by", "created.date"))
    expect_equal(ncol(r4), 2)
    expect_true("created.by" %in% names(r4))
    expect_true("created.date" %in% names(r4))
})

test_that("export after importing modified attributes", {
    gdb.set_readonly_attrs(c("created.by", "created.date", "created.user"))
    attrs <- gdb.get_readonly_attrs()
    gdb.set_readonly_attrs(c(attrs, "testattr1"))
    attrs <- attrs[attrs != "testattr1"]
    gdb.set_readonly_attrs(attrs)

    # Create test tracks
    tmptrack1 <- paste0("temp.tmptrack1_", sample(1:1e9, 1))
    tmptrack2 <- paste0("temp.tmptrack2_", sample(1:1e9, 1))
    gtrack.rm(tmptrack1, force = TRUE)
    gtrack.rm(tmptrack2, force = TRUE)
    withr::defer(gtrack.rm(tmptrack1, force = TRUE))
    withr::defer(gtrack.rm(tmptrack2, force = TRUE))

    gtrack.create_sparse(tmptrack1, "Test track 1", gintervals(c(1, 2)), 1:2)
    gtrack.create_sparse(tmptrack2, "Test track 2", gintervals(c(1, 2)), 3:4)

    # Export attributes and modify testattr1
    r1 <- gtrack.attr.export(tracks = c(tmptrack1, tmptrack2), attrs = c("testattr1", "created.by", "created.date", "created.user", "testattr2"))
    r1$testattr1 <- c("value1", "value2")
    gtrack.attr.import(r1)

    # Export again and verify the imported values
    r2 <- gtrack.attr.export(tracks = c(tmptrack1, tmptrack2), attrs = c("testattr1", "created.by", "created.date", "created.user", "testattr2"))
    expect_equal(r2[tmptrack1, "testattr1"], "value1")
    expect_equal(r2[tmptrack2, "testattr1"], "value2")
})

test_that("import attributes without replacing existing ones", {
    gdb.set_readonly_attrs(c("created.by", "created.date", "created.user"))
    attrs <- gdb.get_readonly_attrs()
    gdb.set_readonly_attrs(c(attrs, "testattr1"))
    attrs <- attrs[attrs != "testattr1"]
    gdb.set_readonly_attrs(attrs)

    # Create test tracks
    tmptrack1 <- paste0("temp.tmptrack1_", sample(1:1e9, 1))
    tmptrack2 <- paste0("temp.tmptrack2_", sample(1:1e9, 1))
    gtrack.rm(tmptrack1, force = TRUE)
    gtrack.rm(tmptrack2, force = TRUE)
    withr::defer(gtrack.rm(tmptrack1, force = TRUE))
    withr::defer(gtrack.rm(tmptrack2, force = TRUE))

    gtrack.create_sparse(tmptrack1, "Test track 1", gintervals(c(1, 2)), 1:2)
    gtrack.create_sparse(tmptrack2, "Test track 2", gintervals(c(1, 2)), 3:4)

    # Set initial testattr1 values
    r1 <- gtrack.attr.export(tracks = c(tmptrack1, tmptrack2), attrs = c("testattr1", "created.by", "created.date", "created.user", "testattr2"))
    r1$testattr1 <- c("initial1", "initial2")
    gtrack.attr.import(r1)

    # Try to import without testattr1 (should not replace existing)
    r1$testattr1 <- NULL
    gtrack.attr.import(r1)

    # Verify testattr1 values are still there
    r2 <- gtrack.attr.export(tracks = c(tmptrack1, tmptrack2), attrs = c("testattr1", "created.by", "created.date", "created.user", "testattr2"))
    expect_equal(r2[tmptrack1, "testattr1"], "initial1")
    expect_equal(r2[tmptrack2, "testattr1"], "initial2")
})

test_that("import attributes with replacing existing ones", {
    gdb.set_readonly_attrs(c("created.by", "created.date", "created.user"))
    attrs <- gdb.get_readonly_attrs()
    gdb.set_readonly_attrs(c(attrs, "testattr1"))
    attrs <- attrs[attrs != "testattr1"]
    gdb.set_readonly_attrs(attrs)

    # Create test tracks
    tmptrack1 <- paste0("temp.tmptrack1_", sample(1:1e9, 1))
    tmptrack2 <- paste0("temp.tmptrack2_", sample(1:1e9, 1))
    gtrack.rm(tmptrack1, force = TRUE)
    gtrack.rm(tmptrack2, force = TRUE)
    withr::defer(gtrack.rm(tmptrack1, force = TRUE))
    withr::defer(gtrack.rm(tmptrack2, force = TRUE))

    gtrack.create_sparse(tmptrack1, "Test track 1", gintervals(c(1, 2)), 1:2)
    gtrack.create_sparse(tmptrack2, "Test track 2", gintervals(c(1, 2)), 3:4)

    # Set initial testattr1 values
    r1 <- gtrack.attr.export(tracks = c(tmptrack1, tmptrack2), attrs = c("testattr1", "created.by", "created.date", "created.user", "testattr2"))
    r1$testattr1 <- c("initial1", "initial2")
    gtrack.attr.import(r1)

    # Replace with new values
    r1$testattr1 <- c("replaced1", "replaced2")
    gtrack.attr.import(r1)

    # Verify values were replaced
    r2 <- gtrack.attr.export(tracks = c(tmptrack1, tmptrack2), attrs = c("testattr1", "created.by", "created.date", "created.user", "testattr2"))
    expect_equal(r2[tmptrack1, "testattr1"], "replaced1")
    expect_equal(r2[tmptrack2, "testattr1"], "replaced2")
})

test_that("import selected track attributes", {
    gdb.set_readonly_attrs(c("created.by", "created.date", "created.user"))
    attrs <- gdb.get_readonly_attrs()
    gdb.set_readonly_attrs(c(attrs, "testattr1"))
    attrs <- attrs[attrs != "testattr1"]
    gdb.set_readonly_attrs(attrs)

    # Create test tracks
    tmptrack1 <- paste0("temp.tmptrack1_", sample(1:1e9, 1))
    tmptrack2 <- paste0("temp.tmptrack2_", sample(1:1e9, 1))
    tmptrack3 <- paste0("temp.tmptrack3_", sample(1:1e9, 1))
    gtrack.rm(tmptrack1, force = TRUE)
    gtrack.rm(tmptrack2, force = TRUE)
    gtrack.rm(tmptrack3, force = TRUE)
    withr::defer(gtrack.rm(tmptrack1, force = TRUE))
    withr::defer(gtrack.rm(tmptrack2, force = TRUE))
    withr::defer(gtrack.rm(tmptrack3, force = TRUE))

    gtrack.create_sparse(tmptrack1, "Test track 1", gintervals(c(1, 2)), 1:2)
    gtrack.create_sparse(tmptrack2, "Test track 2", gintervals(c(1, 2)), 3:4)
    gtrack.create_sparse(tmptrack3, "Test track 3", gintervals(c(1, 2)), 5:6)

    # Export all tracks but only import for selected tracks
    r1 <- gtrack.attr.export(tracks = c(tmptrack1, tmptrack2, tmptrack3), attrs = c("testattr1", "created.by", "created.date", "created.user", "testattr2"))
    r1$testattr1 <- c("value1", "value2", "value3")
    # Only import for first two tracks
    r1_selected <- r1[c(tmptrack1, tmptrack2), , drop = FALSE]
    gtrack.attr.import(r1_selected)

    # Verify only selected tracks got the attribute
    r2 <- gtrack.attr.export(tracks = c(tmptrack1, tmptrack2, tmptrack3), attrs = c("testattr1", "created.by", "created.date", "created.user", "testattr2"))
    expect_equal(r2[tmptrack1, "testattr1"], "value1")
    expect_equal(r2[tmptrack2, "testattr1"], "value2")
    # Third track should not have testattr1 set (or be empty)
    if (tmptrack3 %in% rownames(r2) && "testattr1" %in% names(r2)) {
        expect_true(is.na(r2[tmptrack3, "testattr1"]) || r2[tmptrack3, "testattr1"] == "")
    }
})

test_that("import with testattr1 in readonly attributes", {
    gdb.set_readonly_attrs(c("created.by", "created.date", "created.user"))
    attrs <- gdb.get_readonly_attrs()
    gdb.set_readonly_attrs(c(attrs, "testattr1"))
    attrs <- attrs[attrs != "testattr1"]
    gdb.set_readonly_attrs(c(attrs, "testattr1"))

    # Create test tracks
    tmptrack1 <- paste0("temp.tmptrack1_", sample(1:1e9, 1))
    tmptrack2 <- paste0("temp.tmptrack2_", sample(1:1e9, 1))
    gtrack.rm(tmptrack1, force = TRUE)
    gtrack.rm(tmptrack2, force = TRUE)
    withr::defer(gtrack.rm(tmptrack1, force = TRUE))
    withr::defer(gtrack.rm(tmptrack2, force = TRUE))

    gtrack.create_sparse(tmptrack1, "Test track 1", gintervals(c(1, 2)), 1:2)
    gtrack.create_sparse(tmptrack2, "Test track 2", gintervals(c(1, 2)), 3:4)

    # Try to import readonly attribute - should fail
    r1 <- gtrack.attr.export(tracks = c(tmptrack1, tmptrack2))
    r1$testattr1 <- c("value1", "value2")
    expect_error(gtrack.attr.import(r1))
})
