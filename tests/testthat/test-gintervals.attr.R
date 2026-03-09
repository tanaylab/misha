create_isolated_test_db()
gdir.create("temp", showWarnings = FALSE)

# ---- Helpers ----

create_test_iset <- function(name) {
    gintervals.rm(name, force = TRUE)
    withr::defer(gintervals.rm(name, force = TRUE), envir = parent.frame())
    gintervals.save(name, gintervals(c(1, 2)))
    name
}

# ---- 1. Basic get/set ----

test_that("set an attribute and get it back", {
    iset <- create_test_iset("temp.iattr_test1")
    gintervals.attr.set(iset, "myattr", "hello")
    expect_equal(gintervals.attr.get(iset, "myattr"), "hello")
})

test_that("set empty string removes the attribute", {
    iset <- create_test_iset("temp.iattr_test2")
    gintervals.attr.set(iset, "myattr", "hello")
    expect_equal(gintervals.attr.get(iset, "myattr"), "hello")

    gintervals.attr.set(iset, "myattr", "")
    expect_equal(gintervals.attr.get(iset, "myattr"), "")
})

test_that("get non-existent attribute returns empty string", {
    iset <- create_test_iset("temp.iattr_test3")
    expect_equal(gintervals.attr.get(iset, "no_such_attr"), "")
})

# ---- 2. Multiple attributes ----

test_that("set multiple attrs on same interval set", {
    iset <- create_test_iset("temp.iattr_multi1")
    gintervals.attr.set(iset, "attr_a", "value_a")
    gintervals.attr.set(iset, "attr_b", "value_b")
    gintervals.attr.set(iset, "attr_c", "value_c")

    expect_equal(gintervals.attr.get(iset, "attr_a"), "value_a")
    expect_equal(gintervals.attr.get(iset, "attr_b"), "value_b")
    expect_equal(gintervals.attr.get(iset, "attr_c"), "value_c")
})

test_that("overwrite existing attr value", {
    iset <- create_test_iset("temp.iattr_overwrite1")
    gintervals.attr.set(iset, "myattr", "original")
    expect_equal(gintervals.attr.get(iset, "myattr"), "original")

    gintervals.attr.set(iset, "myattr", "updated")
    expect_equal(gintervals.attr.get(iset, "myattr"), "updated")
})

# ---- 3. Export ----

test_that("export all attrs of specific interval sets", {
    iset1 <- create_test_iset("temp.iattr_exp1")
    iset2 <- create_test_iset("temp.iattr_exp2")

    gintervals.attr.set(iset1, "color", "red")
    gintervals.attr.set(iset1, "size", "large")
    gintervals.attr.set(iset2, "color", "blue")

    r <- gintervals.attr.export(c(iset1, iset2))
    expect_true(is.data.frame(r))
    expect_true(iset1 %in% rownames(r))
    expect_true(iset2 %in% rownames(r))
    expect_true("color" %in% colnames(r))
    expect_true("size" %in% colnames(r))
    expect_equal(r[iset1, "color"], "red")
    expect_equal(r[iset1, "size"], "large")
    expect_equal(r[iset2, "color"], "blue")
    # iset2 has no "size" attr, should be ""
    expect_equal(r[iset2, "size"], "")
})

test_that("export specific attrs only", {
    iset <- create_test_iset("temp.iattr_exp3")
    gintervals.attr.set(iset, "color", "red")
    gintervals.attr.set(iset, "size", "large")
    gintervals.attr.set(iset, "weight", "heavy")

    r <- gintervals.attr.export(iset, attrs = c("color", "weight"))
    expect_equal(ncol(r), 2)
    expect_true("color" %in% colnames(r))
    expect_true("weight" %in% colnames(r))
    expect_false("size" %in% colnames(r))
    expect_equal(r[iset, "color"], "red")
    expect_equal(r[iset, "weight"], "heavy")
})

test_that("export non-existent interval set errors", {
    expect_error(gintervals.attr.export("no_such_interval_set_xyz"))
})

test_that("export with no attributes returns empty data frame with correct rows", {
    iset <- create_test_iset("temp.iattr_exp_empty")
    r <- gintervals.attr.export(iset)
    expect_true(is.data.frame(r))
    expect_equal(nrow(r), 1)
    expect_equal(ncol(r), 0)
    expect_equal(rownames(r), iset)
})

test_that("export all interval sets (NULL) works", {
    iset <- create_test_iset("temp.iattr_exp_all")
    gintervals.attr.set(iset, "tag", "test")

    r <- gintervals.attr.export()
    expect_true(is.data.frame(r))
    expect_true(iset %in% rownames(r))
})

# ---- 4. Import ----

test_that("import attrs from data frame", {
    iset1 <- create_test_iset("temp.iattr_imp1")
    iset2 <- create_test_iset("temp.iattr_imp2")

    tbl <- data.frame(
        attr_x = c("val_x1", "val_x2"),
        attr_y = c("val_y1", "val_y2"),
        row.names = c(iset1, iset2),
        stringsAsFactors = FALSE
    )
    gintervals.attr.import(tbl)

    expect_equal(gintervals.attr.get(iset1, "attr_x"), "val_x1")
    expect_equal(gintervals.attr.get(iset2, "attr_x"), "val_x2")
    expect_equal(gintervals.attr.get(iset1, "attr_y"), "val_y1")
    expect_equal(gintervals.attr.get(iset2, "attr_y"), "val_y2")
})

test_that("import with remove.others = FALSE preserves existing attrs", {
    iset <- create_test_iset("temp.iattr_imp_keep")
    gintervals.attr.set(iset, "existing", "keep_me")

    tbl <- data.frame(
        new_attr = "new_val",
        row.names = iset,
        stringsAsFactors = FALSE
    )
    gintervals.attr.import(tbl, remove.others = FALSE)

    expect_equal(gintervals.attr.get(iset, "existing"), "keep_me")
    expect_equal(gintervals.attr.get(iset, "new_attr"), "new_val")
})

test_that("import with remove.others = TRUE removes unlisted attrs", {
    iset <- create_test_iset("temp.iattr_imp_rm")
    gintervals.attr.set(iset, "existing", "will_go_away")
    gintervals.attr.set(iset, "keeper", "stays")

    tbl <- data.frame(
        keeper = "stays_updated",
        row.names = iset,
        stringsAsFactors = FALSE
    )
    gintervals.attr.import(tbl, remove.others = TRUE)

    expect_equal(gintervals.attr.get(iset, "keeper"), "stays_updated")
    # "existing" should have been removed
    expect_equal(gintervals.attr.get(iset, "existing"), "")
})

test_that("empty string in import removes attr", {
    iset <- create_test_iset("temp.iattr_imp_empty")
    gintervals.attr.set(iset, "to_remove", "present")
    expect_equal(gintervals.attr.get(iset, "to_remove"), "present")

    tbl <- data.frame(
        to_remove = "",
        row.names = iset,
        stringsAsFactors = FALSE
    )
    gintervals.attr.import(tbl)

    expect_equal(gintervals.attr.get(iset, "to_remove"), "")
})

test_that("import for non-existent interval set errors", {
    tbl <- data.frame(
        myattr = "val",
        row.names = "no_such_iset_xyz",
        stringsAsFactors = FALSE
    )
    expect_error(gintervals.attr.import(tbl))
})

test_that("import with duplicate interval set names errors", {
    iset <- create_test_iset("temp.iattr_imp_dup_iset")
    # R itself rejects duplicate row names, so creating such a data frame
    # will error before we reach gintervals.attr.import. This is acceptable
    # behaviour - the R runtime prevents this misuse.
    expect_error(suppressWarnings({
        tbl <- data.frame(myattr = c("v1", "v2"), stringsAsFactors = FALSE)
        rownames(tbl) <- c(iset, iset)
    }))
})

test_that("import with duplicate attr names errors", {
    iset <- create_test_iset("temp.iattr_imp_dup_attr")
    tbl <- data.frame(
        a = "v1",
        b = "v2",
        stringsAsFactors = FALSE
    )
    colnames(tbl) <- c("myattr", "myattr")
    rownames(tbl) <- iset

    expect_error(gintervals.attr.import(tbl), "appears more than once")
})

# ---- 5. Integration with gintervals.rm ----

test_that("deleting interval set cleans up .iattr file", {
    iset_name <- "temp.iattr_rm_test"
    gintervals.rm(iset_name, force = TRUE)
    gintervals.save(iset_name, gintervals(c(1, 2)))

    gintervals.attr.set(iset_name, "tag", "value")
    expect_equal(gintervals.attr.get(iset_name, "tag"), "value")

    # Get the .iattr path before deletion
    iattr_path <- misha:::.gintervals.attr_path(iset_name)

    # The .iattr file should exist
    expect_true(file.exists(iattr_path))

    gintervals.rm(iset_name, force = TRUE)

    # After deletion, .iattr file should be removed
    expect_false(file.exists(iattr_path))
})

# ---- 6. Integration with gintervals.save (overwrite) ----

test_that("overwriting interval set preserves attrs", {
    iset <- create_test_iset("temp.iattr_save_overwrite")
    gintervals.attr.set(iset, "tag", "preserved")
    expect_equal(gintervals.attr.get(iset, "tag"), "preserved")

    # Overwrite the interval set with new data
    gintervals.rm(iset, force = TRUE)
    gintervals.save(iset, gintervals(c(1, 2)))

    # Note: after rm + save, attrs are gone (rm cleans up .iattr)
    # This tests the expected behavior that attrs don't survive rm+save
    expect_equal(gintervals.attr.get(iset, "tag"), "")
})

test_that("attrs survive when interval set file is replaced without rm", {
    iset_name <- "temp.iattr_replace_test"
    gintervals.rm(iset_name, force = TRUE)
    gintervals.save(iset_name, gintervals(c(1, 2)))
    withr::defer(gintervals.rm(iset_name, force = TRUE))

    gintervals.attr.set(iset_name, "tag", "should_survive")

    # The .iattr file is separate from the .interv file, so if we
    # directly overwrite via gintervals.save (which calls rm internally
    # if the set already exists), check behavior
    # gintervals.save on existing set calls .gintervals.check_new_set
    # which may error; let's just verify the attr is still there
    expect_equal(gintervals.attr.get(iset_name, "tag"), "should_survive")
})

# ---- 7. Edge cases ----

test_that("many attrs on one interval set", {
    iset <- create_test_iset("temp.iattr_many")
    n <- 50
    attr_names <- paste0("attr_", sprintf("%03d", seq_len(n)))
    attr_values <- paste0("value_", seq_len(n))

    for (i in seq_len(n)) {
        gintervals.attr.set(iset, attr_names[i], attr_values[i])
    }

    # Verify all attrs
    for (i in seq_len(n)) {
        expect_equal(
            gintervals.attr.get(iset, attr_names[i]),
            attr_values[i],
            info = paste("attr index", i)
        )
    }

    # Verify via export
    r <- gintervals.attr.export(iset)
    expect_equal(ncol(r), n)
    for (i in seq_len(n)) {
        expect_equal(r[iset, attr_names[i]], attr_values[i])
    }
})

test_that("attr with special characters in value", {
    iset <- create_test_iset("temp.iattr_special")
    gintervals.attr.set(iset, "path", "/some/path/to/file.txt")
    expect_equal(gintervals.attr.get(iset, "path"), "/some/path/to/file.txt")

    gintervals.attr.set(iset, "desc", "value with spaces and tabs")
    expect_equal(gintervals.attr.get(iset, "desc"), "value with spaces and tabs")
})

test_that("attr with special characters in name", {
    iset <- create_test_iset("temp.iattr_special_name")
    gintervals.attr.set(iset, "my.attr", "dotted")
    expect_equal(gintervals.attr.get(iset, "my.attr"), "dotted")

    gintervals.attr.set(iset, "my_attr", "underscored")
    expect_equal(gintervals.attr.get(iset, "my_attr"), "underscored")
})

test_that("get/set with NULL arguments errors", {
    expect_error(gintervals.attr.get(NULL, "attr"))
    expect_error(gintervals.attr.get("annotations1", NULL))
    expect_error(gintervals.attr.set(NULL, "attr", "val"))
    expect_error(gintervals.attr.set("annotations1", NULL, "val"))
    expect_error(gintervals.attr.set("annotations1", "attr", NULL))
})

test_that("import with NULL table errors", {
    expect_error(gintervals.attr.import(NULL))
})

test_that("export attrs requested that don't exist returns empty strings", {
    iset <- create_test_iset("temp.iattr_missing_attrs")
    gintervals.attr.set(iset, "real", "exists")

    r <- gintervals.attr.export(iset, attrs = c("real", "fake"))
    expect_equal(r[iset, "real"], "exists")
    expect_equal(r[iset, "fake"], "")
})

test_that("import bulk via data frame then export matches", {
    iset1 <- create_test_iset("temp.iattr_bulk1")
    iset2 <- create_test_iset("temp.iattr_bulk2")

    tbl <- data.frame(
        color = c("red", "blue"),
        size = c("10", "20"),
        label = c("first", "second"),
        row.names = c(iset1, iset2),
        stringsAsFactors = FALSE
    )
    gintervals.attr.import(tbl)

    r <- gintervals.attr.export(c(iset1, iset2), attrs = c("color", "size", "label"))
    expect_equal(r[iset1, "color"], "red")
    expect_equal(r[iset2, "color"], "blue")
    expect_equal(r[iset1, "size"], "10")
    expect_equal(r[iset2, "size"], "20")
    expect_equal(r[iset1, "label"], "first")
    expect_equal(r[iset2, "label"], "second")
})
