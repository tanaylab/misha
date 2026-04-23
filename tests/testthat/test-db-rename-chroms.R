test_that(".misha_rename_normalize_mapping accepts a data.frame", {
    m <- .misha_rename_normalize_mapping(data.frame(
        old = c("chr1", "chr2"),
        new = c("1", "2"),
        stringsAsFactors = FALSE
    ))
    expect_equal(m$old, c("chr1", "chr2"))
    expect_equal(m$new, c("1", "2"))
})

test_that(".misha_rename_normalize_mapping accepts a named character vector", {
    m <- .misha_rename_normalize_mapping(c(chr1 = "1", chr2 = "2"))
    expect_equal(m$old, c("chr1", "chr2"))
    expect_equal(m$new, c("1", "2"))
})

test_that(".misha_rename_normalize_mapping rejects unnamed vectors", {
    expect_error(
        .misha_rename_normalize_mapping(c("1", "2")),
        "named character vector"
    )
})

test_that(".misha_rename_normalize_mapping rejects missing columns", {
    expect_error(
        .misha_rename_normalize_mapping(data.frame(a = "x", b = "y")),
        "columns"
    )
})

test_that(".misha_rename_normalize_mapping rejects non-character inputs", {
    expect_error(.misha_rename_normalize_mapping(list(a = 1, b = 2)))
})

test_that(".misha_rename_validate_mapping rejects empty mapping", {
    m <- data.frame(old = character(), new = character(), stringsAsFactors = FALSE)
    expect_error(
        .misha_rename_validate_mapping(m, existing = c("chr1", "chr2")),
        "empty"
    )
})

test_that(".misha_rename_validate_mapping rejects unknown old names", {
    m <- data.frame(old = c("chrX"), new = c("X"), stringsAsFactors = FALSE)
    expect_error(
        .misha_rename_validate_mapping(m, existing = c("chr1", "chr2")),
        "not present"
    )
})

test_that(".misha_rename_validate_mapping rejects duplicate old names", {
    m <- data.frame(
        old = c("chr1", "chr1"), new = c("A", "B"),
        stringsAsFactors = FALSE
    )
    expect_error(
        .misha_rename_validate_mapping(m, existing = c("chr1", "chr2")),
        "duplicate.*old"
    )
})

test_that(".misha_rename_validate_mapping rejects duplicate new names", {
    m <- data.frame(
        old = c("chr1", "chr2"), new = c("X", "X"),
        stringsAsFactors = FALSE
    )
    expect_error(
        .misha_rename_validate_mapping(m, existing = c("chr1", "chr2")),
        "duplicate.*new"
    )
})

test_that(".misha_rename_validate_mapping rejects collision with un-mapped chrom", {
    m <- data.frame(
        old = c("chr1"), new = c("chr2"),
        stringsAsFactors = FALSE
    )
    expect_error(
        .misha_rename_validate_mapping(m, existing = c("chr1", "chr2", "chr3")),
        "collides"
    )
})

test_that(".misha_rename_validate_mapping accepts a valid swap", {
    m <- data.frame(
        old = c("chr1", "chr2"), new = c("chr2", "chr1"),
        stringsAsFactors = FALSE
    )
    expect_silent(.misha_rename_validate_mapping(m, existing = c("chr1", "chr2", "chr3")))
})

test_that(".misha_rename_validate_mapping accepts a valid partial rename", {
    m <- data.frame(
        old = c("chr1"), new = c("1"),
        stringsAsFactors = FALSE
    )
    expect_silent(.misha_rename_validate_mapping(m, existing = c("chr1", "chr2", "chr3")))
})
