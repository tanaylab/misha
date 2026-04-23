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

test_that(".misha_rename_needs_two_phase detects swaps", {
    m <- data.frame(old = c("a", "b"), new = c("b", "a"), stringsAsFactors = FALSE)
    expect_true(.misha_rename_needs_two_phase(m))
})

test_that(".misha_rename_needs_two_phase detects cycles of length > 2", {
    m <- data.frame(old = c("a", "b", "c"), new = c("b", "c", "a"),
                    stringsAsFactors = FALSE)
    expect_true(.misha_rename_needs_two_phase(m))
})

test_that(".misha_rename_needs_two_phase returns FALSE for non-overlapping rename", {
    m <- data.frame(old = c("a", "b"), new = c("x", "y"), stringsAsFactors = FALSE)
    expect_false(.misha_rename_needs_two_phase(m))
})

test_that(".misha_rename_needs_two_phase returns FALSE for no-op entries", {
    m <- data.frame(old = c("a"), new = c("a"), stringsAsFactors = FALSE)
    expect_false(.misha_rename_needs_two_phase(m))
})

test_that(".misha_rename_normalize_mapping rejects NA names in vector", {
    v <- c("1", "2")
    names(v) <- c("chr1", NA)
    expect_error(
        .misha_rename_normalize_mapping(v),
        "named character vector"
    )
})

test_that(".misha_rename_validate_mapping rejects NA old", {
    m <- data.frame(old = c("chr1", NA), new = c("A", "B"),
                    stringsAsFactors = FALSE)
    expect_error(
        .misha_rename_validate_mapping(m, existing = c("chr1", "chr2")),
        "NA or empty"
    )
})

test_that(".misha_rename_validate_mapping rejects NA new", {
    m <- data.frame(old = c("chr1", "chr2"), new = c("A", NA),
                    stringsAsFactors = FALSE)
    expect_error(
        .misha_rename_validate_mapping(m, existing = c("chr1", "chr2")),
        "NA or empty"
    )
})

test_that(".misha_rename_validate_mapping rejects empty-string new", {
    m <- data.frame(old = c("chr1"), new = c(""),
                    stringsAsFactors = FALSE)
    expect_error(
        .misha_rename_validate_mapping(m, existing = c("chr1", "chr2")),
        "NA or empty"
    )
})

test_that(".misha_rename_remap_factor replaces matching levels and values", {
    f <- factor(c("chr1", "chr2", "chr1"), levels = c("chr1", "chr2", "chr3"))
    result <- .misha_rename_remap_factor(f, old = c("chr1", "chr2"), new = c("A", "B"))
    expect_equal(as.character(result), c("A", "B", "A"))
    expect_equal(levels(result), c("A", "B", "chr3"))
})

test_that(".misha_rename_remap_factor leaves unmapped levels alone", {
    f <- factor(c("chr1", "chr2"), levels = c("chr1", "chr2"))
    result <- .misha_rename_remap_factor(f, old = "chr1", new = "A")
    expect_equal(as.character(result), c("A", "chr2"))
    expect_equal(levels(result), c("A", "chr2"))
})

test_that(".misha_rename_remap_factor is a no-op on NULL or empty", {
    expect_null(.misha_rename_remap_factor(NULL, old = "a", new = "b"))
    empty_f <- factor(character(0), levels = c("chr1", "chr2"))
    result <- .misha_rename_remap_factor(empty_f, old = "chr1", new = "A")
    expect_equal(levels(result), c("A", "chr2"))
})

test_that(".misha_rename_remap_df remaps chrom column in 1D frames", {
    df <- data.frame(
        chrom = factor(c("chr1", "chr2"), levels = c("chr1", "chr2", "chr3")),
        start = c(0L, 0L), end = c(100L, 200L),
        stringsAsFactors = FALSE
    )
    out <- .misha_rename_remap_df(df, old = "chr1", new = "A")
    expect_equal(as.character(out$chrom), c("A", "chr2"))
    expect_equal(levels(out$chrom), c("A", "chr2", "chr3"))
})

test_that(".misha_rename_remap_df remaps chrom1/chrom2 in 2D frames", {
    df <- data.frame(
        chrom1 = factor(c("chr1"), levels = c("chr1", "chr2")),
        start1 = 0L, end1 = 100L,
        chrom2 = factor(c("chr2"), levels = c("chr1", "chr2")),
        start2 = 0L, end2 = 100L,
        stringsAsFactors = FALSE
    )
    out <- .misha_rename_remap_df(df,
        old = c("chr1", "chr2"),
        new = c("A",    "B"))
    expect_equal(as.character(out$chrom1), "A")
    expect_equal(as.character(out$chrom2), "B")
})

test_that(".misha_rename_remap_df tolerates 0-row frames (zerolines)", {
    df <- data.frame(
        chrom = factor(character(0), levels = c("chr1", "chr2")),
        start = integer(0), end = integer(0),
        stringsAsFactors = FALSE
    )
    out <- .misha_rename_remap_df(df, old = "chr1", new = "A")
    expect_equal(levels(out$chrom), c("A", "chr2"))
})

test_that(".misha_rename_remap_df is a no-op on NULL (zeroline path)", {
    expect_null(.misha_rename_remap_df(NULL, old = "chr1", new = "A"))
})

test_that(".misha_rename_atomic_rewrite replaces file atomically", {
    target <- tempfile()
    writeLines("old content", target)
    .misha_rename_atomic_rewrite(target, function(tmp_path) {
        writeLines("new content", tmp_path)
    })
    expect_equal(readLines(target), "new content")
    unlink(target)
})

test_that(".misha_rename_atomic_rewrite leaves original intact on writer error", {
    target <- tempfile()
    writeLines("original", target)
    expect_error(
        .misha_rename_atomic_rewrite(target, function(tmp_path) {
            writeLines("partial", tmp_path)
            stop("boom")
        }),
        "boom"
    )
    expect_equal(readLines(target), "original")
    tmps <- list.files(
        dirname(target),
        pattern = paste0(basename(target), "\\.tmp\\..*"),
        full.names = TRUE
    )
    expect_length(tmps, 0L)
    unlink(target)
})

test_that(".misha_rename_rewrite_meta updates stats + zeroline factor levels", {
    path <- tempfile()
    dir.create(path)
    meta_file <- file.path(path, ".meta")

    stats <- data.frame(
        chrom = factor(c("chr1", "chr2"), levels = c("chr1", "chr2")),
        size = c(100, 200),
        stringsAsFactors = FALSE
    )
    zeroline <- data.frame(
        chrom = factor(character(0), levels = c("chr1", "chr2")),
        start = integer(0), end = integer(0),
        stringsAsFactors = FALSE
    )
    f <- file(meta_file, "wb"); serialize(list(stats = stats, zeroline = zeroline), f); close(f)

    .misha_rename_rewrite_meta(meta_file,
        old = c("chr1", "chr2"), new = c("A", "B"))

    f <- file(meta_file, "rb"); m <- unserialize(f); close(f)
    expect_equal(as.character(m$stats$chrom), c("A", "B"))
    expect_equal(levels(m$zeroline$chrom), c("A", "B"))
    unlink(path, recursive = TRUE)
})

test_that(".misha_rename_rewrite_meta handles 2D stats", {
    path <- tempfile()
    dir.create(path)
    meta_file <- file.path(path, ".meta")

    stats <- data.frame(
        chrom1 = factor(c("chr1"), levels = c("chr1", "chr2")),
        chrom2 = factor(c("chr2"), levels = c("chr1", "chr2")),
        size = 5,
        stringsAsFactors = FALSE
    )
    zeroline <- data.frame(
        chrom1 = factor(character(0), levels = c("chr1", "chr2")),
        chrom2 = factor(character(0), levels = c("chr1", "chr2")),
        start1 = integer(0), end1 = integer(0),
        start2 = integer(0), end2 = integer(0),
        stringsAsFactors = FALSE
    )
    f <- file(meta_file, "wb"); serialize(list(stats = stats, zeroline = zeroline), f); close(f)

    .misha_rename_rewrite_meta(meta_file,
        old = c("chr1", "chr2"), new = c("A", "B"))

    f <- file(meta_file, "rb"); m <- unserialize(f); close(f)
    expect_equal(as.character(m$stats$chrom1), "A")
    expect_equal(as.character(m$stats$chrom2), "B")
    unlink(path, recursive = TRUE)
})

test_that(".misha_rename_rewrite_single_interv updates factor levels and values", {
    path <- tempfile(fileext = ".interv")
    df <- data.frame(
        chrom = factor(c("chr1", "chr2"), levels = c("chr1", "chr2", "chr3")),
        start = c(0L, 0L), end = c(100L, 200L),
        stringsAsFactors = FALSE
    )
    f <- file(path, "wb"); serialize(df, f); close(f)

    .misha_rename_rewrite_single_interv(path,
        old = c("chr1", "chr2"), new = c("A", "B"))

    f <- file(path, "rb"); out <- unserialize(f); close(f)
    expect_equal(as.character(out$chrom), c("A", "B"))
    expect_equal(levels(out$chrom), c("A", "B", "chr3"))
    unlink(path)
})
