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
    m <- data.frame(
        old = c("a", "b", "c"), new = c("b", "c", "a"),
        stringsAsFactors = FALSE
    )
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
    m <- data.frame(
        old = c("chr1", NA), new = c("A", "B"),
        stringsAsFactors = FALSE
    )
    expect_error(
        .misha_rename_validate_mapping(m, existing = c("chr1", "chr2")),
        "NA or empty"
    )
})

test_that(".misha_rename_validate_mapping rejects NA new", {
    m <- data.frame(
        old = c("chr1", "chr2"), new = c("A", NA),
        stringsAsFactors = FALSE
    )
    expect_error(
        .misha_rename_validate_mapping(m, existing = c("chr1", "chr2")),
        "NA or empty"
    )
})

test_that(".misha_rename_validate_mapping rejects empty-string new", {
    m <- data.frame(
        old = c("chr1"), new = c(""),
        stringsAsFactors = FALSE
    )
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
        new = c("A", "B")
    )
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
    f <- file(meta_file, "wb")
    serialize(list(stats = stats, zeroline = zeroline), f)
    close(f)

    .misha_rename_rewrite_meta(meta_file,
        old = c("chr1", "chr2"), new = c("A", "B")
    )

    f <- file(meta_file, "rb")
    m <- unserialize(f)
    close(f)
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
    f <- file(meta_file, "wb")
    serialize(list(stats = stats, zeroline = zeroline), f)
    close(f)

    .misha_rename_rewrite_meta(meta_file,
        old = c("chr1", "chr2"), new = c("A", "B")
    )

    f <- file(meta_file, "rb")
    m <- unserialize(f)
    close(f)
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
    f <- file(path, "wb")
    serialize(df, f)
    close(f)

    .misha_rename_rewrite_single_interv(path,
        old = c("chr1", "chr2"), new = c("A", "B")
    )

    f <- file(path, "rb")
    out <- unserialize(f)
    close(f)
    expect_equal(as.character(out$chrom), c("A", "B"))
    expect_equal(levels(out$chrom), c("A", "B", "chr3"))
    unlink(path)
})

test_that("C_gdb_rewrite_genome_idx updates names without changing offsets", {
    db_dir <- tempfile("misha-rename-")
    dir.create(db_dir, recursive = TRUE)
    gdb.init_examples()
    orig_groot <- get("GROOT", envir = .misha)

    file.copy(orig_groot, db_dir, recursive = TRUE)
    new_groot <- file.path(db_dir, basename(orig_groot))

    gdb.convert_to_indexed(new_groot, force = TRUE, validate = FALSE)

    idx_path <- file.path(new_groot, "seq", "genome.idx")

    # Parse genome.idx to capture the "before" state (names, offsets, lengths).
    parse_idx <- function(path) {
        con <- file(path, "rb")
        on.exit(close(con))
        magic <- rawToChar(readBin(con, "raw", n = 8L))
        stopifnot(magic == "MISHAIDX")
        readBin(con, "integer", n = 1L, size = 4L) # version
        n <- readBin(con, "integer", n = 1L, size = 4L)
        readBin(con, "integer", n = 2L, size = 4L) # checksum (as 2x int32)
        entries <- vector("list", n)
        for (i in seq_len(n)) {
            chromid <- readBin(con, "integer", n = 1L, size = 4L)
            name_len <- readBin(con, "integer", n = 1L, size = 2L, signed = FALSE)
            nm <- if (name_len > 0L) {
                rawToChar(readBin(con, "raw", n = name_len))
            } else {
                ""
            }
            # offset/length/reserved are uint64 -- read as raw and summarize.
            tail_raw <- readBin(con, "raw", n = 24L)
            entries[[i]] <- list(chromid = chromid, name = nm, tail = tail_raw)
        }
        entries
    }

    before <- parse_idx(idx_path)
    before_names <- vapply(before, `[[`, character(1), "name")
    stopifnot("chr1" %in% before_names)

    .gcall("C_gdb_rewrite_genome_idx", idx_path, "chr1", "chrFOO", .misha_env())

    after <- parse_idx(idx_path)
    after_names <- vapply(after, `[[`, character(1), "name")

    # Names for the renamed contig updated; others untouched.
    expect_true("chrFOO" %in% after_names)
    expect_false("chr1" %in% after_names)
    expect_equal(length(before), length(after))
    # Offsets/lengths/chromids untouched.
    expect_equal(
        vapply(before, `[[`, integer(1), "chromid"),
        vapply(after, `[[`, integer(1), "chromid")
    )
    expect_identical(
        lapply(before, `[[`, "tail"),
        lapply(after, `[[`, "tail")
    )

    # Checksum was recomputed correctly -- gdb.init() validates it.
    expect_silent(gdb.init(new_groot))

    unlink(db_dir, recursive = TRUE)
})

test_that(".misha_rename_build_plan enumerates example DB correctly", {
    gdb.init_examples()
    src <- get("GROOT", envir = .misha)
    db_dir <- tempfile("misha-rename-plan-")
    dir.create(db_dir, recursive = TRUE)
    file.copy(src, db_dir, recursive = TRUE)
    groot <- file.path(db_dir, basename(src))

    plan <- .misha_rename_build_plan(
        groot,
        mapping = data.frame(old = "chr1", new = "chrFOO", stringsAsFactors = FALSE)
    )

    # Per-chromosome DB.
    expect_false(plan$is_indexed)
    expect_true(is.na(plan$genome_idx_path))

    # seq rename: chr1.seq -> chrFOO.seq (and only that).
    expect_equal(basename(plan$seq_renames), "chr1.seq")
    expect_equal(basename(plan$seq_renames_new), "chrFOO.seq")

    # Every track dir in the example DB should have a chr1 -> chrFOO rename planned.
    expect_gt(length(plan$track_dir_renames), 0)
    for (d in names(plan$track_dir_renames)) {
        r <- plan$track_dir_renames[[d]]
        # Old files should all contain "chr1" segment; new ones should replace it.
        expect_true(
            all(basename(r$old) %in% c(
                "chr1",
                grep("^chr1", basename(r$old), value = TRUE)
            )) ||
                any(grepl("chr1", basename(r$old), fixed = TRUE)),
            info = sprintf("track dir %s: unexpected rename entries", d)
        )
    }

    # 2D pair rewriting: rects_track.track should contain both chr1-chr1 ->
    # chrFOO-chrFOO AND chr1-chr2 -> chrFOO-chr2 (directional rename, not just
    # one side).
    rects <- plan$track_dir_renames[[grep("rects_track\\.track$",
        names(plan$track_dir_renames),
        value = TRUE
    )]]
    if (!is.null(rects)) {
        old_basenames <- basename(rects$old)
        new_basenames <- basename(rects$new)
        # Find any pair file starting with "chr1-"
        pair_idx <- grep("^chr1-", old_basenames)
        if (length(pair_idx) > 0) {
            expect_true(all(grepl("^chrFOO-", new_basenames[pair_idx])),
                info = "chr1-* pair files should be renamed to chrFOO-*"
            )
        }
    }

    # Single-file interval rewrite: annotations.interv should be in the list.
    expect_true(any(grepl("annotations\\.interv$", plan$single_interv_rewrites)))

    unlink(db_dir, recursive = TRUE)
})

test_that(".misha_rename_apply_plan renames per-chrom DB (no swap)", {
    gdb.init_examples()
    src <- get("GROOT", envir = .misha)
    db_dir <- tempfile("misha-rename-apply-")
    dir.create(db_dir, recursive = TRUE)
    file.copy(src, db_dir, recursive = TRUE)
    groot <- file.path(db_dir, basename(src))

    mapping <- data.frame(old = "chr1", new = "chrFOO", stringsAsFactors = FALSE)
    plan <- .misha_rename_build_plan(groot, mapping)
    .misha_rename_apply_plan(plan, mapping)

    if (!plan$is_indexed) {
        expect_true(file.exists(file.path(groot, "seq", "chrFOO.seq")))
        expect_false(file.exists(file.path(groot, "seq", "chr1.seq")))
    }

    unlink(db_dir, recursive = TRUE)
})

test_that("two-phase plan succeeds for swap", {
    gdb.init_examples()
    src <- get("GROOT", envir = .misha)
    db_dir <- tempfile("misha-rename-swap-")
    dir.create(db_dir, recursive = TRUE)
    file.copy(src, db_dir, recursive = TRUE)
    groot <- file.path(db_dir, basename(src))

    # chrom_sizes.txt stores stripped names (e.g. "1") while seq/ uses the
    # full names ("chr1.seq"); mapping must match the filesystem.
    seq_files <- list.files(file.path(groot, "seq"), pattern = "\\.seq$")
    chroms <- sub("\\.seq$", "", seq_files)
    if (length(chroms) < 2) skip("example DB has <2 chromosomes")
    a <- chroms[1]
    b <- chroms[2]

    mapping <- data.frame(old = c(a, b), new = c(b, a), stringsAsFactors = FALSE)

    plan <- .misha_rename_build_plan(groot, mapping)
    .misha_rename_apply_with_swap(groot, plan, mapping)

    if (!plan$is_indexed) {
        expect_true(file.exists(file.path(groot, "seq", paste0(b, ".seq"))))
        expect_true(file.exists(file.path(groot, "seq", paste0(a, ".seq"))))
    }

    unlink(db_dir, recursive = TRUE)
})

test_that("gdb.rename_chroms round-trips on per-chromosome DB", {
    gdb.init_examples()
    src <- get("GROOT", envir = .misha)
    db_dir <- tempfile("misha-rename-e2e-")
    dir.create(db_dir, recursive = TRUE)
    file.copy(src, db_dir, recursive = TRUE)
    groot <- file.path(db_dir, basename(src))

    # Use ALLGENOME canonical names (post-init form) which match filesystem.
    gdb.init(groot)
    allgenome <- get("ALLGENOME", envir = .misha)[[1]]
    target <- as.character(allgenome$chrom)[1]

    # Capture chrom_sizes.txt state (in its own form, possibly stripped).
    cs <- read.table(file.path(groot, "chrom_sizes.txt"),
        stringsAsFactors = FALSE, col.names = c("chrom", "size")
    )
    stopifnot(nrow(cs) >= 2)

    mapping <- data.frame(
        old = target, new = "RENAMED",
        stringsAsFactors = FALSE
    )

    gdb.rename_chroms(groot = groot, mapping = mapping, force = TRUE)

    # chrom_sizes.txt must reflect the rename (form-preserving -- may strip
    # the prefix if original form did).
    cs2 <- read.table(file.path(groot, "chrom_sizes.txt"),
        stringsAsFactors = FALSE, col.names = c("chrom", "size")
    )
    # Accept either "RENAMED" or a stripped variant in the new row.
    expect_true("RENAMED" %in% cs2$chrom || sub("^chr", "", "RENAMED") %in% cs2$chrom)

    # Database should re-initialize cleanly.
    expect_silent(gdb.init(groot))

    unlink(db_dir, recursive = TRUE)
})

test_that("gdb.rename_chroms round-trips on indexed DB", {
    gdb.init_examples()
    src <- get("GROOT", envir = .misha)
    db_dir <- tempfile("misha-rename-idx-")
    dir.create(db_dir, recursive = TRUE)
    file.copy(src, db_dir, recursive = TRUE)
    groot <- file.path(db_dir, basename(src))

    gdb.convert_to_indexed(groot,
        force = TRUE, validate = FALSE,
        convert_tracks = TRUE, convert_intervals = TRUE
    )
    gdb.init(groot)

    allgenome <- get("ALLGENOME", envir = .misha)[[1]]
    stopifnot(length(allgenome$chrom) >= 1)
    target <- as.character(allgenome$chrom)[1]

    gdb.rename_chroms(
        groot = groot,
        mapping = data.frame(old = target, new = "FOO", stringsAsFactors = FALSE),
        force = TRUE
    )

    cs2 <- read.table(file.path(groot, "chrom_sizes.txt"),
        stringsAsFactors = FALSE, col.names = c("chrom", "size")
    )
    # Accept prefix or stripped form (form-preserving rewrite).
    expect_true("FOO" %in% cs2$chrom || sub("^chr", "", "FOO") %in% cs2$chrom)
    expect_silent(gdb.init(groot))

    unlink(db_dir, recursive = TRUE)
})

test_that("gdb.rename_chroms handles a swap", {
    gdb.init_examples()
    src <- get("GROOT", envir = .misha)
    db_dir <- tempfile("misha-rename-swap2-")
    dir.create(db_dir, recursive = TRUE)
    file.copy(src, db_dir, recursive = TRUE)
    groot <- file.path(db_dir, basename(src))

    # Use ALLGENOME canonical (filesystem-matching) names.
    gdb.init(groot)
    allgenome <- get("ALLGENOME", envir = .misha)[[1]]
    if (length(allgenome$chrom) < 2) skip("example DB has <2 chromosomes")
    a <- as.character(allgenome$chrom)[1]
    b <- as.character(allgenome$chrom)[2]
    size_a <- as.numeric(allgenome$end[1])
    size_b <- as.numeric(allgenome$end[2])

    gdb.rename_chroms(
        groot = groot,
        mapping = data.frame(
            old = c(a, b), new = c(b, a),
            stringsAsFactors = FALSE
        ),
        force = TRUE
    )

    # After the swap, ALLGENOME should reflect the swapped sizes.
    gdb.init(groot)
    allgenome2 <- get("ALLGENOME", envir = .misha)[[1]]
    names2 <- as.character(allgenome2$chrom)
    ends2 <- as.numeric(allgenome2$end)
    expect_equal(ends2[names2 == b], size_a)
    expect_equal(ends2[names2 == a], size_b)
    expect_silent(gdb.init(groot))

    unlink(db_dir, recursive = TRUE)
})

test_that("gdb.rename_chroms dry_run does not modify any file", {
    gdb.init_examples()
    src <- get("GROOT", envir = .misha)
    db_dir <- tempfile("misha-rename-dry-")
    dir.create(db_dir, recursive = TRUE)
    file.copy(src, db_dir, recursive = TRUE)
    groot <- file.path(db_dir, basename(src))

    gdb.init(groot)
    allgenome <- get("ALLGENOME", envir = .misha)[[1]]
    chrom <- as.character(allgenome$chrom)[1]

    before_files <- sort(list.files(groot, recursive = TRUE, full.names = TRUE))
    before_sizes <- file.info(before_files)$size
    before_mtime <- file.info(before_files)$mtime

    mapping <- data.frame(old = chrom, new = "FOO", stringsAsFactors = FALSE)

    out <- capture.output(
        gdb.rename_chroms(groot = groot, mapping = mapping, dry_run = TRUE)
    )
    expect_true(any(grepl("dry run", out)))

    after_files <- sort(list.files(groot, recursive = TRUE, full.names = TRUE))
    after_sizes <- file.info(after_files)$size
    after_mtime <- file.info(after_files)$mtime

    expect_equal(before_files, after_files)
    expect_equal(before_sizes, after_sizes)
    expect_equal(before_mtime, after_mtime)

    unlink(db_dir, recursive = TRUE)
})

test_that("gdb.rename_chroms re-initializes a loaded DB on completion", {
    gdb.init_examples()
    src <- get("GROOT", envir = .misha)
    db_dir <- tempfile("misha-rename-reload-")
    dir.create(db_dir, recursive = TRUE)
    file.copy(src, db_dir, recursive = TRUE)
    groot <- file.path(db_dir, basename(src))

    # Load the copied DB so that `was_loaded = TRUE` in gdb.rename_chroms.
    gdb.init(groot)
    expect_equal(
        normalizePath(get("GROOT", envir = .misha)),
        normalizePath(groot)
    )

    allgenome <- get("ALLGENOME", envir = .misha)[[1]]
    chrom <- as.character(allgenome$chrom)[1]

    gdb.rename_chroms(
        groot = NULL, # use active DB
        mapping = data.frame(
            old = chrom, new = "RELOADED",
            stringsAsFactors = FALSE
        ),
        force = TRUE
    )

    # The active DB should still be the same path, but ALLGENOME should
    # reflect the rename (ALLGENOME is built from chrom_sizes.txt).
    expect_equal(
        normalizePath(get("GROOT", envir = .misha)),
        normalizePath(groot)
    )
    cs2 <- read.table(file.path(groot, "chrom_sizes.txt"),
        stringsAsFactors = FALSE, col.names = c("chrom", "size")
    )
    # Form-preserving rewrite may strip "chr" prefix -- accept either.
    expect_true("RELOADED" %in% cs2$chrom || sub("^chr", "", "RELOADED") %in% cs2$chrom)

    unlink(db_dir, recursive = TRUE)
})

test_that("gdb.rename_chroms preserves sequence and track data across rename", {
    gdb.init_examples()
    src <- get("GROOT", envir = .misha)
    db_dir <- tempfile("misha-rename-integrity-")
    dir.create(db_dir, recursive = TRUE)
    file.copy(src, db_dir, recursive = TRUE)
    groot <- file.path(db_dir, basename(src))

    gdb.init(groot)
    chrom <- as.character(get("ALLGENOME", envir = .misha)[[1]]$chrom[1])
    size <- as.numeric(get("ALLGENOME", envir = .misha)[[1]]$end[
        match(chrom, get("ALLGENOME", envir = .misha)[[1]]$chrom)
    ])
    probe_end <- min(1000L, as.integer(size))

    # Baseline: pull a sequence window and dense-track values.
    seq_before <- gseq.extract(data.frame(
        chrom = chrom, start = 0L, end = probe_end,
        stringsAsFactors = FALSE
    ))

    tracks <- gtrack.ls()
    dense_tracks <- tracks[tracks %in% c("dense_track", "dense_track2", "sparse_track")]
    track_vals_before <- lapply(dense_tracks, function(tr) {
        tryCatch(
            gextract(tr, data.frame(
                chrom = chrom, start = 0L, end = probe_end,
                stringsAsFactors = FALSE
            )),
            error = function(e) NULL
        )
    })

    # Rename.
    new_name <- "RENAMED_X"
    gdb.rename_chroms(
        groot = groot,
        mapping = data.frame(old = chrom, new = new_name, stringsAsFactors = FALSE),
        force = TRUE
    )

    # Verify sequence payload is identical under the new name.
    gdb.init(groot)
    seq_after <- gseq.extract(data.frame(
        chrom = new_name, start = 0L, end = probe_end,
        stringsAsFactors = FALSE
    ))
    expect_equal(seq_after, seq_before)

    # Verify each track returns the same values (ignoring chrom column which
    # legitimately changes).
    for (i in seq_along(dense_tracks)) {
        tr <- dense_tracks[i]
        before <- track_vals_before[[i]]
        if (is.null(before)) next
        after <- gextract(tr, data.frame(
            chrom = new_name, start = 0L, end = probe_end,
            stringsAsFactors = FALSE
        ))
        value_cols <- setdiff(colnames(before), c("chrom", "start", "end", "intervalID"))
        expect_equal(after[, value_cols, drop = FALSE],
            before[, value_cols, drop = FALSE],
            info = sprintf("track '%s' values changed across rename", tr)
        )
    }

    unlink(db_dir, recursive = TRUE)
})


test_that(".misha_rename_check_writable errors on read-only directory", {
    skip_on_os("windows")
    skip_if(
        Sys.getenv("USER") == "root" || Sys.getenv("USERNAME") == "root",
        "permission checks are bypassed under root"
    )

    d <- tempfile("misha-ro-")
    dir.create(d)
    Sys.chmod(d, "0555")
    on.exit(
        {
            Sys.chmod(d, "0755")
            unlink(d, recursive = TRUE)
        },
        add = TRUE
    )

    expect_error(
        .misha_rename_check_writable(d),
        "No write permission|not writable"
    )
})

test_that(".misha_rename_check_writable accepts writable directory", {
    d <- tempfile("misha-rw-")
    dir.create(d)
    on.exit(unlink(d, recursive = TRUE), add = TRUE)
    expect_silent(.misha_rename_check_writable(d))
})

test_that("gdb.rename_chroms bails pre-flight when a dir is read-only, leaving DB intact", {
    skip_on_os("windows")
    skip_if(
        Sys.getenv("USER") == "root" || Sys.getenv("USERNAME") == "root",
        "permission checks are bypassed under root"
    )

    gdb.init_examples()
    src <- get("GROOT", envir = .misha)
    db_dir <- tempfile("misha-rename-ro-")
    dir.create(db_dir, recursive = TRUE)
    file.copy(src, db_dir, recursive = TRUE)
    groot <- file.path(db_dir, basename(src))

    gdb.init(groot)
    chrom <- as.character(get("ALLGENOME", envir = .misha)[[1]]$chrom[1])

    # Find a directory holding a single-file .interv and make it read-only.
    intervs <- list.files(file.path(groot, "tracks"),
        pattern = "\\.interv$", recursive = TRUE, full.names = TRUE,
        include.dirs = TRUE
    )
    single <- intervs[!(file.info(intervs)$isdir %in% TRUE)]
    if (length(single) == 0) skip("no single-file .interv in example DB")
    ro_dir <- dirname(single[1])
    Sys.chmod(ro_dir, "0555")
    on.exit(
        {
            Sys.chmod(ro_dir, "0755")
            unlink(db_dir, recursive = TRUE)
        },
        add = TRUE
    )

    cs_before <- readLines(file.path(groot, "chrom_sizes.txt"))
    idx_mtime_before <- file.info(file.path(groot, "seq", "genome.idx"))$mtime

    expect_error(
        gdb.rename_chroms(
            groot = groot,
            mapping = data.frame(old = chrom, new = "RENAMED", stringsAsFactors = FALSE),
            force = TRUE
        ),
        "write permission"
    )

    # DB untouched: chrom_sizes.txt unchanged, genome.idx mtime unchanged,
    # no breadcrumb.
    expect_equal(readLines(file.path(groot, "chrom_sizes.txt")), cs_before)
    idx_mtime_after <- file.info(file.path(groot, "seq", "genome.idx"))$mtime
    if (!is.na(idx_mtime_before) && !is.na(idx_mtime_after)) {
        expect_equal(idx_mtime_after, idx_mtime_before)
    }
    expect_false(file.exists(file.path(groot, ".rename_interrupted")))
})

test_that("gdb.rename_chroms preserves 2D rectangles-track values across rename", {
    gdb.init_examples()
    src <- get("GROOT", envir = .misha)
    db_dir <- tempfile("misha-rename-2d-")
    dir.create(db_dir, recursive = TRUE)
    file.copy(src, db_dir, recursive = TRUE)
    groot <- file.path(db_dir, basename(src))

    gdb.init(groot)
    chrom <- as.character(get("ALLGENOME", envir = .misha)[[1]]$chrom[1])

    before <- gextract("rects_track", gintervals.2d(chrom, chroms2 = chrom),
        iterator = "rects_track"
    )

    gdb.rename_chroms(
        groot = groot,
        mapping = data.frame(old = chrom, new = "RENAMED_2D", stringsAsFactors = FALSE),
        force = TRUE
    )
    gdb.init(groot)

    after <- gextract("rects_track", gintervals.2d("RENAMED_2D", chroms2 = "RENAMED_2D"),
        iterator = "rects_track"
    )

    value_cols <- setdiff(
        colnames(before),
        c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "intervalID")
    )
    expect_equal(nrow(after), nrow(before))
    if (length(value_cols) > 0 && nrow(before) > 0) {
        expect_equal(after[, value_cols, drop = FALSE], before[, value_cols, drop = FALSE])
    }

    unlink(db_dir, recursive = TRUE)
})

test_that("gdb.rename_chroms partial rename leaves unmapped chroms untouched", {
    gdb.init_examples()
    src <- get("GROOT", envir = .misha)
    db_dir <- tempfile("misha-rename-partial-")
    dir.create(db_dir, recursive = TRUE)
    file.copy(src, db_dir, recursive = TRUE)
    groot <- file.path(db_dir, basename(src))

    # Convert to indexed so chrom_sizes.txt and ALLGENOME agree on the
    # prefixed form -- avoids misha's per-chrom auto-prefix heuristic, which
    # the example DB's mixed form triggers.
    gdb.convert_to_indexed(groot,
        force = TRUE, validate = FALSE,
        convert_tracks = TRUE, convert_intervals = TRUE
    )
    gdb.init(groot)
    chroms <- as.character(get("ALLGENOME", envir = .misha)[[1]]$chrom)
    if (length(chroms) < 2) skip("example DB has <2 chromosomes")
    keep <- chroms[2]
    renamed <- chroms[1]

    # Snapshot a sequence of the un-mapped chrom; it must be byte-identical
    # afterwards.
    keep_size <- get("ALLGENOME", envir = .misha)[[1]]$end[
        match(keep, get("ALLGENOME", envir = .misha)[[1]]$chrom)
    ]
    probe_end <- min(200L, as.integer(keep_size))
    seq_keep_before <- gseq.extract(data.frame(
        chrom = keep, start = 0L, end = probe_end, stringsAsFactors = FALSE
    ))

    gdb.rename_chroms(
        groot = groot,
        mapping = data.frame(old = renamed, new = "ONLY_ONE", stringsAsFactors = FALSE),
        force = TRUE
    )
    gdb.init(groot)

    new_chroms <- as.character(get("ALLGENOME", envir = .misha)[[1]]$chrom)
    expect_true(keep %in% new_chroms)
    expect_true("ONLY_ONE" %in% new_chroms)
    expect_false(renamed %in% new_chroms)

    # Un-mapped chrom data is unchanged.
    seq_keep_after <- gseq.extract(data.frame(
        chrom = keep, start = 0L, end = probe_end, stringsAsFactors = FALSE
    ))
    expect_equal(seq_keep_after, seq_keep_before)

    unlink(db_dir, recursive = TRUE)
})

test_that("gdb.rename_chroms removes breadcrumb on success, blocks re-run without force", {
    gdb.init_examples()
    src <- get("GROOT", envir = .misha)
    db_dir <- tempfile("misha-rename-crumb-")
    dir.create(db_dir, recursive = TRUE)
    file.copy(src, db_dir, recursive = TRUE)
    groot <- file.path(db_dir, basename(src))

    gdb.init(groot)
    chrom <- as.character(get("ALLGENOME", envir = .misha)[[1]]$chrom[1])

    gdb.rename_chroms(
        groot = groot,
        mapping = data.frame(old = chrom, new = "CRUMB_1", stringsAsFactors = FALSE),
        force = TRUE
    )
    # After success, the breadcrumb must be gone.
    expect_false(file.exists(file.path(groot, ".rename_interrupted")))

    # Now plant a stale breadcrumb and verify the function refuses.
    writeLines("stale", file.path(groot, ".rename_interrupted"))
    expect_error(
        gdb.rename_chroms(
            groot = groot,
            mapping = data.frame(old = "CRUMB_1", new = "CRUMB_2", stringsAsFactors = FALSE)
        ),
        "interrupted rename"
    )

    # With force = TRUE the function proceeds and clears the breadcrumb.
    gdb.rename_chroms(
        groot = groot,
        mapping = data.frame(old = "CRUMB_1", new = "CRUMB_2", stringsAsFactors = FALSE),
        force = TRUE
    )
    expect_false(file.exists(file.path(groot, ".rename_interrupted")))

    unlink(db_dir, recursive = TRUE)
})

test_that("gdb.rename_chroms round-trips a rename and its inverse back to the original", {
    gdb.init_examples()
    src <- get("GROOT", envir = .misha)
    db_dir <- tempfile("misha-rename-roundtrip-")
    dir.create(db_dir, recursive = TRUE)
    file.copy(src, db_dir, recursive = TRUE)
    groot <- file.path(db_dir, basename(src))

    # Use indexed format: chrom_sizes.txt and ALLGENOME agree on the prefixed
    # form, making the A -> B -> A round-trip exactly byte-reversible.
    gdb.convert_to_indexed(groot,
        force = TRUE, validate = FALSE,
        convert_tracks = TRUE, convert_intervals = TRUE
    )
    gdb.init(groot)
    chrom <- as.character(get("ALLGENOME", envir = .misha)[[1]]$chrom[1])

    cs_before <- readLines(file.path(groot, "chrom_sizes.txt"))
    seq_dat_sum_before <- tools::md5sum(file.path(groot, "seq", "genome.seq"))

    # Forward rename.
    gdb.rename_chroms(
        groot = groot,
        mapping = data.frame(old = chrom, new = "TMP_NAME", stringsAsFactors = FALSE),
        force = TRUE
    )
    # Reverse.
    gdb.rename_chroms(
        groot = groot,
        mapping = data.frame(old = "TMP_NAME", new = chrom, stringsAsFactors = FALSE),
        force = TRUE
    )

    # chrom_sizes.txt is byte-identical.
    expect_equal(readLines(file.path(groot, "chrom_sizes.txt")), cs_before)
    # Payload files untouched.
    expect_equal(tools::md5sum(file.path(groot, "seq", "genome.seq")), seq_dat_sum_before)

    # DB reloads cleanly with original chroms.
    gdb.init(groot)
    expect_true(chrom %in% as.character(get("ALLGENOME", envir = .misha)[[1]]$chrom))

    unlink(db_dir, recursive = TRUE)
})

test_that("gdb.rename_chroms supports renaming many chromosomes in one call", {
    gdb.init_examples()
    src <- get("GROOT", envir = .misha)
    db_dir <- tempfile("misha-rename-bulk-")
    dir.create(db_dir, recursive = TRUE)
    file.copy(src, db_dir, recursive = TRUE)
    groot <- file.path(db_dir, basename(src))

    gdb.init(groot)
    chroms <- as.character(get("ALLGENOME", envir = .misha)[[1]]$chrom)
    new_names <- paste0("R", seq_along(chroms))

    # Baseline extraction for each chrom.
    before <- lapply(chroms, function(cc) {
        size <- get("ALLGENOME", envir = .misha)[[1]]$end[
            match(cc, get("ALLGENOME", envir = .misha)[[1]]$chrom)
        ]
        gseq.extract(data.frame(
            chrom = cc, start = 0L, end = min(50L, as.integer(size)),
            stringsAsFactors = FALSE
        ))
    })

    gdb.rename_chroms(
        groot = groot,
        mapping = data.frame(old = chroms, new = new_names, stringsAsFactors = FALSE),
        force = TRUE
    )
    gdb.init(groot)

    new_chroms_in_db <- as.character(get("ALLGENOME", envir = .misha)[[1]]$chrom)
    for (nn in new_names) expect_true(nn %in% new_chroms_in_db)
    for (old in chroms) expect_false(old %in% new_chroms_in_db)

    # Sequences come out identical under the new names.
    for (i in seq_along(chroms)) {
        size <- get("ALLGENOME", envir = .misha)[[1]]$end[
            match(new_names[i], get("ALLGENOME", envir = .misha)[[1]]$chrom)
        ]
        after <- gseq.extract(data.frame(
            chrom = new_names[i], start = 0L, end = min(50L, as.integer(size)),
            stringsAsFactors = FALSE
        ))
        expect_equal(after, before[[i]],
            info = sprintf("%s -> %s", chroms[i], new_names[i])
        )
    }

    unlink(db_dir, recursive = TRUE)
})

test_that("gdb.rename_chroms dry_run prints a plan summary with all affected categories", {
    gdb.init_examples()
    src <- get("GROOT", envir = .misha)
    db_dir <- tempfile("misha-rename-dry2-")
    dir.create(db_dir, recursive = TRUE)
    file.copy(src, db_dir, recursive = TRUE)
    groot <- file.path(db_dir, basename(src))

    gdb.init(groot)
    chrom <- as.character(get("ALLGENOME", envir = .misha)[[1]]$chrom[1])

    out <- capture.output(
        gdb.rename_chroms(
            groot = groot,
            mapping = data.frame(old = chrom, new = "DRY", stringsAsFactors = FALSE),
            dry_run = TRUE
        )
    )
    blob <- paste(out, collapse = "\n")

    expect_match(blob, "Database: ", fixed = TRUE)
    expect_match(blob, "Format: per-chromosome")
    expect_match(blob, "sequence")
    expect_match(blob, "single-file \\.interv")
    expect_match(blob, "Chromosomes affected: 1")
    expect_match(blob, "dry run")

    unlink(db_dir, recursive = TRUE)
})
