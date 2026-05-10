test_that(".compute_chrom_aliases handles canonical chr-prefixed names", {
    map <- .compute_chrom_aliases(c("chr1", "chr2", "chrX"))
    # Each canonical name maps to itself, plus an unprefixed alias for each
    expect_equal(map[["chr1"]], "chr1")
    expect_equal(map[["chr2"]], "chr2")
    expect_equal(map[["chrX"]], "chrX")
    expect_equal(map[["1"]], "chr1")
    expect_equal(map[["2"]], "chr2")
    expect_equal(map[["X"]], "chrX")
    expect_setequal(names(map), c("chr1", "chr2", "chrX", "1", "2", "X"))
})

test_that(".compute_chrom_aliases handles canonical unprefixed names", {
    map <- .compute_chrom_aliases(c("1", "2", "X"))
    expect_equal(map[["1"]], "1")
    expect_equal(map[["chr1"]], "1")
    expect_equal(map[["X"]], "X")
    expect_equal(map[["chrX"]], "X")
    expect_setequal(names(map), c("1", "2", "X", "chr1", "chr2", "chrX"))
})

test_that(".compute_chrom_aliases preserves canonical when alias collides", {
    # When both "chr1" and "1" are canonical, neither becomes an alias of the other
    map <- .compute_chrom_aliases(c("chr1", "1"))
    expect_equal(map[["chr1"]], "chr1")
    expect_equal(map[["1"]], "1")
    expect_setequal(names(map), c("chr1", "1"))
})

test_that(".compute_chrom_aliases first-seen-wins for non-canonical alias collisions", {
    # No canonical name is "1", but two distinct canonical names both want "1" as an alias.
    # The first one in input order wins (matching the original loop semantics).
    map <- .compute_chrom_aliases(c("chr1", "chr1_alt"))
    expect_equal(map[["chr1"]], "chr1")
    expect_equal(map[["chr1_alt"]], "chr1_alt")
    expect_equal(map[["1"]], "chr1") # first-seen-wins
    expect_true("1_alt" %in% names(map))
    expect_equal(map[["1_alt"]], "chr1_alt")
})

test_that(".compute_chrom_aliases handles mitochondrial canonical chrM", {
    map <- .compute_chrom_aliases(c("chr1", "chrM"))
    expect_equal(map[["chrM"]], "chrM")
    expect_equal(map[["M"]], "chrM")
    expect_equal(map[["MT"]], "chrM")
})

test_that(".compute_chrom_aliases handles mitochondrial canonical MT", {
    map <- .compute_chrom_aliases(c("chr1", "MT"))
    expect_equal(map[["MT"]], "MT")
    expect_equal(map[["M"]], "MT")
    expect_equal(map[["chrM"]], "MT")
})

test_that(".compute_chrom_aliases handles mitochondrial canonical M", {
    map <- .compute_chrom_aliases(c("chr1", "M"))
    expect_equal(map[["M"]], "M")
    expect_equal(map[["MT"]], "M")
    expect_equal(map[["chrM"]], "M")
})

test_that(".compute_chrom_aliases mitochondrial does not overwrite other canonical", {
    # If both M and MT are canonical, both stay canonical
    map <- .compute_chrom_aliases(c("M", "MT"))
    expect_equal(map[["M"]], "M")
    expect_equal(map[["MT"]], "MT")
    # chrM is added once, pointing to whichever mito chrom appeared first
    expect_equal(map[["chrM"]], "M")
})

test_that(".compute_chrom_aliases handles empty input", {
    map <- .compute_chrom_aliases(character(0))
    expect_equal(length(map), 0)
})

test_that(".compute_chrom_aliases dedups duplicate canonical names", {
    map <- .compute_chrom_aliases(c("chr1", "chr1", "chr2"))
    # Only unique chroms should appear; aliases are correct
    expect_equal(map[["chr1"]], "chr1")
    expect_equal(map[["1"]], "chr1")
    expect_equal(map[["chr2"]], "chr2")
    expect_equal(map[["2"]], "chr2")
})

test_that(".compute_chrom_aliases keeps prefixed aliases for small non-chr non-mito assemblies", {
    # Small contig sets stay below the skip threshold so the Ensembl-style
    # "1"/"chr1" toggle keeps working for any small genome.
    chroms <- c("scaffold_1", "scaffold_2", "contig_42")
    map <- .compute_chrom_aliases(chroms)
    expect_equal(map[["scaffold_1"]], "scaffold_1")
    expect_equal(map[["chrscaffold_1"]], "scaffold_1")
    expect_equal(map[["chrcontig_42"]], "contig_42")
})

test_that(".compute_chrom_aliases skips chr-prefixed aliases on fragmented assemblies", {
    # Phylo447-style: thousands of "scaffold_*" contigs with no chr prefix and no
    # mito hits. Prefixed aliases are pure waste here -- nobody types
    # "chrscaffold_1234567" -- so the alias map should contain canonicals only.
    n <- 1500L
    chroms <- paste0("scaffold_", seq_len(n))
    map <- .compute_chrom_aliases(chroms)
    expect_setequal(names(map), chroms)
    expect_equal(unname(map[chroms]), chroms)
})

test_that(".compute_chrom_aliases keeps prefixed aliases on large Ensembl-style genomes", {
    # Hypothetical 1500-chrom Ensembl-style genome with "MT". Mito presence keeps
    # the prefixed-alias path enabled, so all "chr*" toggles still work.
    chroms <- c(as.character(1:1499), "MT")
    map <- .compute_chrom_aliases(chroms)
    expect_equal(map[["1"]], "1")
    expect_equal(map[["chr1"]], "1")
    expect_equal(map[["MT"]], "MT")
    expect_equal(map[["chrM"]], "MT")
})

test_that(".compute_chrom_aliases scales to large contig counts", {
    # Phylo447-style assembly: 50k contigs should complete in well under 1s.
    # Map should contain only canonical names (no chr-prefixed aliases) since
    # none of the inputs have chr prefix and none are mito-pattern.
    n <- 50000L
    chroms <- paste0("scaffold_", seq_len(n))
    t0 <- Sys.time()
    map <- .compute_chrom_aliases(chroms)
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    expect_lt(elapsed, 2.0)
    expect_equal(length(map), n)
})
