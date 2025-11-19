load_test_db()

# Tests for best_source_cluster policy
# This policy clusters mappings by source overlap:
# - If source chains overlap (duplication): keep all mappings
# - If source chains are disjoint (conflict): keep only the cluster with largest total target length

test_that("best_source_cluster policy retains duplications (overlapping source chains)", {
    local_db_state()

    # Create target genome with enough space
    setup_db(list(">chr1\n", paste(rep("A", 1000), collapse = ""), "\n"))

    chain_file <- new_chain_file()

    # Scenario: Duplication (Chain A and B overlap in source [0-100])
    # Chain A: Source src[0-100] -> Target chr1[0-100] (score 100)
    write_chain_entry(chain_file, "src", 1000, "+", 0, 100, "chr1", 1000, "+", 0, 100, 100)
    # Chain B: Source src[0-100] -> Target chr1[200-300] (score 50, same source region = duplication)
    write_chain_entry(chain_file, "src", 1000, "+", 0, 100, "chr1", 1000, "+", 200, 300, 50)

    # Load chain with keep policy to allow overlapping sources, and best_source_cluster for target policy
    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "best_source_cluster")

    # Input interval that maps through both chains
    intervals <- data.frame(chrom = "src", start = 0, end = 100)

    # With best_source_cluster: both mappings should be retained (overlapping source = duplication)
    result <- gintervals.liftover(intervals, chain)

    expect_equal(nrow(result), 2, info = "Should retain both mappings for overlapping source chains")
    expect_true(0 %in% result$start, info = "Should have mapping to [0-100]")
    expect_true(200 %in% result$start, info = "Should have mapping to [200-300]")
})

test_that("best_source_cluster policy selects best cluster for disjoint source chains", {
    local_db_state()

    # Create target genome
    setup_db(list(">chr1\n", paste(rep("A", 1000), collapse = ""), "\n"))

    chain_file <- new_chain_file()

    # Scenario: Disjoint source chains (conflict)
    # Cluster 1: Source src[0-50] -> Target chr1[0-50] (length 50)
    write_chain_entry(chain_file, "src", 1000, "+", 0, 50, "chr1", 1000, "+", 0, 50, 100)
    # Cluster 2: Source src[100-200] -> Target chr1[200-300] (length 100) - larger mass
    write_chain_entry(chain_file, "src", 1000, "+", 100, 200, "chr1", 1000, "+", 200, 300, 50)

    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "best_source_cluster")

    # Input interval that spans both chain regions (0-200)
    intervals <- data.frame(chrom = "src", start = 0, end = 200)

    # With best_source_cluster: should keep only the cluster with largest total target length
    # Cluster 1: maps [0-50] -> length 50
    # Cluster 2: maps [100-200] -> length 100 (WINS)
    result <- gintervals.liftover(intervals, chain)

    expect_equal(nrow(result), 1, info = "Should keep only the heaviest cluster")
    expect_equal(result$start, 200, info = "Should keep the cluster with larger mass")
    expect_equal(result$end, 300, info = "Should have correct end coordinate")
})

test_that("best_source_cluster policy handles multiple overlapping chains in winning cluster", {
    local_db_state()

    # Create target genome
    setup_db(list(">chr1\n", paste(rep("A", 1000), collapse = ""), "\n"))

    chain_file <- new_chain_file()

    # Cluster 1 (overlapping sources): total mass = 100 + 100 = 200
    # Chain A: src[0-100] -> chr1[0-100]
    write_chain_entry(chain_file, "src", 1000, "+", 0, 100, "chr1", 1000, "+", 0, 100, 100)
    # Chain B: src[50-150] -> chr1[200-300] (overlaps A in source)
    write_chain_entry(chain_file, "src", 1000, "+", 50, 150, "chr1", 1000, "+", 200, 300, 90)

    # Cluster 2 (disjoint): total mass = 50
    # Chain C: src[500-550] -> chr1[500-550]
    write_chain_entry(chain_file, "src", 1000, "+", 500, 550, "chr1", 1000, "+", 500, 550, 80)

    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "best_source_cluster")

    # Input interval spanning all
    intervals <- data.frame(chrom = "src", start = 0, end = 600)

    result <- gintervals.liftover(intervals, chain)

    # Cluster 1 should win (mass 200 > 50), and both chains in cluster should be retained
    expect_equal(nrow(result), 2, info = "Should retain both chains in the winning cluster")
    expect_true(all(result$start %in% c(0, 200)), info = "Should have both mappings from cluster 1")
})

test_that("best_source_cluster policy works with single mapping", {
    local_db_state()

    setup_db(list(">chr1\n", paste(rep("A", 1000), collapse = ""), "\n"))

    chain_file <- new_chain_file()

    # Single chain
    write_chain_entry(chain_file, "src", 1000, "+", 0, 100, "chr1", 1000, "+", 0, 100, 100)

    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "best_source_cluster")

    intervals <- data.frame(chrom = "src", start = 0, end = 100)

    result <- gintervals.liftover(intervals, chain)

    expect_equal(nrow(result), 1)
    expect_equal(result$start, 0)
    expect_equal(result$end, 100)
})

test_that("best_source_cluster policy handles adjacent but non-overlapping sources correctly", {
    local_db_state()

    setup_db(list(">chr1\n", paste(rep("A", 1000), collapse = ""), "\n"))

    chain_file <- new_chain_file()

    # Two chains with adjacent (touching) source regions - NOT overlapping
    # Chain A: src[0-100] -> chr1[0-100]
    write_chain_entry(chain_file, "src", 1000, "+", 0, 100, "chr1", 1000, "+", 0, 100, 100)
    # Chain B: src[100-200] -> chr1[200-300] (starts exactly where A ends)
    write_chain_entry(chain_file, "src", 1000, "+", 100, 200, "chr1", 1000, "+", 200, 300, 100)

    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "best_source_cluster")

    intervals <- data.frame(chrom = "src", start = 0, end = 200)

    result <- gintervals.liftover(intervals, chain)

    # Adjacent but non-overlapping = two separate clusters with equal mass
    # Should keep only one (the first one processed after sorting by start_src)
    expect_equal(nrow(result), 1, info = "Adjacent non-overlapping sources should be separate clusters")
})

test_that("best_source_cluster policy correctly handles mass tie-breaker", {
    local_db_state()

    setup_db(list(">chr1\n", paste(rep("A", 1000), collapse = ""), "\n"))

    chain_file <- new_chain_file()

    # Two clusters with equal mass - first one (by source position) should win
    # Cluster 1: src[0-100] -> chr1[0-100] (mass 100)
    write_chain_entry(chain_file, "src", 1000, "+", 0, 100, "chr1", 1000, "+", 0, 100, 200)
    # Cluster 2: src[500-600] -> chr1[500-600] (mass 100)
    write_chain_entry(chain_file, "src", 1000, "+", 500, 600, "chr1", 1000, "+", 500, 600, 100)

    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "best_source_cluster")

    intervals <- data.frame(chrom = "src", start = 0, end = 700)

    result <- gintervals.liftover(intervals, chain)

    # Equal mass - first cluster should win
    expect_equal(nrow(result), 1)
    # First cluster maps to 0-100
    expect_equal(result$start, 0)
})

# test_that("best_source_cluster policy preserves chain_id and score in output", {
#     local_db_state()

#     setup_db(list(">chr1\n", paste(rep("A", 1000), collapse = ""), "\n"))

#     chain_file <- new_chain_file()

#     # Overlapping sources to get multiple results
#     write_chain_entry(chain_file, "src", 1000, "+", 0, 100, "chr1", 1000, "+", 0, 100, 500)
#     write_chain_entry(chain_file, "src", 1000, "+", 0, 100, "chr1", 1000, "+", 200, 300, 300)

#     chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "best_source_cluster")

#     intervals <- data.frame(chrom = "src", start = 0, end = 100)

#     result <- gintervals.liftover(intervals, chain, include_metadata = TRUE)

#     expect_equal(nrow(result), 2)
#     expect_true("chain_id" %in% names(result))
#     expect_true("score" %in% names(result))
#     expect_true(all(result$score %in% c(500, 300)))
# })

test_that("best_source_cluster policy works with intervals that don't span full chain regions", {
    local_db_state()

    setup_db(list(">chr1\n", paste(rep("A", 1000), collapse = ""), "\n"))

    chain_file <- new_chain_file()

    # Two disjoint chains
    # Chain A: src[0-100] -> chr1[0-100]
    write_chain_entry(chain_file, "src", 1000, "+", 0, 100, "chr1", 1000, "+", 0, 100, 100)
    # Chain B: src[200-400] -> chr1[300-500]
    write_chain_entry(chain_file, "src", 1000, "+", 200, 400, "chr1", 1000, "+", 300, 500, 100)

    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "best_source_cluster")

    # Input interval only partially overlaps each chain
    # Maps through A: [50-100] -> chr1[50-100] (length 50)
    # Maps through B: [200-300] -> chr1[300-400] (length 100) - WINS
    intervals <- data.frame(chrom = "src", start = 50, end = 300)

    result <- gintervals.liftover(intervals, chain)

    expect_equal(nrow(result), 1)
    # Cluster 2 should win with larger mass
    expect_equal(result$start, 300)
    expect_equal(result$end, 400)
})

test_that("best_source_cluster scenario 1: duplication with same source to different chromosomes", {
    local_db_state()

    # Create target genome with chr1 and chr2
    setup_db(list(">chr1\n", paste(rep("A", 1000), collapse = ""), "\n", ">chr2\n", paste(rep("A", 1000), collapse = ""), "\n"))

    chain_file <- new_chain_file()

    # Scenario 1: Same source interval maps to different chromosomes (duplication)
    # 100-200 -> chr1 100-200
    write_chain_entry(chain_file, "src", 1000, "+", 100, 200, "chr1", 1000, "+", 100, 200, 100)
    # 100-200 -> chr2 100-200
    write_chain_entry(chain_file, "src", 1000, "+", 100, 200, "chr2", 1000, "+", 100, 200, 100)

    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "best_source_cluster")

    intervals <- data.frame(chrom = "src", start = 100, end = 200)

    result <- gintervals.liftover(intervals, chain)

    # Result: retain both (duplication - overlapping source chains)
    expect_equal(nrow(result), 2, info = "Should retain both mappings for duplication")
    expect_true(any(result$chrom == "chr1" & result$start == 100 & result$end == 200),
        info = "Should have mapping to chr1 100-200"
    )
    expect_true(any(result$chrom == "chr2" & result$start == 100 & result$end == 200),
        info = "Should have mapping to chr2 100-200"
    )
})

test_that("best_source_cluster scenario 2: adjacent disjoint chains choose higher mass", {
    local_db_state()

    # Create target genome with chr1 and chr2
    setup_db(list(">chr1\n", paste(rep("A", 1000), collapse = ""), "\n", ">chr2\n", paste(rep("A", 1000), collapse = ""), "\n"))

    chain_file <- new_chain_file()

    # Scenario 2: Adjacent but disjoint source chains
    # 100-130 -> chr1 100-130 (mass = 30)
    write_chain_entry(chain_file, "src", 1000, "+", 100, 130, "chr1", 1000, "+", 100, 130, 100)
    # 130-200 -> chr2 130-200 (mass = 70, higher)
    write_chain_entry(chain_file, "src", 1000, "+", 130, 200, "chr2", 1000, "+", 130, 200, 100)

    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "best_source_cluster")

    intervals <- data.frame(chrom = "src", start = 100, end = 200)

    result <- gintervals.liftover(intervals, chain)

    # Result: choose the one with higher mass (130-200 has mass 70 > 30)
    expect_equal(nrow(result), 1, info = "Should keep only the cluster with higher mass")
    expect_equal(as.character(result$chrom), "chr2", info = "Should keep chr2 mapping")
    expect_equal(result$start, 130, info = "Should have correct start coordinate")
    expect_equal(result$end, 200, info = "Should have correct end coordinate")
})

test_that("best_source_cluster scenario 3: overlapping source chains retain both", {
    local_db_state()

    # Create target genome with chr1 and chr2
    setup_db(list(">chr1\n", paste(rep("A", 1000), collapse = ""), "\n", ">chr2\n", paste(rep("A", 1000), collapse = ""), "\n"))

    chain_file <- new_chain_file()

    # Scenario 3: Overlapping source chains
    # 100-200 -> chr1 100-200
    write_chain_entry(chain_file, "src", 1000, "+", 100, 200, "chr1", 1000, "+", 100, 200, 100)
    # 120-220 -> chr2 100-200 (source overlaps with first chain in 120-200)
    write_chain_entry(chain_file, "src", 1000, "+", 120, 220, "chr2", 1000, "+", 100, 200, 100)

    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "best_source_cluster")

    intervals <- data.frame(chrom = "src", start = 100, end = 200)

    result <- gintervals.liftover(intervals, chain)

    # Result: there is overlap (100-200 overlaps with 100-220 in 100-200), so retain both
    expect_equal(nrow(result), 2, info = "Should retain both mappings for overlapping source chains")
    expect_equal(as.character(result$chrom), c("chr1", "chr2"), info = "Should have mappings to both chromosomes")
    expect_equal(result$start, c(100, 100), info = "Should have correct start coordinates")
    expect_equal(result$end, c(200, 180), info = "Should have correct end coordinates")
})

test_that("best_source_cluster scenario 4: Cluster A+B (overlapping) outweighs C (disjoint)", {
    local_db_state()

    # Setup genome
    setup_db(list(">chr1\n", paste(rep("A", 2000), collapse = ""), "\n"))
    chain_file <- new_chain_file()

    # --- Cluster 1: Chains A and B overlap in source ---
    # Chain A: Source 100-200 (Length 100)
    write_chain_entry(chain_file, "src", 2000, "+", 100, 200, "chr1", 2000, "+", 1000, 1100, 1)

    # Chain B: Source 150-250 (Length 100)
    # Overlaps A in source [150-200].
    # Total Mass (A+B) = 100 + 100 = 200.
    write_chain_entry(chain_file, "src", 2000, "+", 150, 250, "chr1", 2000, "+", 1200, 1300, 2)

    # --- Cluster 2: Chain C is disjoint ---
    # Chain C: Source 500-650 (Length 150)
    # Disjoint from A and B.
    # Total Mass (C) = 150.
    write_chain_entry(chain_file, "src", 2000, "+", 500, 650, "chr1", 2000, "+", 1500, 1650, 3)

    # Load chain
    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "best_source_cluster")

    # Query spans all of them
    intervals <- data.frame(chrom = "src", start = 0, end = 1000)

    result <- gintervals.liftover(intervals, chain)

    # Expectation:
    # Cluster 1 (Mass 200) > Cluster 2 (Mass 150).
    # Result should contain Chain A (ID 1) and Chain B (ID 2).
    # Chain C (ID 3) should be discarded.

    expect_equal(nrow(result), 2, info = "Should keep A and B (cluster mass 200) over C (mass 150)")

    # Verify IDs
    result_ids <- sort(result$chain_id)
    expect_equal(result_ids, c(1, 2))

    # Verify coordinates
    expect_true(any(result$start == 1000), info = "Chain A present")
    expect_true(any(result$start == 1200), info = "Chain B present")
    expect_false(any(result$start == 1500), info = "Chain C discarded")
})

test_that("best_source_cluster scenario 5: C (disjoint) outweighs A+B (overlapping)", {
    local_db_state()

    # Setup genome
    setup_db(list(">chr1\n", paste(rep("A", 3000), collapse = ""), "\n"))
    chain_file <- new_chain_file()

    # --- Cluster 1: Chains A and B overlap in source ---
    # Total Mass = 200 (100 + 100)

    # Chain A: Source 100-200 (Length 100)
    write_chain_entry(chain_file, "src", 3000, "+", 100, 200, "chr1", 3000, "+", 1000, 1100, 1)

    # Chain B: Source 150-250 (Length 100)
    # Overlaps A in source [150-200].
    write_chain_entry(chain_file, "src", 3000, "+", 150, 250, "chr1", 3000, "+", 1200, 1300, 2)

    # --- Cluster 2: Chain C is disjoint and larger ---
    # Mass = 250

    # Chain C: Source 500-750 (Length 250)
    # Disjoint from A and B.
    write_chain_entry(chain_file, "src", 3000, "+", 500, 750, "chr1", 3000, "+", 2000, 2250, 3)

    # Load chain
    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "best_source_cluster")

    # Query spans all of them
    intervals <- data.frame(chrom = "src", start = 0, end = 1000)

    result <- gintervals.liftover(intervals, chain)

    # Expectation:
    # Cluster 2 (Mass 250) > Cluster 1 (Mass 200).
    # Result should contain ONLY Chain C (ID 3).
    # Chains A and B should be discarded.

    expect_equal(nrow(result), 1, info = "Should keep C (mass 250) over A+B (mass 200)")
    expect_equal(result$chain_id, 3, info = "Should be chain C")
    expect_equal(result$start, 2000, info = "Should be at chain C target start")
    expect_equal(result$end, 2250, info = "Should be at chain C target end")
})

test_that("best_source_cluster scenario 6: Transitive overlap (Bridge) A-B-C outweighs D", {
    local_db_state()

    # Setup genome
    setup_db(list(">chr1\n", paste(rep("A", 3000), collapse = ""), "\n"))
    chain_file <- new_chain_file()

    # --- Cluster 1: Transitive Chain A -> B -> C ---

    # Chain A: Source 100-200 (Target Len 100)
    write_chain_entry(chain_file, "src", 3000, "+", 100, 200, "chr1", 3000, "+", 1000, 1100, 1)

    # Chain B: Source 150-250 (Target Len 100)
    # Overlaps A (150-200)
    write_chain_entry(chain_file, "src", 3000, "+", 150, 250, "chr1", 3000, "+", 1200, 1300, 2)

    # Chain C: Source 240-340 (Target Len 100)
    # Overlaps B (240-250), but disjoint from A (ends at 200)
    write_chain_entry(chain_file, "src", 3000, "+", 240, 340, "chr1", 3000, "+", 1400, 1500, 3)

    # Total Mass of Cluster A-B-C = 300

    # --- Cluster 2: Disjoint Distractor D ---
    # Chain D: Source 500-750 (Target Len 250)
    # Mass = 250
    write_chain_entry(chain_file, "src", 3000, "+", 500, 750, "chr1", 3000, "+", 2000, 2250, 4)

    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "best_source_cluster")

    # Query covers all
    intervals <- data.frame(chrom = "src", start = 0, end = 1000)

    result <- gintervals.liftover(intervals, chain)

    # Expectation:
    # The transitive cluster (mass 300) beats the single large distractor (mass 250).
    # If transitivity failed, A+B (200) vs C (100) vs D (250) -> D would win.

    expect_equal(nrow(result), 3, info = "Should keep transitive cluster A, B, C")

    result_ids <- sort(result$chain_id)
    expect_equal(result_ids, c(1, 2, 3), info = "Chain IDs should be A, B, C")

    # Ensure D is gone
    expect_false(4 %in% result$chain_id, info = "Chain D should be discarded")
})

test_that("best_source_cluster scenario 7: Mixed strand cluster (Overlap detected regardless of strand)", {
    local_db_state()

    setup_db(list(">chr1\n", paste(rep("A", 2000), collapse = ""), "\n"))
    chain_file <- new_chain_file()

    # Chain A (+): Source 100-200 -> Target 1000-1100
    write_chain_entry(chain_file, "src", 2000, "+", 100, 200, "chr1", 2000, "+", 1000, 1100, 1)

    # Chain B (-): Source 150-250 -> Target 1200-1300
    # Overlaps A in Source [150-200].
    # Note: Strand is negative, but start_src/end_src should be genomic coordinates for clustering.
    write_chain_entry(chain_file, "src", 2000, "+", 150, 250, "chr1", 2000, "-", 1200, 1300, 2)

    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "best_source_cluster")
    intervals <- data.frame(chrom = "src", start = 0, end = 500)
    result <- gintervals.liftover(intervals, chain)

    # Expectation: Overlap detected -> Cluster formed -> Both kept.
    expect_equal(nrow(result), 2, info = "Should cluster mixed-strand overlapping chains")
    expect_true(all(c(1, 2) %in% result$chain_id))
})

test_that("best_source_cluster scenario 8: Tie-breaking (First source position wins)", {
    local_db_state()

    setup_db(list(">chr1\n", paste(rep("A", 2000), collapse = ""), "\n"))
    chain_file <- new_chain_file()

    # Chain A: Source 100-200 (Mass 100)
    write_chain_entry(chain_file, "src", 2000, "+", 100, 200, "chr1", 2000, "+", 100, 200, 1)

    # Chain B: Source 500-600 (Mass 100) - Disjoint, Same Mass
    write_chain_entry(chain_file, "src", 2000, "+", 500, 600, "chr1", 2000, "+", 500, 600, 2)

    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "best_source_cluster")
    intervals <- data.frame(chrom = "src", start = 0, end = 1000)
    result <- gintervals.liftover(intervals, chain)

    # Expectation: Tie logic uses `>` (strictly greater).
    # Since A comes first (sorted by source start), it sets the 'best'.
    # B is equal, not greater, so A remains the winner.
    expect_equal(nrow(result), 1, info = "Should break ties deterministically")
    expect_equal(result$chain_id, 1, info = "Should keep the first chain (by source position) in a tie")
})
