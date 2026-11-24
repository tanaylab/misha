# gseq.shuffle: exact k-mer preserving shuffles (shared C++ core with pymisha).
# Re-roots into the example database for the real-sequence test; hand this
# process's overlay back to the next file in the parallel worker.
restore_groot_on_exit()

kmer_table <- function(s, k) {
    s <- toupper(s)
    L <- nchar(s)
    if (L < k) {
        return(table(character(0)))
    }
    table(substring(s, 1:(L - k + 1), k:L))
}

expect_same_kmers <- function(a, b, k) {
    expect_identical(kmer_table(a, k), kmer_table(b, k))
}

# every distinct sequence with the k-mer counts and the first k-1 bases of s
# (s holds only A/C/G/T), by depth-first search over the k-mers left to place
enumerate_shuffles <- function(s, k) {
    L <- nchar(s)
    if (k == 1) {
        chars <- strsplit(s, "")[[1]]
        perms <- function(v) {
            if (length(v) <= 1) {
                return(list(v))
            }
            out <- list()
            for (x in unique(v)) {
                rest <- v[-match(x, v)]
                out <- c(out, lapply(perms(rest), function(p) c(x, p)))
            }
            out
        }
        return(sort(vapply(perms(chars), paste, "", collapse = "")))
    }
    counts <- kmer_table(s, k)
    left <- as.integer(counts)
    names(left) <- names(counts)
    res <- character(0)
    rec <- function(cur, left) {
        if (nchar(cur) == L) {
            res <<- c(res, cur)
            return(invisible())
        }
        w <- substr(cur, nchar(cur) - k + 2, nchar(cur))
        for (km in names(left)) {
            if (left[[km]] > 0 && substr(km, 1, k - 1) == w) {
                left[[km]] <- left[[km]] - 1L
                rec(paste0(cur, substr(km, k, k)), left)
                left[[km]] <- left[[km]] + 1L
            }
        }
    }
    rec(substr(s, 1, k - 1), left)
    sort(unique(res))
}

# The shortcut of DeepLIFT / tangermeme: every (k-1)-mer keeps the edge it was
# last left through in s, and only its other out-edges are permuted. Valid, but
# not uniform. Written here from that description as the negative control.
shortcut_shuffle <- function(s, k, n) {
    L <- nchar(s)
    node <- substring(s, 1:(L - k + 2), (k - 1):L)
    m <- L - k + 1
    nexts <- split(seq_len(m) + 1L, node[seq_len(m)])
    vapply(seq_len(n), function(i) {
        ord <- lapply(nexts, function(v) {
            if (length(v) <= 2) {
                return(v)
            }
            c(sample(v[-length(v)]), v[length(v)])
        })
        ptr <- setNames(rep(1L, length(ord)), names(ord))
        out <- character(L)
        out[seq_len(k - 1)] <- strsplit(node[1], "")[[1]]
        cur <- node[1]
        for (p in k:L) {
            nx <- ord[[cur]][ptr[[cur]]]
            ptr[[cur]] <- ptr[[cur]] + 1L
            cur <- node[nx]
            out[p] <- substr(cur, k - 1, k - 1)
        }
        paste(out, collapse = "")
    }, "")
}

uniformity_p <- function(draws, valid) {
    obs <- table(factor(draws, levels = valid))
    stats::chisq.test(as.integer(obs), p = rep(1 / length(valid), length(valid)))$p.value
}

random_dna <- function(L, alphabet = c("A", "C", "G", "T")) {
    paste(sample(alphabet, L, replace = TRUE), collapse = "")
}

test_that("gseq.shuffle keeps the exact k-mer counts for k = 1..6 on random sequences", {
    set.seed(17)
    seqs <- c(
        vapply(c(2, 5, 10, 37, 100, 500, 2000), random_dna, ""),
        random_dna(300, c("A", "A", "A", "C", "G", "T")), # skewed
        random_dna(400, c("A", "T")), # two letters
        "ACGTACGTACGTACGT", # periodic
        paste(rep("CG", 50), collapse = "")
    )
    for (k in 1:6) {
        sh <- gseq.shuffle(seqs, k = k, n = 3, seed = k)
        expect_equal(dim(sh), c(length(seqs), 3L))
        for (i in seq_along(seqs)) {
            for (j in 1:3) {
                s <- sh[i, j]
                expect_equal(nchar(s), nchar(seqs[i]))
                expect_same_kmers(s, seqs[i], k)
                if (k > 1 && nchar(seqs[i]) >= k) {
                    expect_identical(substr(s, 1, k - 1), substr(seqs[i], 1, k - 1))
                    L <- nchar(s)
                    expect_identical(substr(s, L - k + 2, L), substr(seqs[i], L - k + 2, L))
                }
            }
        }
    }
    # it does move things around
    long <- seqs[6]
    expect_false(identical(gseq.shuffle(long, seed = 1), long))
})

test_that("gseq.shuffle keeps the k-mer counts of real genomic sequence", {
    local_db_state()
    gdb.init_examples()
    intervs <- gintervals(c(1, 1, 2), c(10000, 50000, 20000), c(12000, 55000, 21000))
    seqs <- gseq.extract(intervs)
    for (k in 1:6) {
        sh <- gseq.shuffle(seqs, k = k, seed = 42)
        for (i in seq_along(seqs)) {
            expect_same_kmers(sh[i], seqs[i], k)
        }
        expect_false(any(sh == seqs))
    }
    # an intervals data frame is read with gseq.extract
    expect_identical(gseq.shuffle(intervs, k = 3, seed = 5), gseq.shuffle(seqs, k = 3, seed = 5))
})

test_that("gseq.shuffle draws uniformly among all valid shuffles (k = 2)", {
    s <- "ACAGTCATGACTG"
    valid <- enumerate_shuffles(s, 2)
    expect_gt(length(valid), 20)
    draws <- as.vector(gseq.shuffle(s, k = 2, n = 30000, seed = 2024))
    expect_true(all(draws %in% valid))
    expect_setequal(unique(draws), valid)
    expect_gt(uniformity_p(draws, valid), 1e-3)
})

test_that("gseq.shuffle draws uniformly among all valid shuffles (k = 1, 3)", {
    s1 <- "AACCGT"
    valid1 <- enumerate_shuffles(s1, 1)
    expect_length(valid1, 180) # 6! / (2! 2!)
    draws1 <- as.vector(gseq.shuffle(s1, k = 1, n = 36000, seed = 3))
    expect_setequal(unique(draws1), valid1)
    expect_gt(uniformity_p(draws1, valid1), 1e-3)

    s3 <- "AACAACATAACATCATAA"
    valid3 <- enumerate_shuffles(s3, 3)
    expect_gt(length(valid3), 5)
    draws3 <- as.vector(gseq.shuffle(s3, k = 3, n = 20000, seed = 4))
    expect_true(all(draws3 %in% valid3))
    expect_setequal(unique(draws3), valid3)
    expect_gt(uniformity_p(draws3, valid3), 1e-3)
})

test_that("the uniformity test rejects the fixed-last-edge shortcut", {
    s <- "ACAGTCATGACTG"
    valid <- enumerate_shuffles(s, 2)
    withr::local_seed(99)
    draws <- shortcut_shuffle(s, 2, 5000)
    # what it draws is valid ...
    expect_true(all(draws %in% valid))
    # ... but some valid sequences can never come out, and the test sees it
    expect_lt(length(unique(draws)), length(valid))
    expect_lt(uniformity_p(draws, valid), 1e-10)
})

test_that("gseq.shuffle is reproducible and independent of R's RNG state given a seed", {
    seqs <- c("ACGTTGCAACGTAGCTAGCTTTACG", "GGGCCCATATCGCGATTA")
    a <- gseq.shuffle(seqs, seed = 11)
    expect_identical(gseq.shuffle(seqs, seed = 11), a)
    expect_false(identical(gseq.shuffle(seqs, seed = 12), a))

    # a seed neither reads nor advances .Random.seed
    set.seed(1)
    before <- .Random.seed
    gseq.shuffle(seqs, seed = 11)
    expect_identical(.Random.seed, before)

    # without a seed, set.seed() reproduces the call
    set.seed(5)
    b1 <- gseq.shuffle(seqs, n = 4)
    set.seed(5)
    b2 <- gseq.shuffle(seqs, n = 4)
    expect_identical(b1, b2)

    # each sequence has its own stream: earlier shuffles do not depend on later
    # sequences, and the first of n shuffles is the n = 1 shuffle
    expect_identical(gseq.shuffle(seqs[1], seed = 11), a[1])
    m <- gseq.shuffle(seqs, n = 5, seed = 11)
    expect_identical(m[, 1], a)
    # seeds are whole numbers up to 2^53, negative ones included
    expect_type(gseq.shuffle(seqs, seed = -3), "character")
    expect_type(gseq.shuffle(seqs, seed = 2^53), "character")
})

test_that("gseq.shuffle gives the same sequences as pymisha for the same seed", {
    # Golden values, asserted byte for byte by pymisha's tests/test_gseq_shuffle.py
    # too: the core and its generator are shared, so the outputs must agree.
    seqs <- c("ACGTTGCAACGTAGCTAGCTTTACG", "ccgNNacgtacgtTTGACGatcg", "AAAAAAAAAA", "ACGT")
    expect_identical(
        gseq.shuffle(seqs, k = 2, seed = 42),
        c("AACACTACGTTTAGTTAGCTGCGCG", "ccgNNacgattcgTGACGTtacg", "AAAAAAAAAA", "ACGT")
    )
    expect_identical(
        gseq.shuffle(seqs[1], k = 3, n = 3, seed = -7),
        matrix(c("ACGTAGCTTGCTAGCAACGTTTACG", "ACGTAGCAACGTTGCTTTAGCTACG", "ACGTAGCTTGCTACGTTTAGCAACG"), nrow = 1)
    )
    expect_identical(gseq.shuffle(seqs[1], k = 1, seed = 2^40), "ACTTTTCCAAGGACTTAGGAGCCTG")
})

test_that("gseq.shuffle keeps N's and other characters in place and shuffles between them", {
    set.seed(3)
    s <- paste0(
        random_dna(40), "NNNN", random_dna(30), "N", random_dna(3), "R",
        random_dna(50), "-", random_dna(1), "N"
    )
    sh <- gseq.shuffle(s, k = 2, n = 50, seed = 8)
    fixed <- gregexpr("[^ACGT]", s)[[1]]
    seg_bounds <- function(x) strsplit(x, "[^ACGT]+")[[1]]
    for (x in sh) {
        expect_identical(substring(x, fixed, fixed), substring(s, fixed, fixed))
        expect_same_kmers(x, s, 2)
        # each segment keeps its own counts
        a <- seg_bounds(x)
        b <- seg_bounds(s)
        for (i in seq_along(b)) expect_same_kmers(a[i], b[i], 2)
    }
    expect_gt(length(unique(sh)), 40)

    # N is never shuffled in, even with k = 1
    sh1 <- gseq.shuffle("ACGTNNACGT", k = 1, n = 20, seed = 1)
    expect_true(all(substr(sh1, 5, 6) == "NN"))
})

test_that("gseq.shuffle compares bases ignoring case and keeps each position's case", {
    s <- "ACGTacgtAAccGGttACGTAcgTacGT"
    sh <- gseq.shuffle(s, k = 2, n = 30, seed = 1)
    is_lower <- function(x) strsplit(x, "")[[1]] %in% letters
    for (x in sh) {
        expect_identical(is_lower(x), is_lower(s))
        expect_same_kmers(x, s, 2)
    }
    expect_identical(toupper(gseq.shuffle(s, k = 3, seed = 9)), gseq.shuffle(toupper(s), k = 3, seed = 9))
})

test_that("gseq.shuffle edge cases", {
    expect_identical(gseq.shuffle(character(0), seed = 1), character(0))
    expect_identical(dim(gseq.shuffle(character(0), n = 3, seed = 1)), c(0L, 3L))
    expect_identical(gseq.shuffle("", seed = 1), "")
    expect_identical(gseq.shuffle(NA_character_, seed = 1), NA_character_)
    expect_identical(gseq.shuffle(c("ACGT", NA), n = 2, seed = 1)[2, ], c(NA_character_, NA_character_))
    expect_identical(gseq.shuffle(NA, seed = 1), NA_character_)
    # at most k bases: one arrangement
    expect_identical(gseq.shuffle(c("A", "AC", "ACG"), k = 3, seed = 1), c("A", "AC", "ACG"))
    expect_identical(gseq.shuffle("ACGTACGTAC", k = 10, seed = 1), "ACGTACGTAC")
    expect_identical(gseq.shuffle("ACGTACGTAC", k = 32, seed = 1), "ACGTACGTAC")
    # a single letter, and sequences with no other arrangement
    expect_identical(gseq.shuffle("AAAAAAAAAAAA", k = 2, n = 3, seed = 1), matrix("AAAAAAAAAAAA", 1, 3))
    expect_identical(gseq.shuffle("ACGTTGCA", k = 4, seed = 1), "ACGTTGCA")
    expect_identical(gseq.shuffle("NNNNN", seed = 1), "NNNNN")
    # names are kept
    expect_named(gseq.shuffle(c(a = "ACGTAC", b = "GGTTAA"), seed = 1), c("a", "b"))
    expect_identical(rownames(gseq.shuffle(c(a = "ACGTAC", b = "GGTTAA"), n = 2, seed = 1)), c("a", "b"))
    # factors
    expect_identical(gseq.shuffle(factor("ACGTAC"), seed = 1), gseq.shuffle("ACGTAC", seed = 1))

    expect_error(gseq.shuffle("ACGT", k = 0), "k must be")
    expect_error(gseq.shuffle("ACGT", k = 33), "k must be")
    expect_error(gseq.shuffle("ACGT", k = 2.5), "k must be")
    expect_error(gseq.shuffle("ACGT", n = 0), "n must be")
    expect_error(gseq.shuffle("ACGT", seed = 1.5), "seed must be")
    expect_error(gseq.shuffle("ACGT", seed = 2^54), "seed must be")
    expect_error(gseq.shuffle("ACGT", seed = NA), "seed must be")
    expect_error(gseq.shuffle(1:3), "seqs must be")
})

test_that("gseq.shuffle handles long sequences quickly", {
    set.seed(2)
    s <- random_dna(2e5)
    t0 <- proc.time()[["elapsed"]]
    sh <- gseq.shuffle(s, k = 6, seed = 1)
    expect_lt(proc.time()[["elapsed"]] - t0, 30)
    expect_same_kmers(sh, s, 6)
})
