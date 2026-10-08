load_test_db()
test_that("indexed database preserves chrom_sizes.txt order for chromid mapping", {
    # Regression test for bug where chromosome identifiers resolved to wrong
    # chromosomes due to sorting mismatch between chrom_sizes.txt and ALLGENOME
    #
    # BUG SCENARIO: When chrom_sizes.txt is unsorted but genome.idx is sorted,
    # chromosome lookups would use the wrong chromid

    local_db_state()

    # Create a FASTA file with chromosomes in NON-alphabetical order
    test_fasta <- tempfile(fileext = ".fasta")
    cat(
        ">chr15\n", "AAAAAAAAAA", "\n",
        ">chr10\n", "CCCCCCCCCC", "\n",
        ">chr17_random\n", "GGGGGGGGGG", "\n",
        ">chr1\n", "TTTTTTTTTT", "\n",
        sep = "",
        file = test_fasta
    )

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    # Create database in indexed format
    withr::with_options(list(gmulticontig.indexed_format = TRUE), {
        suppressMessages(gdb.create(groot = test_db, fasta = test_fasta, verbose = FALSE))
    })

    # Initialize the database
    suppressMessages(gdb.init(test_db))

    # Check that chrom_sizes.txt is sorted (C++ import sorts it)
    chrom_sizes <- read.table(file.path(test_db, "chrom_sizes.txt"),
        header = FALSE, stringsAsFactors = FALSE
    )
    colnames(chrom_sizes) <- c("chrom", "size")

    # After C++ import, chromosomes should be sorted alphabetically
    expect_equal(chrom_sizes$chrom, sort(chrom_sizes$chrom))

    # Get ALLGENOME - should match chrom_sizes.txt order (NO SORTING in R for indexed)
    allgenome <- .misha$ALLGENOME[[1]]
    expect_equal(as.character(allgenome$chrom), chrom_sizes$chrom)

    # KEY TEST: Extract sequences from each chromosome
    # Before the fix, chromid mismatches would cause wrong sequences to be returned

    seq_chr1 <- gseq.extract(gintervals("chr1", 0, 10))
    seq_chr10 <- gseq.extract(gintervals("chr10", 0, 10))
    seq_chr15 <- gseq.extract(gintervals("chr15", 0, 10))
    seq_chr17 <- gseq.extract(gintervals("chr17_random", 0, 10))

    # Each chromosome should have its expected sequence
    expect_equal(seq_chr1, "TTTTTTTTTT")
    expect_equal(seq_chr10, "CCCCCCCCCC")
    expect_equal(seq_chr15, "AAAAAAAAAA")
    expect_equal(seq_chr17, "GGGGGGGGGG")

    # Verify ALLGENOME lengths are correct
    expect_equal(allgenome$end[allgenome$chrom == "chr1"], 10)
    expect_equal(allgenome$end[allgenome$chrom == "chr10"], 10)
    expect_equal(allgenome$end[allgenome$chrom == "chr15"], 10)
    expect_equal(allgenome$end[allgenome$chrom == "chr17_random"], 10)
})

test_that("a per-chromosome database gets the same chrom ids in a C and an en_US locale", {
    local_db_state()
    td <- withr::local_tempdir()
    db <- create_db_with_unsorted_chrom_sizes(file.path(td, "db"))
    cs <- utils::read.csv(file.path(db, "chrom_sizes.txt"),
        sep = "\t", header = FALSE, col.names = c("chrom", "size"), colClasses = c("character", "numeric")
    )
    chroms_by_id <- function(collate) {
        withr::with_collate(collate, {
            chrom_order <- misha:::.gdb.chrom_order(db, cs)
            gsetroot(db)
            list(order = chrom_order$names[chrom_order$id_order], allgenome = as.character(gintervals.all()$chrom))
        })
    }

    # the order of an en_US.UTF-8 session (ICU root collation), where a C locale's byte
    # order would put chr10 before chr1_KI270706v1_random
    expected <- c("chr1", "chr1_KI270706v1_random", "chr10", "chr2", "chrX")
    in_c <- chroms_by_id("C")
    expect_equal(in_c$order, expected)
    expect_equal(in_c$allgenome, expected)

    collate_available <- function(locale) {
        old <- Sys.getlocale("LC_COLLATE")
        withr::defer(Sys.setlocale("LC_COLLATE", old))
        nzchar(suppressWarnings(Sys.setlocale("LC_COLLATE", locale)))
    }
    en_us <- Filter(collate_available, c("en_US.UTF-8", "en_US.utf8"))
    skip_if(!length(en_us), "no en_US.UTF-8 locale")
    expect_equal(chroms_by_id(en_us[[1]]), in_c)
})

test_that("the chrom sort key orders names as order() does under ICU root collation", {
    skip_if_not(capabilities("ICU"), "R built without ICU")
    names <- c(
        "chr1", "chr10", "chr1_KI270706v1_random", "chr2", "chrX", "chrY", "chrM", "chrEBV",
        "chrUn_GL000220v1", "chrUn_KI270302v1", "chr6_GL000250v2_alt", "chr2_KI270894v1_alt",
        "HLA-A*01:01:01:01", "HLA-A*01:01:01:02N", "HLA-B*07:02:01", "HLA-DRB1*15:01:01:01",
        "NC_000001.11", "NC_000010.11", "NT_187361.1", "chr1.1", "chr1-1", "chr1_1", "chr1 1",
        "Chr1", "CHR1", "chr1a", "chr1A", "scaffold_10", "scaffold-10", "scaffold.10", "Scaffold_9",
        "super-scaffold_1", "Super_Scaffold_1", "MT", "mt", "contig#2", "contig@2", "contig(2)"
    )
    # and random printable-ASCII names
    set.seed(5)
    printable <- intToUtf8(32:126, multiple = TRUE)
    names <- unique(c(names, vapply(1:500, function(i) {
        paste(sample(printable, sample(1:12, 1), replace = TRUE), collapse = "")
    }, character(1))))

    old_collate <- icuGetCollate()
    withr::defer(icuSetCollate(locale = if (identical(old_collate, "ICU not in use")) "ASCII" else old_collate))
    icuSetCollate(locale = "root")
    expect_equal(names[order(misha:::.gdb.chrom_sort_key(names), method = "radix")], names[order(names)])
})
