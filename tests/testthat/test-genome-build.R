# Unit tests for gdb.build_genome() and friends. All offline — no network,
# no UCSC binary required. End-to-end build tests are deferred (would need
# skip_on_cran + skip_if_offline + a real upstream).

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

write_temp_yaml <- function(text) {
    f <- tempfile(fileext = ".yaml")
    writeLines(text, f)
    f
}

# ---------------------------------------------------------------------------
# Registry parsing / normalization / validation
# ---------------------------------------------------------------------------

test_that(".normalize_recipe accepts string (legacy local form)", {
    r <- misha:::.normalize_recipe("/tmp/groot", "hg19")
    expect_equal(r$source, "local")
    expect_equal(r$path, "/tmp/groot")
})

test_that(".normalize_recipe preserves a dict and requires 'source'", {
    r <- misha:::.normalize_recipe(list(source = "ucsc", assembly = "hg38"), "hg38")
    expect_equal(r$source, "ucsc")
    expect_equal(r$assembly, "hg38")

    expect_error(
        misha:::.normalize_recipe(list(assembly = "hg38"), "hg38"),
        "missing required 'source:'"
    )
    expect_error(
        misha:::.normalize_recipe(list(source = "made_up"), "x"),
        "unknown source"
    )
    expect_error(
        misha:::.normalize_recipe(42, "x"),
        "invalid registry entry"
    )
})

test_that(".validate_recipe enforces per-source required fields", {
    expect_error(misha:::.validate_recipe(list(source = "ucsc"), "x"), "assembly")
    expect_error(misha:::.validate_recipe(list(source = "ncbi"), "x"), "accession")
    expect_error(
        misha:::.validate_recipe(list(source = "ncbi", accession = "garbage"), "x"),
        "does not match"
    )
    expect_error(misha:::.validate_recipe(list(source = "local"), "x"), "path")
    expect_error(misha:::.validate_recipe(list(source = "manual"), "x"), "fasta")

    # OK paths
    expect_silent(misha:::.validate_recipe(list(source = "ucsc", assembly = "hg38"), "x"))
    expect_silent(misha:::.validate_recipe(list(source = "ncbi", accession = "GCF_000001405.40"), "x"))
    expect_silent(misha:::.validate_recipe(list(source = "manual", fasta = "url"), "x"))

    # NCBI chrom_naming
    expect_error(
        misha:::.validate_recipe(list(
            source = "ncbi", accession = "GCF_000001405.40",
            chrom_naming = "bad"
        ), "x"),
        "chrom_naming"
    )
    expect_silent(misha:::.validate_recipe(
        list(source = "ncbi", accession = "GCF_000001405.40", chrom_naming = "ucsc"), "x"
    ))
    expect_silent(misha:::.validate_recipe(
        list(source = "ncbi", accession = "GCF_000001405.40", chrom_naming = "sequence_name"), "x"
    ))
})

test_that(".parse_genome_registry handles valid + invalid YAML", {
    f <- write_temp_yaml(c(
        "genome:",
        "  hg38:",
        "    source: ucsc",
        "    assembly: hg38",
        "  custom: /tmp/local"
    ))
    on.exit(unlink(f))
    entries <- misha:::.parse_genome_registry(f)
    expect_named(entries, c("hg38", "custom"), ignore.order = TRUE)
    expect_equal(entries$custom, "/tmp/local")
    expect_equal(entries$hg38$source, "ucsc")

    # No genome: key -> empty list, not error
    f2 <- write_temp_yaml(c("other: stuff"))
    on.exit(unlink(f2), add = TRUE)
    expect_equal(length(misha:::.parse_genome_registry(f2)), 0)

    # Missing file
    expect_error(
        misha:::.parse_genome_registry("/nonexistent/registry.yaml"),
        "does not exist"
    )
})

# ---------------------------------------------------------------------------
# Resolution chain
# ---------------------------------------------------------------------------

test_that(".resolve_genome explicit registry overrides built-in", {
    f <- write_temp_yaml(c(
        "genome:",
        "  hg38:",
        "    source: manual",
        "    fasta: file:///tmp/custom.fa"
    ))
    on.exit(unlink(f))
    res <- misha:::.resolve_genome("hg38", registry = f)
    expect_equal(res$recipe$source, "manual")
    expect_match(res$resolved_from, "registry arg")
})

test_that(".resolve_genome falls back to built-in registry", {
    res <- misha:::.resolve_genome("mm10")
    expect_equal(res$recipe$source, "ucsc")
    expect_equal(res$recipe$assembly, "mm10")
    expect_match(res$resolved_from, "built-in")
})

test_that(".resolve_genome pattern fallback for GC[FA]_* accessions routes to ucsc-hub", {
    res <- misha:::.resolve_genome("GCF_009806435.1")
    expect_equal(res$recipe$source, "ucsc-hub")
    expect_equal(res$recipe$accession, "GCF_009806435.1")
    expect_match(res$resolved_from, "pattern fallback")

    res2 <- misha:::.resolve_genome("GCA_000001405.29")
    expect_equal(res2$recipe$source, "ucsc-hub")
})

test_that(".resolve_genome errors clearly for unknown name", {
    expect_error(misha:::.resolve_genome("totally_made_up_genome"), "not found")
})

# ---------------------------------------------------------------------------
# UCSC URL builder
# ---------------------------------------------------------------------------

test_that(".ucsc_urls constructs the right URLs for hg38", {
    urls <- misha:::.ucsc_urls("hg38",
        annotations = c("genes", "rmsk", "cpgIsland", "cytoband")
    )
    expect_match(urls$fasta, "/goldenPath/hg38/bigZips/hg38\\.fa\\.gz$")
    expect_match(urls$genes, "/database/ncbiRefSeq\\.txt\\.gz$")
    expect_match(urls$annots, "/database/ncbiRefSeqLink\\.txt\\.gz$")
    expect_match(urls$rmsk, "/database/rmsk\\.txt\\.gz$")
    expect_match(urls$cpgIsland, "/database/cpgIslandExt\\.txt\\.gz$")
    expect_match(urls$cytoband, "/database/cytoBandIdeo\\.txt\\.gz$")
})

test_that(".ucsc_urls omits annotations not requested", {
    urls <- misha:::.ucsc_urls("mm10", annotations = c("genes"))
    expect_true(!is.null(urls$genes))
    expect_null(urls$rmsk)
    expect_null(urls$cpgIsland)
    expect_null(urls$cytoband)
})

# ---------------------------------------------------------------------------
# Architecture detection
# ---------------------------------------------------------------------------

test_that(".detect_arch returns one of the known platforms on this host", {
    arch <- misha:::.detect_arch()
    expect_true(arch %in% c("linux_x86_64", "macos_x86_64", "macos_arm64"))
})

# ---------------------------------------------------------------------------
# Binary path resolver
# ---------------------------------------------------------------------------

test_that(".gff3_to_genepred_path honors MISHA_GFF3_TO_GENEPRED", {
    skip_if_not_installed("withr")
    stub <- tempfile(fileext = ".sh")
    writeLines("#!/bin/sh\necho 'stub'", stub)
    Sys.chmod(stub, "0755")
    on.exit(unlink(stub))

    withr::with_envvar(
        c(MISHA_GFF3_TO_GENEPRED = stub),
        {
            p <- misha:::.gff3_to_genepred_path()
            expect_equal(normalizePath(p), normalizePath(stub))
        }
    )
})

test_that(".gff3_to_genepred_path errors if env var points at nonexistent file", {
    skip_if_not_installed("withr")
    withr::with_envvar(
        c(MISHA_GFF3_TO_GENEPRED = "/nonexistent/binary"),
        {
            expect_error(misha:::.gff3_to_genepred_path(), "non-existent")
        }
    )
})

test_that(".gff3_to_genepred_path returns NULL when nothing is cached and no env var", {
    skip_if_not_installed("withr")
    # Guard: make sure no env var is set.
    withr::with_envvar(
        c(MISHA_GFF3_TO_GENEPRED = ""),
        {
            # Move any existing cache out of the way for the duration.
            cached <- misha:::.gff3_to_genepred_cache_path()
            if (file.exists(cached)) {
                bak <- paste0(cached, ".bak.", Sys.getpid())
                file.rename(cached, bak)
                on.exit(file.rename(bak, cached), add = TRUE)
            }
            expect_null(misha:::.gff3_to_genepred_path())
        }
    )
})

# ---------------------------------------------------------------------------
# UCSC table parsers
# ---------------------------------------------------------------------------

test_that(".parse_ucsc_rmsk handles a minimal record", {
    f <- tempfile(fileext = ".txt")
    on.exit(unlink(f))
    # 17 cols, two rows
    writeLines(c(
        paste0(
            "0\t1000\t100\t10\t10\t", # bin/swScore/milliDiv/milliDel/milliIns
            "chr1\t100\t200\t-1000\t+\t", # genoName/genoStart/genoEnd/genoLeft/strand
            "L1HS\tLINE\tL1\t", # repName/repClass/repFamily
            "1\t100\t-50\t1"
        ), # repStart/repEnd/repLeft/id
        paste0(
            "0\t800\t100\t10\t10\t",
            "chr2\t500\t800\t-500\t-\t",
            "AluSx\tSINE\tAlu\t",
            "1\t300\t-100\t2"
        )
    ), f)
    df <- misha:::.parse_ucsc_rmsk(f)
    expect_equal(nrow(df), 2)
    expect_equal(df$chrom, c("chr1", "chr2"))
    expect_equal(df$strand, c(1L, -1L))
    expect_equal(df$class, c("LINE", "SINE"))
    expect_equal(df$family, c("L1", "Alu"))
})

test_that(".parse_ucsc_cpg_island handles a minimal record", {
    f <- tempfile(fileext = ".txt")
    on.exit(unlink(f))
    # 11 cols
    writeLines("0\tchr1\t100\t300\tCpG: 25\t200\t25\t150\t25.0\t75.0\t0.85", f)
    df <- misha:::.parse_ucsc_cpg_island(f)
    expect_equal(nrow(df), 1)
    expect_equal(df$chrom, "chr1")
    expect_equal(df$start, 100)
    expect_equal(df$end, 300)
})

test_that(".parse_ucsc_cytoband handles a minimal record", {
    f <- tempfile(fileext = ".txt")
    on.exit(unlink(f))
    # 5 cols
    writeLines(c(
        "chr1\t0\t2300000\tp36.33\tgneg",
        "chr1\t2300000\t5400000\tp36.32\tgpos25"
    ), f)
    df <- misha:::.parse_ucsc_cytoband(f)
    expect_equal(nrow(df), 2)
    expect_equal(df$stain, c("gneg", "gpos25"))
})

# ---------------------------------------------------------------------------
# gdb.list_genomes / gdb.genome_info
# ---------------------------------------------------------------------------

test_that("gdb.list_genomes includes built-in entries", {
    df <- gdb.list_genomes()
    expect_true("hg38" %in% df$name)
    expect_true("mm10" %in% df$name)
    expect_equal(unique(df$source[df$name %in% c("hg38", "mm10")]), "ucsc")
})

test_that("gdb.list_genomes accepts an explicit registry", {
    f <- write_temp_yaml(c(
        "genome:",
        "  fake: {source: manual, fasta: file:///tmp/x.fa}"
    ))
    on.exit(unlink(f))
    df <- gdb.list_genomes(registry = f)
    expect_true("fake" %in% df$name)
    expect_equal(df$source[df$name == "fake"], "manual")
})

test_that("gdb.genome_info returns a recipe", {
    info <- gdb.genome_info("hg38")
    expect_equal(info$recipe$source, "ucsc")
    expect_equal(info$recipe$assembly, "hg38")
})

test_that("gdb.genome_info on GC[FA]_* accession uses pattern fallback to ucsc-hub", {
    info <- gdb.genome_info("GCF_009806435.1")
    expect_equal(info$recipe$source, "ucsc-hub")
    expect_equal(info$recipe$accession, "GCF_009806435.1")
})

# ---------------------------------------------------------------------------
# UCSC TSV preprocessing helpers
# ---------------------------------------------------------------------------

test_that(".heal_ucsc_tsv_escapes collapses backslash-newline / -tab / stray CR", {
    in_path <- tempfile(fileext = ".txt")
    out_path <- tempfile(fileext = ".txt")
    on.exit(unlink(c(in_path, out_path)))
    # Three pathological lines:
    #   1. clean line
    #   2. line with backslash-newline split (Windows-style: \\\r\n)
    #   3. line with backslash-tab split
    f <- file(in_path, "wb")
    writeBin(charToRaw("clean\trow\tone\n"), f)
    # row 2: cols [a, b, c\<CR><LF>more, d]
    writeBin(charToRaw("a\tb\tc\\\r\nmore\td\n"), f)
    # row 3: cols [x, y, z\<TAB>extra, w]
    writeBin(charToRaw("x\ty\tz\\\textra\tw\n"), f)
    close(f)

    misha:::.heal_ucsc_tsv_escapes(in_path, out_path)
    healed <- readLines(out_path, warn = FALSE)
    expect_equal(length(healed), 3)
    expect_equal(healed[1], "clean\trow\tone")
    expect_match(healed[2], "^a\tb\tc more\td$")
    expect_match(healed[3], "^x\ty\tz extra\tw$")
})

test_that(".normalize_ucsc_tsv pads short rows and joins overflow", {
    in_path <- tempfile(fileext = ".txt")
    out_path <- tempfile(fileext = ".txt")
    on.exit(unlink(c(in_path, out_path)))
    writeLines(c(
        "a\tb\tc", # short: NF=3
        "1\t2\t3\t4", # exact: NF=4
        "x\ty\tz\tw\textra" # overflow: NF=5
    ), in_path)
    misha:::.normalize_ucsc_tsv(in_path, out_path, n_cols = 4)
    out <- readLines(out_path, warn = FALSE)
    expect_equal(length(out), 3)
    # Each line must have exactly 3 tabs (4 fields). Use char count instead
    # of strsplit (which drops trailing empties).
    n_tabs <- vapply(out, function(s) sum(charToRaw(s) == as.raw(0x09)), integer(1))
    expect_equal(unname(n_tabs), c(3L, 3L, 3L))
    # short row padded with empty tail
    expect_equal(out[1], "a\tb\tc\t")
    # exact row preserved
    expect_equal(out[2], "1\t2\t3\t4")
    # overflow joined into last col
    expect_equal(out[3], "x\ty\tz\tw extra")
})

test_that(".normalize_ucsc_genepred handles 15-col (gff3ToGenePred) and 16-col (UCSC table)", {
    in15 <- tempfile(fileext = ".txt")
    in16 <- tempfile(fileext = ".txt")
    out15 <- tempfile(fileext = ".txt")
    out16 <- tempfile(fileext = ".txt")
    on.exit(unlink(c(in15, in16, out15, out16)))

    # 15-col: gff3ToGenePred output (no bin column)
    writeLines(paste(c(
        "NM_001", "chr1", "+", "100", "200", "120", "180", "1",
        "100,", "200,", "0", "GENE1", "cmpl", "cmpl", "0,"
    ), collapse = "\t"), in15)
    misha:::.normalize_ucsc_genepred(in15, out15)
    fields <- strsplit(readLines(out15, warn = FALSE), "\t", fixed = TRUE)[[1]]
    expect_equal(length(fields), 12)
    expect_equal(fields[1], "NM_001") # ID = name
    expect_equal(fields[2], "chr1") # CHROM
    expect_equal(fields[3], "+") # STRAND

    # 16-col: UCSC table (with bin column)
    writeLines(paste(c(
        "585", "NM_001", "chr1", "+", "100", "200", "120", "180",
        "1", "100,", "200,", "0", "GENE1", "cmpl", "cmpl", "0,"
    ), collapse = "\t"), in16)
    misha:::.normalize_ucsc_genepred(in16, out16)
    fields16 <- strsplit(readLines(out16, warn = FALSE), "\t", fixed = TRUE)[[1]]
    expect_equal(length(fields16), 12)
    expect_equal(fields16[1], "NM_001") # ID = name (bin dropped)
    expect_equal(fields16[2], "chr1")
    expect_equal(fields16[3], "+")
})

test_that(".normalize_ucsc_genepred rejects unsupported column count", {
    in_path <- tempfile(fileext = ".txt")
    out_path <- tempfile(fileext = ".txt")
    on.exit(unlink(c(in_path, out_path)))
    writeLines("a\tb\tc\td", in_path) # 4 cols, unsupported
    expect_error(
        misha:::.normalize_ucsc_genepred(in_path, out_path),
        "unsupported column count 4"
    )
})

# ---------------------------------------------------------------------------
# NCBI sequence_report parsing + rename-map building + stream rewrites
# ---------------------------------------------------------------------------

write_temp_seqrep <- function() {
    f <- tempfile(fileext = ".jsonl")
    writeLines(c(
        '{"refseqAccession":"NC_067374.1","genbankAccession":"CM028932.1","chrName":"1","sequenceName":"contig_2989","role":"assembled-molecule","length":208594813}',
        '{"refseqAccession":"NC_001913.1","genbankAccession":"CM_NA","chrName":"MT","sequenceName":"mito","role":"assembled-molecule","length":17245}',
        '{"refseqAccession":"NW_026256937.1","genbankAccession":"VIYN02000001.1","chrName":"Un","sequenceName":"contig_1","role":"unplaced-scaffold","length":1333029}'
    ), f)
    f
}

test_that(".parse_ncbi_sequence_report extracts the right fields", {
    f <- write_temp_seqrep()
    on.exit(unlink(f))
    df <- misha:::.parse_ncbi_sequence_report(f)
    expect_equal(nrow(df), 3)
    expect_equal(df$refseqAccession, c("NC_067374.1", "NC_001913.1", "NW_026256937.1"))
    expect_equal(df$chrName, c("1", "MT", "Un"))
    expect_equal(df$role, c("assembled-molecule", "assembled-molecule", "unplaced-scaffold"))
    expect_equal(df$length[1], 208594813)
})

test_that(".build_ncbi_rename_map dispatches per chrom_naming", {
    f <- write_temp_seqrep()
    on.exit(unlink(f))
    seqrep <- misha:::.parse_ncbi_sequence_report(f)

    seq_map <- misha:::.build_ncbi_rename_map(seqrep, "sequence_name")
    expect_equal(unname(seq_map[["NC_067374.1"]]), "1")
    expect_equal(unname(seq_map[["NC_001913.1"]]), "MT")
    # Unplaced falls back to refseq accession (chrName=Un would collide).
    expect_equal(unname(seq_map[["NW_026256937.1"]]), "NW_026256937.1")

    ucsc_map <- misha:::.build_ncbi_rename_map(seqrep, "ucsc")
    expect_equal(unname(ucsc_map[["NC_067374.1"]]), "chr1")
    expect_equal(unname(ucsc_map[["NC_001913.1"]]), "chrM") # not chrMT
    expect_equal(unname(ucsc_map[["NW_026256937.1"]]), "chrUn_NW026256937v1")

    acc_map <- misha:::.build_ncbi_rename_map(seqrep, "accession")
    expect_equal(unname(acc_map[["NC_067374.1"]]), "NC_067374.1")
})

test_that(".build_ncbi_rename_map errors on duplicate target names", {
    # Two assembled molecules with the same chrName would collide.
    f <- tempfile(fileext = ".jsonl")
    on.exit(unlink(f))
    writeLines(c(
        '{"refseqAccession":"NC_001.1","chrName":"1","role":"assembled-molecule","length":1}',
        '{"refseqAccession":"NC_002.1","chrName":"1","role":"assembled-molecule","length":1}'
    ), f)
    seqrep <- misha:::.parse_ncbi_sequence_report(f)
    expect_error(
        misha:::.build_ncbi_rename_map(seqrep, "sequence_name"),
        "duplicate names"
    )
})

test_that(".rename_fasta_headers rewrites ID and preserves body", {
    in_path <- tempfile(fileext = ".fa")
    out_path <- tempfile(fileext = ".fa")
    on.exit(unlink(c(in_path, out_path)))
    writeLines(c(
        ">NC_067374.1 chr 1 some description",
        "ACGTACGT",
        "GGGGCCCC",
        ">NC_001913.1",
        "ATGCATGC",
        ">NW_unmapped.1 unmapped contig",
        "TTTT"
    ), in_path)
    rmap <- c("NC_067374.1" = "1", "NC_001913.1" = "MT")
    misha:::.rename_fasta_headers(in_path, out_path, rmap, verbose = FALSE)
    out <- readLines(out_path, warn = FALSE)
    expect_equal(out[1], ">1") # description dropped
    expect_equal(out[4], ">MT")
    expect_equal(out[6], ">NW_unmapped.1") # unmapped: ID kept, description dropped
    expect_equal(out[2], "ACGTACGT") # body preserved
    expect_equal(out[3], "GGGGCCCC")
})

test_that(".rename_gff3_seqids rewrites col 1 and preserves comments", {
    in_path <- tempfile(fileext = ".gff")
    out_path <- tempfile(fileext = ".gff")
    on.exit(unlink(c(in_path, out_path)))
    writeLines(c(
        "##gff-version 3",
        "#!processor RefSeq",
        "NC_067374.1\tRefSeq\tgene\t100\t200\t.\t+\t.\tID=g1",
        "NC_067374.1\tRefSeq\texon\t100\t200\t.\t+\t.\tID=e1;Parent=g1",
        "NC_unmapped.1\tRefSeq\tgene\t1\t10\t.\t+\t.\tID=g2"
    ), in_path)
    rmap <- c("NC_067374.1" = "1")
    misha:::.rename_gff3_seqids(in_path, out_path, rmap, verbose = FALSE)
    out <- readLines(out_path, warn = FALSE)
    expect_equal(out[1], "##gff-version 3") # comments preserved
    expect_equal(out[2], "#!processor RefSeq")
    expect_match(out[3], "^1\t") # renamed
    expect_match(out[4], "^1\t")
    expect_match(out[5], "^NC_unmapped.1\t") # unmapped passed through
})

# ---------------------------------------------------------------------------
# .compute_chrom_aliases TSV augmentation
# ---------------------------------------------------------------------------

test_that(".compute_chrom_aliases without groot still does chr-prefix + MT toggles", {
    a <- misha:::.compute_chrom_aliases(c("1", "X", "MT"))
    # Each canonical maps to itself
    expect_equal(unname(a[["1"]]), "1")
    expect_equal(unname(a[["X"]]), "X")
    expect_equal(unname(a[["MT"]]), "MT")
    # chr-prefix aliases
    expect_equal(unname(a[["chr1"]]), "1")
    expect_equal(unname(a[["chrX"]]), "X")
    # mitochondrial aliases
    expect_equal(unname(a[["M"]]), "MT")
    expect_equal(unname(a[["chrM"]]), "MT")
})

test_that(".compute_chrom_aliases adds refseq/genbank/sequenceName/chrName when chrom_aliases.tsv is present", {
    g <- tempfile()
    dir.create(g)
    on.exit(unlink(g, recursive = TRUE))
    writeLines(c(
        "canonical\trefseqAccession\tgenbankAccession\tsequenceName\tchrName\trole\tlength",
        "1\tNC_067374.1\tCM028932.1\tcontig_2989\t1\tassembled-molecule\t208594813",
        "MT\tNC_001913.1\tAJ001588.1\tMT\tMT\tassembled-molecule\t17245"
    ), file.path(g, "chrom_aliases.tsv"))

    a <- misha:::.compute_chrom_aliases(c("1", "MT"), groot = g)

    # Existing aliases still resolve.
    expect_equal(unname(a[["1"]]), "1")
    expect_equal(unname(a[["chr1"]]), "1")
    expect_equal(unname(a[["chrM"]]), "MT")
    # New TSV-driven aliases.
    expect_equal(unname(a[["NC_067374.1"]]), "1")
    expect_equal(unname(a[["CM028932.1"]]), "1")
    expect_equal(unname(a[["contig_2989"]]), "1")
    expect_equal(unname(a[["NC_001913.1"]]), "MT")
})

test_that(".compute_chrom_aliases ignores TSV rows whose canonical is not among `chroms`", {
    # Defensive: a stale TSV referring to renamed chroms must not be applied.
    g <- tempfile()
    dir.create(g)
    on.exit(unlink(g, recursive = TRUE))
    writeLines(c(
        "canonical\trefseqAccession\tgenbankAccession\tsequenceName\tchrName\trole\tlength",
        "ANCIENT_NAME\tNC_999999.1\t.\t.\t.\tassembled\t1"
    ), file.path(g, "chrom_aliases.tsv"))
    a <- misha:::.compute_chrom_aliases(c("1"), groot = g)
    expect_false("NC_999999.1" %in% names(a))
})

test_that(".compute_chrom_aliases skips alias names that already exist (e.g. canonical itself)", {
    # If the TSV's refseqAccession happens to equal the canonical (e.g. when
    # chrom_naming = 'accession'), don't enter it as a self-mapping.
    g <- tempfile()
    dir.create(g)
    on.exit(unlink(g, recursive = TRUE))
    writeLines(c(
        "canonical\trefseqAccession\tgenbankAccession\tsequenceName\tchrName\trole\tlength",
        "NC_067374.1\tNC_067374.1\tCM028932.1\tcontig_2989\t1\tassembled-molecule\t208594813"
    ), file.path(g, "chrom_aliases.tsv"))
    a <- misha:::.compute_chrom_aliases(c("NC_067374.1"), groot = g)
    expect_equal(unname(a[["NC_067374.1"]]), "NC_067374.1")
    expect_equal(unname(a[["CM028932.1"]]), "NC_067374.1")
    expect_equal(unname(a[["contig_2989"]]), "NC_067374.1")
    expect_equal(unname(a[["1"]]), "NC_067374.1")
})

# ---------------------------------------------------------------------------
# Stub-based dispatch sanity check
# ---------------------------------------------------------------------------

test_that("gdb.build_genome refuses to overwrite existing path", {
    d <- tempfile()
    dir.create(d)
    on.exit(unlink(d, recursive = TRUE))
    expect_error(gdb.build_genome("hg38", path = d), "already exists")
})
