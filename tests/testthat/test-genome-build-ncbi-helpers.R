# Pure helpers added to support NCBI parity for gdb.install_intervals /
# gdb.build_genome:
#   - .ncbi_datasets_zip_url(accession, include)
#   - .ncbi_seqrep_to_alias_df(seqrep)
#   - .ncbi_ftp_assembly_dir(accession, assembly_name)
# All three are pure (no IO) so unit-tested directly. Wiring into
# .ncbi_fetch_assets / .build_seq_ncbi is exercised by the live-integration
# tests in test-genome-build-integration-zoonomia-smoke.R.

test_that(".ncbi_datasets_zip_url uses the requested include list", {
    u <- .ncbi_datasets_zip_url("GCF_000001405.40",
        include = c("GENOME_GFF", "SEQUENCE_REPORT")
    )
    expect_match(u, "/genome/accession/GCF_000001405.40/download\\?")
    expect_match(u, "include_annotation_type=GENOME_GFF,SEQUENCE_REPORT")
    expect_false(grepl("GENOME_FASTA", u))
})

test_that(".ncbi_datasets_zip_url accepts a single-element include", {
    u <- .ncbi_datasets_zip_url("GCF_X", include = "SEQUENCE_REPORT")
    expect_match(u, "include_annotation_type=SEQUENCE_REPORT$")
})

test_that(".ncbi_datasets_zip_url errors on empty include", {
    expect_error(.ncbi_datasets_zip_url("GCF_X", include = character(0)))
})

# Synthetic sequence_report covering the three roles + an unplaced row that's
# missing its GenBank twin (real assemblies do have those).
.mock_seqrep <- function() {
    data.frame(
        refseqAccession = c(
            "NC_000001.11", "NC_000023.11", "NC_012920.1",
            "NW_011332687.1", "NW_009646195.1"
        ),
        genbankAccession = c(
            "CM000663.2", "CM000685.2", "J01415.2",
            "KQ090024.1", ""
        ),
        chrName = c(
            "1", "X", "MT",
            "1", "Un"
        ),
        sequenceName = c(
            "1", "X", "MT",
            "HSCHR1_CTG3", "HSCHRUN_X1"
        ),
        role = c(
            "assembled-molecule", "assembled-molecule",
            "assembled-molecule",
            "unlocalized-scaffold", "unplaced-scaffold"
        ),
        length = c(248956422L, 156040895L, 16569L, 1000000L, 50000L),
        stringsAsFactors = FALSE
    )
}

test_that(".ncbi_seqrep_to_alias_df returns the expected five columns", {
    df <- .ncbi_seqrep_to_alias_df(.mock_seqrep())
    expect_named(df, c("accession", "genbank", "sequence_name", "chr_name", "ucsc"))
    expect_equal(nrow(df), 5L)
})

test_that(".ncbi_seqrep_to_alias_df preserves RefSeq + GenBank accessions verbatim", {
    df <- .ncbi_seqrep_to_alias_df(.mock_seqrep())
    expect_equal(df$accession[1], "NC_000001.11")
    expect_equal(df$genbank[1], "CM000663.2")
    # Empty GenBank cell stays empty - .detect_alias_column tolerates partial cols.
    expect_equal(df$genbank[5], "")
})

test_that(".ncbi_seqrep_to_alias_df chr_name == bare chrName from NCBI", {
    df <- .ncbi_seqrep_to_alias_df(.mock_seqrep())
    expect_equal(df$chr_name, c("1", "X", "MT", "1", "Un"))
})

test_that(".ncbi_seqrep_to_alias_df sequence_name uses chrName for assembled, sequenceName otherwise", {
    df <- .ncbi_seqrep_to_alias_df(.mock_seqrep())
    expect_equal(
        df$sequence_name,
        c("1", "X", "MT", "HSCHR1_CTG3", "HSCHRUN_X1")
    )
})

test_that(".ncbi_seqrep_to_alias_df ucsc column uses .ncbi_to_ucsc_name", {
    df <- .ncbi_seqrep_to_alias_df(.mock_seqrep())
    # Assembled: chr<N>; MT becomes chrM (UCSC convention, not chrMT).
    expect_equal(df$ucsc[1:3], c("chr1", "chrX", "chrM"))
    # Unlocalized: chr<chrName>_<encoded>_random.
    expect_equal(df$ucsc[4], "chr1_NW011332687v1_random")
    # Unplaced: chrUn_<encoded>.
    expect_equal(df$ucsc[5], "chrUn_NW009646195v1")
})

test_that(".ncbi_ftp_assembly_dir builds the GCF path for GRCh38.p14", {
    expect_equal(
        .ncbi_ftp_assembly_dir("GCF_000001405.40", "GRCh38.p14"),
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14"
    )
})

test_that(".ncbi_ftp_assembly_dir builds the GCA path", {
    expect_equal(
        .ncbi_ftp_assembly_dir("GCA_000001405.15", "GRCh38"),
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38"
    )
})

test_that(".ncbi_ftp_assembly_dir handles mouse GRCm39", {
    expect_equal(
        .ncbi_ftp_assembly_dir("GCF_000001635.27", "GRCm39"),
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39"
    )
})

test_that(".ncbi_ftp_assembly_dir handles non-zero leading digits", {
    # GCF_009914755.1 -> 009/914/755 (T2T-CHM13v2.0)
    expect_equal(
        .ncbi_ftp_assembly_dir("GCF_009914755.1", "T2T-CHM13v2.0"),
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0"
    )
})

test_that(".ncbi_parse_annotation_info exposes assembly_name when present", {
    rep <- list(reports = list(list(
        accession = "GCF_000001635.27",
        assembly_info = list(assembly_name = "GRCm39"),
        organism = list(organism_name = "Mus musculus", tax_id = 10090),
        annotation_info = list(name = "X", provider = "NCBI RefSeq")
    )))
    info <- .ncbi_parse_annotation_info(rep)
    expect_equal(info$assembly_name, "GRCm39")
})

test_that(".ncbi_parse_annotation_info exposes empty assembly_name when missing", {
    rep <- list(reports = list(list(
        accession = "GCA_X",
        organism = list(organism_name = "X"),
        annotation_info = list()
    )))
    info <- .ncbi_parse_annotation_info(rep)
    expect_equal(info$assembly_name, "")
})

test_that(".ncbi_ftp_assembly_name_from_dir extracts the right suffix", {
    # NCBI's parent FTP dir lists every assembly_name suffix for an accession's
    # version family. Helper picks the suffix matching a specific accession --
    # used as a fallback when /dataset_report is suppressed (older accessions).
    listing <- c(
        "<a href=\"GCF_000001635.20_GRCm38.p4/\">GCF_000001635.20_GRCm38.p4/</a>",
        "<a href=\"GCF_000001635.26_GRCm38.p6/\">GCF_000001635.26_GRCm38.p6/</a>",
        "<a href=\"GCF_000001635.27_GRCm39/\">GCF_000001635.27_GRCm39/</a>"
    )
    expect_equal(
        .ncbi_ftp_assembly_name_from_dir("GCF_000001635.26", listing),
        "GRCm38.p6"
    )
    expect_equal(
        .ncbi_ftp_assembly_name_from_dir("GCF_000001635.27", listing),
        "GRCm39"
    )
})

test_that(".ncbi_ftp_assembly_name_from_dir returns empty when not found", {
    listing <- c(
        "<a href=\"GCF_000001635.20_GRCm38.p4/\">GCF_000001635.20_GRCm38.p4/</a>"
    )
    expect_equal(
        .ncbi_ftp_assembly_name_from_dir("GCF_999999999.1", listing),
        ""
    )
})
