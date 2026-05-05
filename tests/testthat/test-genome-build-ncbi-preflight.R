# Pre-flight annotation availability check for NCBI source.
#
# Background: NCBI Datasets ships some assemblies (community submissions like
# Sanger TOL / DToL — about 80% of Phylo447) without any annotation. The
# Datasets zip then arrives with FASTA + sequence_report and no GFF, and our
# old code path warned and silently skipped 'genes' AFTER spending 800+ MB
# of download. The pre-flight calls /dataset_report (a few KB) first and
# trims 'genes' from the request before any heavy download.

# Mock dataset_report shapes covering the three real-world cases:
#   - RefSeq-curated GCF (annotation_info populated)
#   - GenBank-only community GCA (annotation_info empty/missing)
#   - GenBank with submitter annotation (annotation_info has provider != NCBI RefSeq)
.mock_report_annotated <- function() {
    list(reports = list(list(
        accession = "GCF_000001635.27",
        organism = list(organism_name = "Mus musculus", tax_id = 10090),
        annotation_info = list(
            name = "GCF_000001635.27-RS_2024_02",
            provider = "NCBI RefSeq",
            status = "Updated annotation",
            release_date = "2024-02-15"
        )
    )))
}

.mock_report_unannotated <- function() {
    list(reports = list(list(
        accession = "GCA_963575185.1",
        organism = list(organism_name = "Gorilla beringei", tax_id = 499232),
        annotation_info = list()
    )))
}

.mock_report_submitter_annotated <- function() {
    list(reports = list(list(
        accession = "GCA_047834645.1",
        organism = list(organism_name = "Nasalis larvatus", tax_id = 43780),
        annotation_info = list(
            name = "Annotation submitted by Northwest University",
            provider = "Northwest University"
        )
    )))
}

test_that(".ncbi_parse_annotation_info detects annotated assemblies", {
    info <- .ncbi_parse_annotation_info(.mock_report_annotated())
    expect_true(info$has_annotation)
    expect_equal(info$provider, "NCBI RefSeq")
    expect_equal(info$organism_tax_id, 10090)
})

test_that(".ncbi_parse_annotation_info detects unannotated assemblies", {
    info <- .ncbi_parse_annotation_info(.mock_report_unannotated())
    expect_false(info$has_annotation)
    expect_equal(info$provider, "")
    expect_equal(info$organism_tax_id, 499232)
})

test_that(".ncbi_parse_annotation_info accepts non-RefSeq submitter annotations", {
    # 'Has annotation' should be true whenever provider is non-empty, even if
    # it's a community submission rather than NCBI RefSeq. The user can still
    # use the GFF.
    info <- .ncbi_parse_annotation_info(.mock_report_submitter_annotated())
    expect_true(info$has_annotation)
    expect_equal(info$provider, "Northwest University")
})

test_that(".ncbi_resolve_sets_with_preflight passes through when annotated", {
    info <- .ncbi_parse_annotation_info(.mock_report_annotated())
    res <- .ncbi_resolve_sets_with_preflight(c("genes", "rmsk"), info, accession = "GCF_X")
    expect_equal(res$sets, c("genes", "rmsk"))
    expect_equal(length(res$warnings), 0L)
})

test_that(".ncbi_resolve_sets_with_preflight drops 'genes' when no annotation", {
    info <- .ncbi_parse_annotation_info(.mock_report_unannotated())
    res <- .ncbi_resolve_sets_with_preflight(c("genes", "rmsk"), info,
        accession = "GCA_963575185.1"
    )
    expect_false("genes" %in% res$sets)
    expect_true("rmsk" %in% res$sets)
    expect_equal(length(res$warnings), 1L)
    expect_match(res$warnings[[1L]], "no annotation")
    expect_match(res$warnings[[1L]], "GCA_963575185.1")
})

test_that(".ncbi_resolve_sets_with_preflight no-ops when 'genes' not requested", {
    info <- .ncbi_parse_annotation_info(.mock_report_unannotated())
    res <- .ncbi_resolve_sets_with_preflight(c("rmsk", "cgi"), info, accession = "GCA_X")
    expect_equal(res$sets, c("rmsk", "cgi"))
    expect_equal(length(res$warnings), 0L)
})

test_that(".ncbi_resolve_sets_with_preflight includes hint in warning when supplied", {
    info <- .ncbi_parse_annotation_info(.mock_report_unannotated())
    res <- .ncbi_resolve_sets_with_preflight(c("genes"), info,
        accession = "GCA_X",
        hint = "Try GCF_111.1 (NCBI RefSeq Annotation Release 100) for the same taxon."
    )
    expect_match(res$warnings[[1L]], "GCF_111.1")
    expect_match(res$warnings[[1L]], "annotated alternative", ignore.case = TRUE)
})

# ---- Live integration tests (gated; require network + NCBI Datasets API) ----

test_that(".ncbi_dataset_report + parse round-trips against live NCBI (annotated)", {
    skip_if_offline("api.ncbi.nlm.nih.gov")
    skip_on_cran()
    rep <- tryCatch(.ncbi_dataset_report("GCF_000001635.27", timeout = 30),
        error = function(e) skip(paste("NCBI Datasets fetch failed:", conditionMessage(e)))
    )
    info <- .ncbi_parse_annotation_info(rep)
    expect_true(info$has_annotation)
    expect_equal(info$organism_name, "Mus musculus")
    expect_equal(info$organism_tax_id, 10090L)
})

test_that(".ncbi_dataset_report + parse round-trips against live NCBI (unannotated)", {
    skip_if_offline("api.ncbi.nlm.nih.gov")
    skip_on_cran()
    rep <- tryCatch(.ncbi_dataset_report("GCA_963575185.1", timeout = 30),
        error = function(e) skip(paste("NCBI Datasets fetch failed:", conditionMessage(e)))
    )
    info <- .ncbi_parse_annotation_info(rep)
    expect_false(info$has_annotation)
    expect_equal(info$organism_name, "Gorilla beringei")
})

test_that(".ncbi_suggest_annotated_alternative finds a RefSeq companion when one exists", {
    skip_if_offline("api.ncbi.nlm.nih.gov")
    skip_on_cran()
    # Sus scrofa (tax 9823): GCA_000003025.7 is unannotated; GCF_000003025.6
    # carries NCBI RefSeq Annotation Release 106. The suggestion endpoint
    # should surface a GCF_* accession.
    hint <- .ncbi_suggest_annotated_alternative(9823L, "GCA_000003025.7")
    expect_match(hint, "GCF_")
})
