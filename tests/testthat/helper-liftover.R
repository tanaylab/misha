write_chain_entry <- function(con, srcName, srcSize, srcStrand, srcStart, srcEnd,
                              tgtName, tgtSize, tgtStrand, tgtStart, tgtEnd, id) {
    # "chain 1000 source1 100 + 10 30 chr1 44 + 0 20 1"
    # where "source*" is the source (query) side and "chr*" is the target side.
    # Use %.0f for large integers to avoid %d overflow
    cat(
        sprintf(
            "chain 1000 %s %.0f %s %.0f %.0f %s %.0f %s %.0f %.0f %d\n",
            srcName, as.numeric(srcSize), srcStrand, as.numeric(srcStart), as.numeric(srcEnd),
            tgtName, as.numeric(tgtSize), tgtStrand, as.numeric(tgtStart), as.numeric(tgtEnd), id
        ),
        file = con, append = TRUE
    )
    # single block size line (tgtEnd - tgtStart), then blank line
    cat(sprintf("%.0f\n\n", as.numeric(tgtEnd - tgtStart)), file = con, append = TRUE)
}

setup_db <- function(chrom_defs) {
    # chrom_defs = list(chr=">chr1\nAAAA...", ... as strings to write)
    target_fasta <- tempfile(fileext = ".fasta")
    for (i in seq_along(chrom_defs)) {
        cat(chrom_defs[[i]], file = target_fasta, append = (i > 1))
    }
    target_db <- tempfile()
    gdb.create(groot = target_db, fasta = target_fasta, verbose = FALSE)
    gdb.init(target_db)
    withr::defer(
        {
            unlink(target_db, recursive = TRUE)
            unlink(target_fasta)
        },
        testthat::teardown_env()
    )
    target_fasta
}

new_chain_file <- function() {
    f <- tempfile(fileext = ".chain")
    withr::defer(unlink(f), testthat::teardown_env())
    f
}
