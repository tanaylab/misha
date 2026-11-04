# Regression tests for the track import and creation findings of the 2026-11 review of 5.14.4.

restore_groot_on_exit()

# A small database of its own (chr1, chr2, chrX, chrM), in the given format, made the root for the
# calling test. Removed, and the previous root restored, when the test ends.
import_test_db <- function(format = "indexed", env = parent.frame()) {
    local_db_state(env)
    dir <- tempfile("import_review_")
    dir.create(dir)
    withr::defer(unlink(dir, recursive = TRUE), envir = env)
    sizes <- c(chr1 = 10000, chr2 = 10000, chrX = 5000, chrM = 500)
    records <- paste0(">", names(sizes), "\n", strrep("A", sizes), "\n")
    # an indexed database from one multi-FASTA, a per-chromosome one from a FASTA per chromosome
    fa <- if (format == "indexed") file.path(dir, "genome.fa") else file.path(dir, paste0(names(sizes), ".fa"))
    if (format == "indexed") {
        cat(records, sep = "", file = fa)
    } else {
        for (i in seq_along(sizes)) cat(records[i], file = fa[i])
    }
    db <- file.path(dir, "db")
    suppressMessages(gdb.create(groot = db, fasta = fa, format = format, verbose = FALSE))
    suppressMessages(gdb.init(db))
    db
}

write_file <- function(text, ext) {
    path <- tempfile(fileext = ext)
    cat(text, file = path)
    path
}

for (format in c("indexed", "per-chromosome")) {
    test_that(sprintf("a tab-delimited file without a trailing newline keeps its last line (%s)", format), {
        import_test_db(format)

        f <- write_file("chrom\tstart\tend\tvalue\nchr1\t0\t10\t1\nchr1\t20\t30\t2\nchrX\t5\t9\t3", ".tsv")
        gtrack.import("nonl", "", f, 0)
        res <- gextract("nonl", gintervals.all())
        expect_equal(as.character(res$chrom), c("chr1", "chr1", "chrX"))
        expect_equal(res$start, c(0, 20, 5))
        expect_equal(res$nonl, c(1, 2, 3))

        f <- write_file("chrom\tstart\tend\ta\tb\nchr1\t0\t10\t1\t2\nchr1\t20\t30\t3\t4", ".tsv")
        gtrack.array.import("arr_nonl", "", f)
        res <- gtrack.array.extract("arr_nonl", NULL, gintervals.all())
        expect_equal(res$start, c(0, 20))
        expect_equal(res$a, c(1, 3))
        expect_equal(res$b, c(2, 4))
    })
}
