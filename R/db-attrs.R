# Attribute functions
#' Sets read-only track attributes
#'
#' Sets read-only track attributes.
#'
#' This function sets the list of read-only track attributes. The specified
#' attributes may or may not already exist in the tracks.
#'
#' If 'attrs' is 'NULL' the list of read-only attributes is emptied.
#'
#' @param attrs a vector of read-only attributes names or 'NULL'
#' @return None.
#' @seealso \code{\link{gdb.get_readonly_attrs}},
#' \code{\link{gtrack.attr.get}}, \code{\link{gtrack.attr.set}}
#' @keywords ~attr ~attribute
#' @export gdb.set_readonly_attrs
gdb.set_readonly_attrs <- function(attrs) {
    .gcheckroot()

    filename <- paste(get("GROOT", envir = .misha), ".ro_attributes", sep = "/")

    if (is.null(attrs)) {
        unlink(filename)
    } else {
        attrs <- as.character(attrs)

        idx <- which(duplicated(attrs))[1]
        if (!is.na(idx)) {
            stop(sprintf("Attribute %s appears more than once", attrs[idx]), call. = FALSE)
        }

        idx <- which(attrs == "")[1]
        if (!is.na(idx)) {
            stop("Attribute name cannot be an empty string", call. = FALSE)
        }

        f <- file(filename, "wb")
        serialize(attrs, f)
        close(f)
    }
    retv <- 0 # suppress return value
}

#' Creates a new Genomic Database
#'
#' Creates a new Genomic Database.
#'
#' This function creates a new Genomic Database at the location specified by
#' 'groot'. FASTA files are converted to 'Seq' format and appropriate
#' 'chrom_sizes.txt' file is generated (see "User Manual" for more details).
#'
#' Two database formats are supported:
#' \itemize{
#'   \item \strong{indexed}: Single genome.seq + genome.idx (default). Recommended for
#'         genomes with many contigs. Provides better performance and scalability.
#'   \item \strong{per-chromosome}: Separate .seq file per contig.
#' }
#'
#' If 'genes.file' is not 'NULL' four sets of intervals are created in the
#' database: \code{tss}, \code{exons}, \code{utr3} and \code{utr5}. See
#' \link{gintervals.import_genes} for more details about importing genes
#' intervals.
#'
#' 'fasta', 'genes.file' and 'annots.file' can be either a file path or URL in
#' a form of 'ftp://[address]/[file]'. 'fasta' can also contain wildcards to
#' indicate multiple files. Files that these arguments point to can be zipped
#' or unzipped.
#'
#' See the 'Genomes' vignette for details on how to create a database from common
#' genome sources.
#'
#' @param groot path to newly created database
#' @param fasta an array of names or URLs of FASTA files. Can contain wildcards
#' for multiple files
#' @param genes.file name or URL of file that contains genes. If 'NULL' no
#' genes are imported
#' @param annots.file name of URL file that contains annotations. If 'NULL' no
#' annotations are imported
#' @param annots.names annotations names
#' @param format database format: "indexed" (default, single genome.seq + genome.idx)
#' or "per-chromosome" (separate .seq file per contig). If NULL, uses the value from
#' \code{getOption("gmulticontig.indexed_format", TRUE)}
#' @param verbose if TRUE, prints verbose messages
#' @return None.
#' @seealso \code{\link{gdb.init}}, \code{\link{gdb.reload}},
#' \code{\link{gintervals.import_genes}}
#' @keywords ~database ~create ~genes
#' @examples
#' \donttest{
#' # ftp <- "ftp://hgdownload.soe.ucsc.edu/goldenPath/mm10"
#' # mm10_dir <- file.path(tempdir(), "mm10")
#' # # only a single chromosome is loaded in this example
#' # # see "Genomes" vignette how to download all of them and how
#' # # to download other genomes
#' # gdb.create(
#' #     mm10_dir,
#' #     paste(ftp, "chromosomes", paste0(
#' #         "chr", c("X"),
#' #         ".fa.gz"
#' #     ), sep = "/"),
#' #     paste(ftp, "database/knownGene.txt.gz", sep = "/"),
#' #     paste(ftp, "database/kgXref.txt.gz", sep = "/"),
#' #     c(
#' #         "kgID", "mRNA", "spID", "spDisplayID", "geneSymbol",
#' #         "refseq", "protAcc", "description", "rfamAcc",
#' #         "tRnaName"
#' #     )
#' # )
#' # gdb.init(mm10_dir)
#' # gintervals.ls()
#' # gintervals.all()
#' }
#'
#' @export gdb.create
