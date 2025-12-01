# Database creation functions
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
gdb.create <- function(groot = NULL, fasta = NULL, genes.file = NULL, annots.file = NULL, annots.names = NULL, format = NULL, verbose = FALSE) {
    if (is.null(groot) || is.null(fasta)) {
        stop("Usage: gdb.create(groot, fasta, genes.file = NULL, annots.file = NULL, annots.names = NULL, format = 'indexed', verbose = FALSE)", call. = FALSE)
    }

    # Determine format: parameter overrides option
    if (is.null(format)) {
        # Check for per-chromosome option for backward compatibility
        use_indexed_format <- getOption("gmulticontig.indexed_format", TRUE)
        format <- if (use_indexed_format) "indexed" else "per-chromosome"
    } else {
        # Validate format parameter
        format <- match.arg(format, c("indexed", "per-chromosome"))
    }
    use_indexed_format <- (format == "indexed")

    if (file.exists(groot)) {
        stop(sprintf("Directory %s already exists", groot), call. = FALSE)
    }

    success <- FALSE
    allgenome.old <- NULL
    groot.old <- NULL
    chrom_alias.old <- NULL
    db_per_chromosome.old <- NULL
    if (exists("ALLGENOME", envir = .misha)) {
        allgenome.old <- get("ALLGENOME", envir = .misha)
    }
    if (exists("GROOT", envir = .misha)) {
        groot.old <- get("GROOT", envir = .misha)
    }
    if (exists("CHROM_ALIAS", envir = .misha)) {
        chrom_alias.old <- get("CHROM_ALIAS", envir = .misha)
    }
    if (exists("DB_IS_PER_CHROMOSOME", envir = .misha)) {
        db_per_chromosome.old <- get("DB_IS_PER_CHROMOSOME", envir = .misha)
    }

    tryCatch(
        {
            assign("CHROM_ALIAS", NULL, envir = .misha)
            assign("DB_IS_PER_CHROMOSOME", FALSE, envir = .misha)
            dir.create(groot, showWarnings = FALSE, recursive = TRUE, mode = "0777")
            dir.create(paste(groot, "pssms", sep = "/"), showWarnings = FALSE, recursive = TRUE, mode = "0777")
            dir.create(paste(groot, "seq", sep = "/"), showWarnings = FALSE, recursive = TRUE, mode = "0777")
            dir.create(paste(groot, "tracks", sep = "/"), showWarnings = FALSE, recursive = TRUE, mode = "0777")

            # Detect import mode: multi-FASTA (indexed) vs per-chromosome
            is_single_file <- (length(fasta) == 1 && file.exists(fasta))

            chroms <- NULL
            chrom.sizes <- NULL

            if (use_indexed_format && is_single_file) {
                # Multi-FASTA import to indexed format
                # C++ function parses FASTA, creates genome.seq and genome.idx,
                # and returns a data frame with columns: name, size (sorted alphabetically)
                if (verbose) message("Creating indexed genome format...")
                contig_info <- .gseq.import_multifasta(groot, fasta, verbose)

                if (nrow(contig_info) == 0) {
                    stop("No contigs were imported from multi-FASTA file", call. = FALSE)
                }

                # Use the data frame directly (sorted by name to match genome.idx)
                chrom.sizes <- data.frame(chrom = contig_info$name, size = contig_info$size)
            } else {
                # Per-chromosome import
                if (verbose) message("Creating per-chromosome genome format...")
                chroms <- .gseq.import(groot, fasta)

                if (!length(chroms)) {
                    stop("No FASTA files were imported", call. = FALSE)
                }

                # Use chromosome names as-is (no chr prefix)
                seq.files <- paste(chroms, ".seq", sep = "")
                seq.files <- paste(paste(groot, "seq", sep = "/"), seq.files, sep = "/")
                chrom.sizes <- data.frame(chrom = chroms, size = file.info(seq.files)$size)
            }

            utils::write.table(chrom.sizes, paste(groot, "chrom_sizes.txt", sep = "/"), quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)

            # Compute chromosome aliases early so we can add them to factor levels
            canonical_names <- as.character(chrom.sizes$chrom)
            alias_map <- .compute_chrom_aliases(canonical_names)
            all_chrom_names <- unique(c(canonical_names, names(alias_map)))

            # before calling gintervals.import_genes new ALLGENOME must be set
            intervals <- data.frame(
                chrom = canonical_names,
                start = 0, end = as.numeric(chrom.sizes$size)
            )
            intervals$chrom <- factor(intervals$chrom, levels = all_chrom_names)
            # For indexed databases, chrom_sizes and intervals are already sorted by the C++ import
            # Don't sort again to maintain consistency with genome.idx chromid assignments
            # For per-chromosome databases, sort alphabetically for consistency
            if (!use_indexed_format) {
                intervals <- intervals[order(intervals$chrom), ]
            }
            rownames(intervals) <- 1:nrow(intervals)

            # Lazy 2D generation: only materialize for small genomes
            contig_threshold <- getOption("gmulticontig.2d.threshold", 100)

            if (nrow(intervals) <= contig_threshold) {
                # Small genome: materialize full 2D grid
                cartesian <- expand.grid(1:nrow(intervals), 1:nrow(intervals))
                intervals2d <- cbind(intervals[cartesian[, 2], ], intervals[cartesian[, 1], ])
                names(intervals2d) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2")
                # Ensure chrom1 and chrom2 have the same factor levels as intervals$chrom (including aliases)
                intervals2d$chrom1 <- factor(intervals2d$chrom1, levels = all_chrom_names)
                intervals2d$chrom2 <- factor(intervals2d$chrom2, levels = all_chrom_names)
                rownames(intervals2d) <- 1:nrow(intervals2d)
            } else {
                # Large genome: defer 2D generation
                intervals2d <- .create_deferred_2d(nrow(intervals))
                if (verbose) {
                    message(sprintf(
                        "Deferring 2D genome generation (%d contigs > threshold %d). Use gintervals.2d() for specific pairs.",
                        nrow(intervals), contig_threshold
                    ))
                }
            }

            assign("ALLGENOME", list(intervals, intervals2d), envir = .misha)
            assign("GROOT", groot, envir = .misha)
            .store_chrom_aliases(levels(intervals$chrom))

            if (!is.null(genes.file)) {
                intervs <- gintervals.import_genes(genes.file, annots.file, annots.names)
                if (!is.null(intervs)) {
                    for (i in 1:length(intervs)) {
                        if (!is.null(intervs$tss)) {
                            .gcall_noninteractive(.gintervals.save_file, sprintf("%s/tracks/%s.interv", groot, names(intervs)[i]), intervs[[i]])
                        }
                    }
                }
            }

            # write read-only attributes
            f <- file(paste(groot, ".ro_attributes", sep = "/"), "wb")
            serialize(c("created.by", "created.date", "created.user"), f)
            close(f)

            if (verbose) message("Database was successfully created")
            success <- TRUE
        },
        finally = {
            assign("ALLGENOME", allgenome.old, envir = .misha)
            assign("GROOT", groot.old, envir = .misha)
            assign("CHROM_ALIAS", chrom_alias.old, envir = .misha)
            assign("DB_IS_PER_CHROMOSOME", db_per_chromosome.old, envir = .misha)
            if (!success) {
                unlink(groot, recursive = TRUE)
            }
        }
    )
    retv <- 0 # suppress return value
}

#' Create and Load a Genome Database
#'
#' This function downloads, extracts, and loads a misha genome database for the specified genome.
#'
#' @param genome A character string specifying the genome to download. Supported genomes are "mm9", "mm10", "mm39", "hg19", and "hg38".
#' @param path A character string specifying the directory where the genome will be extracted. Defaults to genome name (e.g. "mm10") in the current working directory.
#' @param tmpdir A character string specifying the directory for storing temporary files. This is used for storing the downloaded genome file.
#'
#' @details
#' The function checks if the specified genome is available. If tmpdir, it constructs the download URL, downloads the genome file,
#' extracts it to the specified directory, and loads the genome database using \code{gsetroot}. The function also calls \code{gdb.reload} to reload the genome database.
#'
#' @return None.
#'
#' @examples
#' \donttest{
#' mm10_dir <- tempdir()
#' gdb.create_genome("mm10", path = mm10_dir)
#' list.files(file.path(mm10_dir, "mm10"))
#' gsetroot(file.path(mm10_dir, "mm10"))
#' gintervals.ls()
#' }
#'
#' @export
gdb.create_genome <- function(genome, path = getwd(), tmpdir = tempdir()) {
    supported_genomes <- c("mm9", "mm10", "mm39", "hg19", "hg38")

    if (!genome %in% supported_genomes) {
        stop(paste("The genome", genome, "is not available yet. Available genomes are:", paste(supported_genomes, collapse = ", ")))
    }

    base_url <- "https://misha-genome.s3.eu-west-1.amazonaws.com/"
    url <- paste0(base_url, genome, ".tar.gz")

    temp_file <- tempfile(fileext = ".tar.gz", tmpdir = tmpdir)

    withr::local_options(list(timeout = 60 * 60 * 2))
    message("Downloading ", genome, " genome...")
    utils::download.file(url, temp_file, mode = "wb")

    dir.create(path, showWarnings = FALSE, recursive = TRUE)

    message("Extracting ", genome, " genome...")
    utils::untar(temp_file, exdir = path)

    unlink(temp_file)

    message("Loading misha root...")
    gdb.init(file.path(path, genome), rescan = TRUE)

    message(genome, " genome has been successfully downloaded and extracted to ", file.path(path, genome))
}


#' Returns a list of read-only track attributes
#'
#' Returns a list of read-only track attributes.
#'
#' This function returns a list of read-only track attributes. These attributes
#' are not allowed to be modified or deleted.
#'
#' If no attributes are marked as read-only a 'NULL' is returned.
#'
#' @return A list of read-only track attributes.
#' @seealso \code{\link{gdb.set_readonly_attrs}},
#' \code{\link{gtrack.attr.get}}, \code{\link{gtrack.attr.set}}
#' @keywords ~attr ~attribute
#' @export gdb.get_readonly_attrs
gdb.get_readonly_attrs <- function() {
    .gcheckroot()

    filename <- paste(get("GROOT", envir = .misha), ".ro_attributes", sep = "/")
    attrs <- NULL
    if (file.exists(filename)) {
        f <- file(filename, "rb")
        attrs <- unserialize(f)
        close(f)
        if (!is.character(attrs)) {
            stop(sprintf("Invalid format of read-only atrributes file %s", filename), call. = FALSE)
        }

        attrs <- unique(attrs)
        attrs <- attrs[attrs != ""]
    }
    attrs
}


#' Initializes connection with Genomic Database
#'
#' Initializes connection with Genomic Database: loads the list of tracks,
#' intervals, etc.
#'
#' 'gdb.init' initializes the connection with the Genomic Database. It is
#' typically called first prior to any other function. When the package is
#' attached it internally calls to 'gdb.init.examples' which opens the
#' connection with the database located at 'PKGDIR/trackdb/test' directory,
#' where 'PKGDIR' is the directory where the package is installed.
#'
#' The current working directory inside the Genomic Database is set to 'dir'.
#' If 'dir' is 'NULL', the current working directory is set to 'GROOT/tracks'.
#'
#' If 'rescan' is 'TRUE', the list of tracks and intervals is achieved by
#' rescanning directory structure under the current current working directory.
#' Otherwise 'gdb.init' attempts to use the cached list that resides in
#' 'groot/.db.cache' file.
#'
#' Upon completion the connection is established with the database. If
#' auto-completion mode is switched on (see 'gset_input_method') the list of
#' tracks and intervals sets is loaded and added as variables to the global
#' environment allowing auto-completion of object names with <TAB> key. Also a
#' few variables are defined at an environment called \code{.misha}, and can be
#' accessed using \code{.misha$variable}, e.g. \code{.misha$ALLGENOME}.
#' These variables should not be modified by user.
#'
#' \tabular{ll}{ GROOT \tab Root directory of Genomic Database\cr GWD \tab
#' Current working directory inside Genomic Database\cr GTRACKS \tab List of
#' all available tracks\cr GINTERVS \tab List of all available intervals\cr
#' GVTRACKS \tab List of all available virtual tracks\cr ALLGENOME \tab List of
#' all chromosomes and their sizes\cr GITERATOR.INTERVALS \tab A set of
#' iterator intervals for which the track expression is evaluated\cr }
#'
#' When option 'gmulticontig.indexed_format' is set to TRUE, the function
#' loads a database with "indexed" track format.
#'
#' @aliases gdb.init gdb.init.examples gsetroot
#' @param groot the root directory of the Genomic Database
#' @param dir the current working directory inside the Genomic Database
#' @param rescan indicates whether the file structure should be rescanned
#' @return None.
#' @seealso \code{\link{gdb.reload}}, \code{\link{gdb.create}},
#' \code{\link{gdir.cd}}, \code{\link{gtrack.ls}}, \code{\link{gintervals.ls}},
#' \code{\link{gvtrack.ls}}
#' @keywords ~db ~data ~database
#' @export gdb.init
gdb.init <- function(groot = NULL, dir = NULL, rescan = FALSE) {
    if (is.null(groot)) {
        stop("Usage: gdb.init(groot, dir = NULL, rescan = FALSE)", call. = FALSE)
    }
    gsetroot(groot, dir, rescan)
}

#' @rdname gdb.init
#' @export
gdb.init_examples <- function() {
    db_dir <- tempdir()
    test_path <- file.path(db_dir, "trackdb/test")
    unlink(test_path, recursive = TRUE)
    utils::untar(system.file("testdb.tar.gz", package = "misha"), exdir = db_dir)
    gsetroot(test_path)
    if (getOption("gmulticontig.indexed_format", FALSE)) {
        gdb.convert_to_indexed(test_path, convert_tracks = TRUE, remove_old_files = TRUE, convert_intervals = TRUE, verbose = FALSE, force = TRUE)
        gsetroot(test_path)
    }
}

#' Reloads database from the disk
#'
#' Reloads database from disk: list of tracks, intervals, etc.
#'
#' Reloads Genomic Database from disk: list of tracks, intervals, etc. Use this
#' function if you manually add tracks or if for any reason the database
#' becomes corrupted. If 'rescan' is 'TRUE', the list of tracks and intervals
#' is achieved by rescanning directory structure under the current current
#' working directory. Otherwise 'gdb.reload' attempts to use the cached list
#' that resides in 'GROOT/.db.cache' file.
#'
#' @param rescan indicates whether the file structure should be rescanned
#' @seealso \code{\link{gdb.init}}, \code{\link{gdb.create}},
#' \code{\link{gdir.cd}},
#' @keywords ~db
#' @return No return value, called for side effects.
#' @export gdb.reload
