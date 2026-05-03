# Build a misha genome from a name.
#
# Public surface (exported):
#   gdb.build_genome()           — build from registry-resolved recipe
#   gdb.list_genomes()           — list resolvable genome names
#   gdb.genome_info()            — show resolved recipe without building
#   gdb.install_gff3_converter() — pre-install UCSC's gff3ToGenePred binary
#
# Spec: dev/notes/specs/2026-04-30_genome-build-design.md

# ---------------------------------------------------------------------------
# Internal constants
# ---------------------------------------------------------------------------

.MISHA_GENOME_SOURCES <- c("ucsc", "ncbi", "s3", "manual", "local")

.UCSC_GOLDENPATH <- "https://hgdownload.soe.ucsc.edu/goldenPath"
.NCBI_DATASETS_API <- "https://api.ncbi.nlm.nih.gov/datasets/v2"

# NCBI backend chromosome naming options. See ?gdb.build_genome.
.NCBI_CHROM_NAMING_VALUES <- c("sequence_name", "ucsc", "accession")
.NCBI_DEFAULT_CHROM_NAMING <- "sequence_name"

# RefSeqLink columns we expose as gene annotations (used by UCSC backend).
# Schema: 19 columns as ncbiRefSeqLink.txt.gz ships them. The C++ importer
# requires the count to match exactly.
.UCSC_NCBI_REFSEQ_LINK_COLS <- c(
    "id", "status", "name", "product", "mrnaAcc", "protAcc", "locusLinkId",
    "omimId", "hgnc", "genbank", "pseudo", "gbkey", "source", "gene_biotype",
    "gene_synonym", "ncrna_class", "note", "description", "externalId"
)

# ---------------------------------------------------------------------------
# Registry parsing and resolution
# ---------------------------------------------------------------------------

# Normalize a raw registry entry to a recipe list with $source set.
# Accepts:
#   - character scalar  -> {source: local, path: <string>}  (legacy form)
#   - named list with $source
.normalize_recipe <- function(entry, name) {
    if (is.character(entry) && length(entry) == 1) {
        return(list(source = "local", path = entry))
    }
    if (!is.list(entry)) {
        stop(sprintf(
            "Genome '%s' has invalid registry entry: expected string or mapping, got %s",
            name, class(entry)[[1]]
        ), call. = FALSE)
    }
    if (is.null(entry$source) || !is.character(entry$source) || length(entry$source) != 1) {
        stop(sprintf("Genome '%s' is missing required 'source:' field", name), call. = FALSE)
    }
    if (!entry$source %in% .MISHA_GENOME_SOURCES) {
        stop(sprintf(
            "Genome '%s' has unknown source '%s'. Valid sources: %s",
            name, entry$source, paste(.MISHA_GENOME_SOURCES, collapse = ", ")
        ), call. = FALSE)
    }
    entry
}

.validate_recipe <- function(recipe, name) {
    src <- recipe$source
    miss <- function(field) {
        stop(sprintf("Genome '%s' (source: %s) is missing required field '%s'", name, src, field), call. = FALSE)
    }
    if (src == "ucsc") {
        if (is.null(recipe$assembly)) miss("assembly")
    } else if (src == "ncbi") {
        if (is.null(recipe$accession)) miss("accession")
        if (!grepl("^GC[FA]_[0-9]+\\.[0-9]+$", recipe$accession)) {
            stop(sprintf(
                "Genome '%s': accession '%s' does not match GC[FA]_<digits>.<digits>",
                name, recipe$accession
            ), call. = FALSE)
        }
        if (!is.null(recipe$chrom_naming) &&
            !recipe$chrom_naming %in% .NCBI_CHROM_NAMING_VALUES) {
            stop(sprintf(
                "Genome '%s': chrom_naming '%s' invalid. Valid values: %s",
                name, recipe$chrom_naming,
                paste(.NCBI_CHROM_NAMING_VALUES, collapse = ", ")
            ), call. = FALSE)
        }
    } else if (src == "s3") {
        if (is.null(recipe$assembly)) miss("assembly")
    } else if (src == "local") {
        if (is.null(recipe$path)) miss("path")
    } else if (src == "manual") {
        if (is.null(recipe$fasta)) miss("fasta")
    }
    invisible(recipe)
}

# Read a single registry YAML file. Returns a named list <name> -> recipe (raw,
# not yet validated).
.parse_genome_registry <- function(path) {
    if (!file.exists(path)) {
        stop(sprintf("Registry file does not exist: %s", path), call. = FALSE)
    }
    y <- tryCatch(
        yaml::read_yaml(path),
        error = function(e) stop(sprintf("Failed to parse YAML registry %s: %s", path, conditionMessage(e)), call. = FALSE)
    )
    if (is.null(y$genome)) {
        return(list())
    }
    if (!is.list(y$genome) || is.null(names(y$genome))) {
        stop(sprintf("Registry %s: 'genome:' must be a named mapping", path), call. = FALSE)
    }
    y$genome
}

# Walk up from getwd() to git root looking for a misha.yaml.
.find_project_misha_yaml <- function() {
    dir <- normalizePath(getwd(), mustWork = FALSE)
    while (TRUE) {
        candidate <- file.path(dir, "misha.yaml")
        if (file.exists(candidate)) {
            return(candidate)
        }
        if (file.exists(file.path(dir, ".git"))) {
            return(NULL)
        }
        parent <- dirname(dir)
        if (parent == dir) {
            return(NULL)
        }
        dir <- parent
    }
}

.builtin_registry_path <- function() {
    system.file("genomes.yaml", package = "misha")
}

# Resolve `name` through the chain. Returns a list:
#   list(recipe = <list>, resolved_from = <character>)
.resolve_genome <- function(name, registry = NULL) {
    if (!is.character(name) || length(name) != 1 || !nzchar(name)) {
        stop("name must be a non-empty string", call. = FALSE)
    }

    sources <- list()
    if (!is.null(registry)) {
        sources[[length(sources) + 1]] <- list(label = sprintf("registry arg (%s)", registry), path = registry)
    }
    opt <- getOption("misha.genome_registry")
    if (!is.null(opt)) {
        sources[[length(sources) + 1]] <- list(label = sprintf("getOption('misha.genome_registry') (%s)", opt), path = opt)
    }
    proj <- .find_project_misha_yaml()
    if (!is.null(proj)) {
        sources[[length(sources) + 1]] <- list(label = sprintf("project misha.yaml (%s)", proj), path = proj)
    }
    builtin <- .builtin_registry_path()
    if (nzchar(builtin)) {
        sources[[length(sources) + 1]] <- list(label = "built-in (inst/genomes.yaml)", path = builtin)
    }

    for (src in sources) {
        entries <- tryCatch(.parse_genome_registry(src$path), error = function(e) {
            stop(sprintf("Error reading %s: %s", src$label, conditionMessage(e)), call. = FALSE)
        })
        if (name %in% names(entries)) {
            recipe <- .normalize_recipe(entries[[name]], name)
            .validate_recipe(recipe, name)
            return(list(recipe = recipe, resolved_from = src$label))
        }
    }

    # Pattern fallback for NCBI accessions.
    if (grepl("^GC[FA]_[0-9]+\\.[0-9]+$", name)) {
        recipe <- list(source = "ncbi", accession = name)
        return(list(recipe = recipe, resolved_from = "pattern fallback (NCBI accession)"))
    }

    stop(sprintf(
        "Genome '%s' not found in any registry. Searched: %s. To define it, add an entry to a misha.yaml or use gdb.create() directly.",
        name, paste(vapply(sources, `[[`, character(1), "label"), collapse = "; ")
    ), call. = FALSE)
}

# ---------------------------------------------------------------------------
# Helpers for downloads + post-build annotation loading
# ---------------------------------------------------------------------------

# Download URL to a destination (binary mode). Returns dest path.
.download_to <- function(url, dest, verbose = TRUE) {
    if (verbose) message(sprintf("Downloading %s ...", url))
    utils::download.file(url, dest, mode = "wb", quiet = !verbose)
    dest
}

# Heal UCSC TSV escape artifacts in description/note fields:
#   - backslash + (CR? + LF)  -> space    (line continuation)
#   - backslash + tab          -> space    (intra-field tab escape)
#   - stray CR (carriage return) -> nothing  (Windows line endings R's
#     readLines() would otherwise treat as line terminators)
#
# Loads the whole file into memory; UCSC tables are small (a few MB) so
# this is fine in practice and lets us do the regex fixups on the full
# string without chunk-boundary worries.
.heal_ucsc_tsv_escapes <- function(in_path, out_path) {
    con_in <- if (grepl("\\.gz$", in_path)) gzfile(in_path, "rb") else file(in_path, "rb")
    on.exit(close(con_in), add = TRUE)
    bytes <- raw(0)
    repeat {
        chunk <- readBin(con_in, raw(), n = 1024L * 1024L)
        if (!length(chunk)) break
        bytes <- c(bytes, chunk)
    }
    s <- rawToChar(bytes)
    # Order matters: do CRLF before LF so we consume the whole sequence.
    s <- gsub("\\\\\r\n", " ", s, perl = TRUE)
    s <- gsub("\\\\\n", " ", s, perl = TRUE)
    s <- gsub("\\\\\t", " ", s, perl = TRUE)
    s <- gsub("\r", "", s, perl = TRUE, fixed = FALSE)
    con_out <- file(out_path, "wb")
    on.exit(close(con_out), add = TRUE)
    writeBin(charToRaw(s), con_out)
    out_path
}

# Normalize a UCSC TSV table to exactly N columns per row. UCSC sometimes
# ships rows with stray embedded tabs (e.g. in description fields) — these
# lines have NF != expected count and would crash the C++ importer's strict
# column-count check. Short rows are padded with empty strings; long rows
# have their trailing extras joined back into the last column with " ".
.normalize_ucsc_tsv <- function(in_path, out_path, n_cols) {
    con_in <- if (grepl("\\.gz$", in_path)) gzfile(in_path, "rt") else file(in_path, "rt")
    on.exit(close(con_in), add = TRUE)
    con_out <- file(out_path, "wt")
    on.exit(close(con_out), add = TRUE)
    repeat {
        lines <- readLines(con_in, n = 50000L, warn = FALSE)
        if (!length(lines)) break
        fields <- strsplit(lines, "\t", fixed = TRUE)
        out <- vapply(fields, function(x) {
            if (length(x) == n_cols) {
                paste(x, collapse = "\t")
            } else if (length(x) < n_cols) {
                paste(c(x, rep("", n_cols - length(x))), collapse = "\t")
            } else {
                # Join overflow into the last column so the row has exactly n_cols fields.
                paste(c(x[seq_len(n_cols - 1)], paste(x[n_cols:length(x)], collapse = " ")),
                    collapse = "\t"
                )
            }
        }, character(1))
        writeLines(out, con_out)
    }
    out_path
}

# Trim a genePred-format file to the classic 12-col layout the C++
# gintervals_import_genes expects. Two input shapes are handled:
#   - 16 cols: UCSC extended genePred with leading bin column (ncbiRefSeq,
#     refGene, knownGene from goldenPath/database/). Take cols 2-13.
#   - 15 cols: extended genePred without bin (gff3ToGenePred output, NCBI).
#     Take cols 1-12.
# In both cases the resulting 12 cols are:
#   [name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount,
#    exonStarts, exonEnds, score, name2]
# Cols 11-12 occupy the C++ importer's PROTEINID/ALIGNID slots — read but
# unused, so populating with score/name2 is harmless.
#
# Accepts .txt or .txt.gz, writes a .txt with 12-col content. Streams the
# input so memory stays bounded.
.normalize_ucsc_genepred <- function(in_path, out_path) {
    con_in <- if (grepl("\\.gz$", in_path)) gzfile(in_path, "rt") else file(in_path, "rt")
    on.exit(close(con_in), add = TRUE)
    con_out <- file(out_path, "wt")
    on.exit(close(con_out), add = TRUE)

    chunk_size <- 50000L
    take_range <- NULL
    repeat {
        lines <- readLines(con_in, n = chunk_size, warn = FALSE)
        if (!length(lines)) break
        fields <- strsplit(lines, "\t", fixed = TRUE)
        nf <- vapply(fields, length, integer(1))
        if (is.null(take_range)) {
            shape <- nf[[1]]
            take_range <- if (shape == 16L) {
                2:13
            } else if (shape == 15L) {
                1:12
            } else if (shape == 12L) {
                1:12
            } else {
                stop(sprintf(
                    "genePred file %s: unsupported column count %d (expected 12, 15, or 16)",
                    in_path, shape
                ), call. = FALSE)
            }
        }
        bad <- nf < max(take_range)
        if (any(bad)) {
            stop(sprintf(
                "genePred file %s: row %d has only %d fields, need at least %d",
                in_path, which(bad)[1], nf[which(bad)[1]], max(take_range)
            ), call. = FALSE)
        }
        trimmed <- vapply(
            fields, function(x) paste(x[take_range], collapse = "\t"),
            character(1)
        )
        writeLines(trimmed, con_out)
    }
    out_path
}

# Decompress a .gz file (file -> file with .gz removed). Returns the new path.
.gunzip_to_file <- function(gz_path, out_path = sub("\\.gz$", "", gz_path)) {
    con_in <- gzfile(gz_path, "rb")
    on.exit(close(con_in), add = TRUE)
    con_out <- file(out_path, "wb")
    on.exit(close(con_out), add = TRUE)
    repeat {
        chunk <- readBin(con_in, raw(), n = 1024 * 1024)
        if (length(chunk) == 0) break
        writeBin(chunk, con_out)
    }
    out_path
}

# Parse UCSC rmsk.txt(.gz) — 17 columns. Returns a data.frame with intervals
# columns plus name/class/family. Strand is normalized to numeric (1/-1/0).
.parse_ucsc_rmsk <- function(file) {
    cols <- c(
        "bin", "swScore", "milliDiv", "milliDel", "milliIns",
        "genoName", "genoStart", "genoEnd", "genoLeft", "strand",
        "repName", "repClass", "repFamily", "repStart", "repEnd",
        "repLeft", "id"
    )
    df <- utils::read.table(file,
        sep = "\t", header = FALSE, col.names = cols,
        quote = "", comment.char = "", stringsAsFactors = FALSE,
        na.strings = character(0)
    )
    data.frame(
        chrom = df$genoName,
        start = df$genoStart,
        end = df$genoEnd,
        strand = ifelse(df$strand == "+", 1L, ifelse(df$strand == "-", -1L, 0L)),
        name = df$repName,
        class = df$repClass,
        family = df$repFamily,
        stringsAsFactors = FALSE
    )
}

.parse_ucsc_cpg_island <- function(file) {
    cols <- c(
        "bin", "chrom", "chromStart", "chromEnd", "name",
        "length", "cpgNum", "gcNum", "perCpg", "perGc", "obsExp"
    )
    df <- utils::read.table(file,
        sep = "\t", header = FALSE, col.names = cols,
        quote = "", comment.char = "", stringsAsFactors = FALSE,
        na.strings = character(0)
    )
    data.frame(
        chrom = df$chrom,
        start = df$chromStart,
        end = df$chromEnd,
        name = df$name,
        length = df$length,
        cpgNum = df$cpgNum,
        perCpg = df$perCpg,
        perGc = df$perGc,
        obsExp = df$obsExp,
        stringsAsFactors = FALSE
    )
}

.parse_ucsc_cytoband <- function(file) {
    cols <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")
    df <- utils::read.table(file,
        sep = "\t", header = FALSE, col.names = cols,
        quote = "", comment.char = "", stringsAsFactors = FALSE,
        na.strings = character(0)
    )
    data.frame(
        chrom = df$chrom,
        start = df$chromStart,
        end = df$chromEnd,
        name = df$name,
        stain = df$gieStain,
        stringsAsFactors = FALSE
    )
}

# Save a parsed-table data.frame as a misha intervals set (post-build,
# requires gdb.init has already been called on the new groot).
# Drops rows whose chrom is not in ALLGENOME (different sources sometimes
# include haplotype/alt contigs not present in the assembled FASTA).
.save_post_build_intervals <- function(intervals_set_name, df, verbose = TRUE) {
    if (!exists("ALLGENOME", envir = .misha)) {
        stop(".save_post_build_intervals called without an active GROOT", call. = FALSE)
    }
    chroms <- as.character(get("ALLGENOME", envir = .misha)[[1]]$chrom)
    keep <- df$chrom %in% chroms
    if (any(!keep) && verbose) {
        message(sprintf(
            "  %s: dropped %d/%d rows on contigs not in ALLGENOME",
            intervals_set_name, sum(!keep), nrow(df)
        ))
    }
    df <- df[keep, , drop = FALSE]
    if (nrow(df) == 0) {
        if (verbose) message(sprintf("  %s: no rows remained after filtering, skipping", intervals_set_name))
        return(invisible(NULL))
    }
    df$chrom <- factor(df$chrom, levels = chroms)
    gintervals.save(intervals_set_name, df)
    if (verbose) message(sprintf("  %s: saved %d intervals", intervals_set_name, nrow(df)))
    invisible(NULL)
}

# Write genome_info.yaml — a record of where this groot came from.
.write_genome_info <- function(groot, name, recipe, annotations, files = list()) {
    info <- list(
        name = name,
        source = recipe$source,
        downloaded_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
        misha_version = as.character(utils::packageVersion("misha")),
        annotations = as.list(annotations),
        recipe = recipe,
        files = files
    )
    yaml::write_yaml(info, file.path(groot, "genome_info.yaml"))
    invisible(NULL)
}

# ---------------------------------------------------------------------------
# UCSC backend
# ---------------------------------------------------------------------------

# Returns a named list of URLs for an assembly + requested annotations.
.ucsc_urls <- function(assembly, annotations, single_fasta = TRUE) {
    base <- sprintf("%s/%s", .UCSC_GOLDENPATH, assembly)
    out <- list()
    if (single_fasta) {
        out$fasta <- sprintf("%s/bigZips/%s.fa.gz", base, assembly)
    } # per-chromosome list is built dynamically (see .ucsc_per_chrom_fasta)
    if ("genes" %in% annotations) {
        out$genes <- sprintf("%s/database/ncbiRefSeq.txt.gz", base)
        out$annots <- sprintf("%s/database/ncbiRefSeqLink.txt.gz", base)
    }
    if ("rmsk" %in% annotations) {
        out$rmsk <- sprintf("%s/database/rmsk.txt.gz", base)
    }
    if ("cpgIsland" %in% annotations) {
        out$cpgIsland <- sprintf("%s/database/cpgIslandExt.txt.gz", base)
    }
    if ("cytoband" %in% annotations) {
        out$cytoband <- sprintf("%s/database/cytoBandIdeo.txt.gz", base)
    }
    out
}

.run_ucsc_backend <- function(name, recipe, path, annotations, format, verbose) {
    assembly <- recipe$assembly
    urls <- .ucsc_urls(assembly, annotations, single_fasta = TRUE)

    if (verbose) message(sprintf("UCSC backend: assembly=%s -> %s", assembly, path))

    # Pre-download FASTA + genes + annots locally. The underlying gdb.create()
    # / .gseq.import path only handles ftp:// or local files; passing https://
    # URLs straight through would fail. By downloading first, we get
    # protocol-agnostic behavior and let gdb.create() use its indexed-format
    # multi-FASTA fast path (which requires file.exists()).
    workdir <- tempfile("misha_ucsc_")
    dir.create(workdir, recursive = TRUE)
    on.exit(unlink(workdir, recursive = TRUE), add = TRUE)

    local_fasta <- file.path(workdir, sprintf("%s.fa.gz", assembly))
    .download_to(urls$fasta, local_fasta, verbose = verbose)

    local_genes <- NULL
    if (!is.null(urls$genes)) {
        raw_genes <- file.path(workdir, "ncbiRefSeq.txt.gz")
        .download_to(urls$genes, raw_genes, verbose = verbose)
        local_genes <- file.path(workdir, "ncbiRefSeq.12col.txt")
        if (verbose) message("  Trimming extended-genePred (16 cols) to classic 12-col format ...")
        .normalize_ucsc_genepred(raw_genes, local_genes)
    }
    local_annots <- NULL
    annots_names <- NULL
    if (!is.null(urls$annots)) {
        raw_annots <- file.path(workdir, "ncbiRefSeqLink.raw.txt.gz")
        .download_to(urls$annots, raw_annots, verbose = verbose)
        healed <- file.path(workdir, "ncbiRefSeqLink.healed.txt")
        if (verbose) message("  Healing UCSC escape artifacts in ncbiRefSeqLink ...")
        .heal_ucsc_tsv_escapes(raw_annots, healed)
        local_annots <- file.path(workdir, "ncbiRefSeqLink.normalized.txt")
        if (verbose) message("  Normalizing ncbiRefSeqLink to fixed column count ...")
        .normalize_ucsc_tsv(healed, local_annots,
            n_cols = length(.UCSC_NCBI_REFSEQ_LINK_COLS)
        )
        annots_names <- .UCSC_NCBI_REFSEQ_LINK_COLS
    }

    gdb.create(
        groot        = path,
        fasta        = local_fasta,
        genes.file   = local_genes,
        annots.file  = local_annots,
        annots.names = annots_names,
        format       = format,
        verbose      = verbose
    )

    # Switch to the new groot for post-build annotation loading.
    gdb.init(path, rescan = TRUE)

    files_record <- list(
        fasta = list(url = urls$fasta),
        genes = list(url = urls$genes %||% NA_character_),
        annots = list(url = urls$annots %||% NA_character_)
    )

    if ("rmsk" %in% annotations && !is.null(urls$rmsk)) {
        tryCatch(
            {
                tmp <- tempfile(fileext = ".gz")
                on.exit(unlink(tmp), add = TRUE)
                .download_to(urls$rmsk, tmp, verbose = verbose)
                unz <- .gunzip_to_file(tmp, sub("\\.gz$", "", tmp))
                df <- .parse_ucsc_rmsk(unz)
                .save_post_build_intervals("rmsk", df, verbose = verbose)
                files_record$rmsk <- list(url = urls$rmsk)
            },
            error = function(e) {
                warning(sprintf("rmsk load failed for %s: %s", assembly, conditionMessage(e)), call. = FALSE)
            }
        )
    }

    if ("cpgIsland" %in% annotations && !is.null(urls$cpgIsland)) {
        tryCatch(
            {
                tmp <- tempfile(fileext = ".gz")
                on.exit(unlink(tmp), add = TRUE)
                .download_to(urls$cpgIsland, tmp, verbose = verbose)
                unz <- .gunzip_to_file(tmp, sub("\\.gz$", "", tmp))
                df <- .parse_ucsc_cpg_island(unz)
                .save_post_build_intervals("cpgIsland", df, verbose = verbose)
                files_record$cpgIsland <- list(url = urls$cpgIsland)
            },
            error = function(e) {
                warning(sprintf("cpgIsland load failed for %s: %s", assembly, conditionMessage(e)), call. = FALSE)
            }
        )
    }

    if ("cytoband" %in% annotations && !is.null(urls$cytoband)) {
        tryCatch(
            {
                tmp <- tempfile(fileext = ".gz")
                on.exit(unlink(tmp), add = TRUE)
                .download_to(urls$cytoband, tmp, verbose = verbose)
                unz <- .gunzip_to_file(tmp, sub("\\.gz$", "", tmp))
                df <- .parse_ucsc_cytoband(unz)
                .save_post_build_intervals("cytoband", df, verbose = verbose)
                files_record$cytoband <- list(url = urls$cytoband)
            },
            error = function(e) {
                warning(sprintf("cytoband load failed for %s: %s", assembly, conditionMessage(e)), call. = FALSE)
            }
        )
    }

    .write_genome_info(path, name, recipe, annotations, files_record)
    invisible(NULL)
}

# `%||%` is defined in db-core.R (shared utility).

# ---------------------------------------------------------------------------
# NCBI Datasets backend
# ---------------------------------------------------------------------------

.ncbi_datasets_zip_url <- function(accession) {
    sprintf(
        "%s/genome/accession/%s/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,SEQUENCE_REPORT",
        .NCBI_DATASETS_API, accession
    )
}

# Parse NCBI Datasets sequence_report.jsonl into a data.frame.
# Returns: data.frame(refseqAccession, genbankAccession, chrName, sequenceName,
# role, length).
.parse_ncbi_sequence_report <- function(path) {
    if (!file.exists(path)) {
        stop(sprintf("Sequence report not found: %s", path), call. = FALSE)
    }
    lines <- readLines(path, warn = FALSE)
    lines <- lines[nzchar(lines)]
    parsed <- lapply(lines, function(l) {
        # Tiny JSON parser via yaml (yaml is a JSON superset; already a dep).
        yaml::yaml.load(l)
    })
    pull <- function(field, default = NA_character_) {
        vapply(parsed, function(r) {
            v <- r[[field]]
            if (is.null(v) || !length(v)) default else as.character(v)[[1]]
        }, character(1))
    }
    pull_int <- function(field) {
        vapply(parsed, function(r) {
            v <- r[[field]]
            if (is.null(v) || !length(v)) NA_real_ else as.numeric(v)[[1]]
        }, numeric(1))
    }
    data.frame(
        refseqAccession  = pull("refseqAccession"),
        genbankAccession = pull("genbankAccession"),
        chrName          = pull("chrName"),
        sequenceName     = pull("sequenceName"),
        role             = pull("role"),
        length           = pull_int("length"),
        stringsAsFactors = FALSE
    )
}

# UCSC-style canonical name from a sequence_report row.
# Assembled molecules: chr1, chrX, chrM (UCSC uses chrM not chrMT).
# Unplaced: chrUn_<accession_underscore_dot_to_v>, e.g. NW_026256937.1 -> chrUn_NW026256937v1.
# Unlocalized: chr<chrom>_<acc>_random.
.ncbi_to_ucsc_name <- function(refseq_acc, chr_name, role) {
    # Build the UCSC-friendly suffix from an accession: drop underscore + dot,
    # encode version as v<N>: "NW_026256937.1" -> "NW026256937v1".
    encode_acc <- function(acc) {
        s <- gsub("_", "", acc, fixed = TRUE)
        s <- sub("\\.([0-9]+)$", "v\\1", s)
        s
    }
    if (role == "assembled-molecule") {
        if (chr_name %in% c("MT", "M")) {
            return("chrM")
        }
        return(paste0("chr", chr_name))
    }
    if (role == "unlocalized-scaffold") {
        return(paste0("chr", chr_name, "_", encode_acc(refseq_acc), "_random"))
    }
    paste0("chrUn_", encode_acc(refseq_acc))
}

# Build a named character vector: original FASTA/GFF id -> target canonical
# chrom name, given the desired chrom_naming and a parsed sequence report.
# Always indexed by the refseqAccession used in the FASTA/GFF (the field NCBI
# uses as the seqid).
.build_ncbi_rename_map <- function(seqrep, chrom_naming) {
    chrom_naming <- match.arg(chrom_naming, .NCBI_CHROM_NAMING_VALUES)
    targets <- if (chrom_naming == "accession") {
        seqrep$refseqAccession
    } else if (chrom_naming == "sequence_name") {
        # Use NCBI chrName for assembled molecules; refseq accession for
        # everything else (chrName="Un" would collide).
        ifelse(seqrep$role == "assembled-molecule",
            seqrep$chrName,
            seqrep$refseqAccession
        )
    } else { # ucsc
        mapply(.ncbi_to_ucsc_name,
            seqrep$refseqAccession, seqrep$chrName, seqrep$role,
            USE.NAMES = FALSE
        )
    }
    if (anyDuplicated(targets)) {
        dups <- unique(targets[duplicated(targets)])
        stop(sprintf(
            "chrom_naming='%s' produced duplicate names: %s",
            chrom_naming, paste(utils::head(dups, 5), collapse = ", ")
        ), call. = FALSE)
    }
    setNames(targets, seqrep$refseqAccession)
}

# Stream-rewrite FASTA: replace headers ">acc ..." with ">new_name". Works
# on plain or .gz input; output is plain.
.rename_fasta_headers <- function(in_path, out_path, rename_map, verbose = TRUE) {
    con_in <- if (grepl("\\.gz$", in_path)) gzfile(in_path, "rt") else file(in_path, "rt")
    on.exit(close(con_in), add = TRUE)
    con_out <- file(out_path, "wt")
    on.exit(close(con_out), add = TRUE)
    n_renamed <- 0L
    n_unmapped <- 0L
    repeat {
        lines <- readLines(con_in, n = 100000L, warn = FALSE)
        if (!length(lines)) break
        is_header <- startsWith(lines, ">")
        if (any(is_header)) {
            headers <- lines[is_header]
            ids <- sub("^>([^[:space:]]+).*$", "\\1", headers, perl = TRUE)
            new <- rename_map[ids]
            unmapped <- is.na(new)
            new[unmapped] <- ids[unmapped]
            n_renamed <- n_renamed + sum(!unmapped)
            n_unmapped <- n_unmapped + sum(unmapped)
            lines[is_header] <- paste0(">", new)
        }
        writeLines(lines, con_out)
    }
    if (verbose) {
        message(sprintf(
            "  Renamed %d FASTA contigs (%d unmapped, kept original id).",
            n_renamed, n_unmapped
        ))
    }
    invisible(out_path)
}

# Stream-rewrite GFF3: replace seqid (col 1) using rename_map. Comment lines
# preserved unchanged. Lines whose seqid isn't in the map are kept as-is.
.rename_gff3_seqids <- function(in_path, out_path, rename_map, verbose = TRUE) {
    con_in <- if (grepl("\\.gz$", in_path)) gzfile(in_path, "rt") else file(in_path, "rt")
    on.exit(close(con_in), add = TRUE)
    con_out <- file(out_path, "wt")
    on.exit(close(con_out), add = TRUE)
    n_renamed <- 0L
    n_unmapped <- 0L
    repeat {
        lines <- readLines(con_in, n = 100000L, warn = FALSE)
        if (!length(lines)) break
        is_data <- nzchar(lines) & !startsWith(lines, "#")
        if (any(is_data)) {
            data_lines <- lines[is_data]
            tab_pos <- regexpr("\t", data_lines, fixed = TRUE)
            seqids <- ifelse(tab_pos > 0, substring(data_lines, 1, tab_pos - 1), data_lines)
            rest <- ifelse(tab_pos > 0, substring(data_lines, tab_pos), "")
            new <- rename_map[seqids]
            unmapped <- is.na(new)
            new[unmapped] <- seqids[unmapped]
            n_renamed <- n_renamed + sum(!unmapped)
            n_unmapped <- n_unmapped + sum(unmapped)
            lines[is_data] <- paste0(new, rest)
        }
        writeLines(lines, con_out)
    }
    if (verbose) {
        message(sprintf(
            "  Rewrote %d GFF3 records' seqids (%d unmapped passed through).",
            n_renamed, n_unmapped
        ))
    }
    invisible(out_path)
}

# Persist the full sequence-report mapping (canonical, refseq, genbank,
# sequenceName, role, length) as a TSV at <groot>/chrom_aliases.tsv. Provides
# a self-describing aliases table for downstream tooling and future hooks
# into .compute_chrom_aliases.
.write_chrom_aliases_tsv <- function(groot, seqrep, rename_map) {
    df <- data.frame(
        canonical        = unname(rename_map[seqrep$refseqAccession]),
        refseqAccession  = seqrep$refseqAccession,
        genbankAccession = seqrep$genbankAccession,
        sequenceName     = seqrep$sequenceName,
        chrName          = seqrep$chrName,
        role             = seqrep$role,
        length           = seqrep$length,
        stringsAsFactors = FALSE
    )
    utils::write.table(
        df, file.path(groot, "chrom_aliases.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
    )
    invisible(NULL)
}

.run_ncbi_backend <- function(name, recipe, path, annotations, format, verbose) {
    accession <- recipe$accession
    chrom_naming <- recipe$chrom_naming %||% .NCBI_DEFAULT_CHROM_NAMING
    converter <- .gff3_to_genepred_resolve_or_install()

    if (verbose) {
        message(sprintf(
            "NCBI backend: accession=%s -> %s (chrom_naming=%s)",
            accession, path, chrom_naming
        ))
    }

    workdir <- tempfile("misha_ncbi_")
    dir.create(workdir, recursive = TRUE)
    on.exit(unlink(workdir, recursive = TRUE), add = TRUE)

    zip_path <- file.path(workdir, "datasets.zip")
    .download_to(.ncbi_datasets_zip_url(accession), zip_path, verbose = verbose)

    extract_dir <- file.path(workdir, "extract")
    dir.create(extract_dir)
    utils::unzip(zip_path, exdir = extract_dir)

    fasta_files <- list.files(extract_dir,
        pattern = "\\.(fna|fasta|fa)(\\.gz)?$",
        recursive = TRUE, full.names = TRUE
    )
    gff_files <- list.files(extract_dir,
        pattern = "\\.gff(\\.gz)?$",
        recursive = TRUE, full.names = TRUE
    )
    seqrep_files <- list.files(extract_dir,
        pattern = "sequence_report\\.jsonl$",
        recursive = TRUE, full.names = TRUE
    )

    if (!length(fasta_files)) {
        stop(sprintf("No FASTA file found in NCBI Datasets payload for %s", accession), call. = FALSE)
    }
    fasta_file <- fasta_files[[1]]
    if (grepl("\\.gz$", fasta_file)) {
        fasta_file <- .gunzip_to_file(fasta_file)
    }

    # Build the rename map from the sequence report. If absent (rare for
    # modern assemblies), fall back to "accession" naming with a warning.
    seqrep <- NULL
    rename_map <- NULL
    if (length(seqrep_files)) {
        seqrep <- .parse_ncbi_sequence_report(seqrep_files[[1]])
        rename_map <- .build_ncbi_rename_map(seqrep, chrom_naming)
        if (verbose) {
            message(sprintf(
                "  Sequence report: %d contigs (%d assembled, %d other)",
                nrow(seqrep),
                sum(seqrep$role == "assembled-molecule"),
                sum(seqrep$role != "assembled-molecule")
            ))
        }
    } else if (chrom_naming != "accession") {
        warning(sprintf(
            "No sequence_report.jsonl in NCBI payload for %s; falling back to chrom_naming='accession'.",
            accession
        ), call. = FALSE)
        chrom_naming <- "accession"
    }

    if (!is.null(rename_map) && chrom_naming != "accession") {
        renamed_fasta <- file.path(workdir, "renamed.fna")
        .rename_fasta_headers(fasta_file, renamed_fasta, rename_map, verbose = verbose)
        fasta_file <- renamed_fasta
    }

    genepred_file <- NULL
    if ("genes" %in% annotations && length(gff_files)) {
        gff_file <- gff_files[[1]]
        if (grepl("\\.gz$", gff_file)) {
            gff_file <- .gunzip_to_file(gff_file)
        }
        if (!is.null(rename_map) && chrom_naming != "accession") {
            renamed_gff <- file.path(workdir, "renamed.gff")
            .rename_gff3_seqids(gff_file, renamed_gff, rename_map, verbose = verbose)
            gff_file <- renamed_gff
        }
        raw_genepred <- file.path(workdir, "annot.raw.genePred")
        if (verbose) message(sprintf("  Running gff3ToGenePred on %s", basename(gff_file)))
        ret <- system2(converter,
            args = c(shQuote(gff_file), shQuote(raw_genepred)),
            stdout = if (verbose) "" else FALSE,
            stderr = if (verbose) "" else FALSE
        )
        if (ret != 0) {
            stop(sprintf("gff3ToGenePred failed (exit %d) on %s", ret, gff_file), call. = FALSE)
        }
        # gff3ToGenePred emits 15-col extended genePred — trim to the 12 cols
        # the C++ importer expects.
        genepred_file <- file.path(workdir, "annot.genePred")
        .normalize_ucsc_genepred(raw_genepred, genepred_file)
    } else if ("genes" %in% annotations) {
        warning(sprintf("No GFF3 found for %s; building seq-only genome (no genes/exons/utr*)", accession), call. = FALSE)
    }

    gdb.create(
        groot       = path,
        fasta       = fasta_file,
        genes.file  = genepred_file,
        annots.file = NULL,
        format      = format,
        verbose     = verbose
    )

    gdb.init(path, rescan = TRUE)

    if (!is.null(seqrep) && !is.null(rename_map)) {
        .write_chrom_aliases_tsv(path, seqrep, rename_map)
        if (verbose) message(sprintf("  Wrote chrom_aliases.tsv (%d rows).", nrow(seqrep)))
    }

    if ("cpgIsland" %in% annotations) {
        warning("cpgIsland not available from NCBI Datasets; skipping.", call. = FALSE)
    }
    if ("cytoband" %in% annotations) {
        warning("cytoband not available from NCBI Datasets; skipping.", call. = FALSE)
    }
    if ("rmsk" %in% annotations) {
        warning("rmsk loading from NCBI is not implemented in v1; skipping. Use the UCSC backend or the manual recipe to load repeats.", call. = FALSE)
    }

    .write_genome_info(
        path, name, recipe, annotations,
        list(
            fasta = list(name = basename(fasta_file)),
            genes = list(name = basename(genepred_file %||% ""))
        )
    )
    invisible(NULL)
}

# ---------------------------------------------------------------------------
# s3 / manual / local backends
# ---------------------------------------------------------------------------

.run_s3_backend <- function(name, recipe, path, annotations, format, verbose) {
    if (verbose) message(sprintf("S3 backend: aliasing to gdb.create_genome(%s)", recipe$assembly))
    parent_dir <- dirname(path)
    if (!dir.exists(parent_dir)) {
        dir.create(parent_dir, recursive = TRUE)
    }
    if (basename(path) != recipe$assembly) {
        warning(sprintf(
            "S3 backend extracts to <path>/%s; got path=%s. Final groot at %s/%s.",
            recipe$assembly, path, path, recipe$assembly
        ), call. = FALSE)
    }
    gdb.create_genome(recipe$assembly, path = parent_dir)
    invisible(NULL)
}

.run_local_backend <- function(name, recipe, path, annotations, format, verbose) {
    src <- recipe$path
    if (!dir.exists(src)) {
        stop(sprintf("Local groot does not exist: %s", src), call. = FALSE)
    }
    if (!file.exists(file.path(src, "chrom_sizes.txt"))) {
        stop(sprintf("Path %s does not look like a misha groot (no chrom_sizes.txt)", src), call. = FALSE)
    }
    if (verbose) message(sprintf("Local backend: gdb.init(%s)", src))
    gdb.init(src, rescan = TRUE)
    invisible(NULL)
}

.run_manual_backend <- function(name, recipe, path, annotations, format, verbose) {
    fasta <- recipe$fasta
    genes_file <- recipe$genes
    genes_format <- recipe$genes_format %||% "genepred"

    if (!is.null(genes_file) && genes_format == "gff3") {
        converter <- .gff3_to_genepred_resolve_or_install()
        workdir <- tempfile("misha_manual_")
        dir.create(workdir, recursive = TRUE)
        on.exit(unlink(workdir, recursive = TRUE), add = TRUE)

        local_gff <- file.path(workdir, "input.gff")
        .download_to(genes_file, local_gff, verbose = verbose)
        if (grepl("\\.gz$", genes_file)) {
            local_gff <- .gunzip_to_file(
                local_gff,
                file.path(workdir, "input.gff3")
            )
        }
        local_gp <- file.path(workdir, "input.genePred")
        ret <- system2(converter,
            args = c(shQuote(local_gff), shQuote(local_gp)),
            stdout = if (verbose) "" else FALSE,
            stderr = if (verbose) "" else FALSE
        )
        if (ret != 0) {
            stop(sprintf("gff3ToGenePred failed (exit %d)", ret), call. = FALSE)
        }
        genes_file <- local_gp
    } else if (!is.null(genes_file) && genes_format != "genepred") {
        stop(sprintf(
            "Unsupported genes_format '%s'. Supported: genepred, gff3.",
            genes_format
        ), call. = FALSE)
    }

    annots_file <- recipe$annots_file
    annots_names <- recipe$annots_names

    gdb.create(
        groot        = path,
        fasta        = fasta,
        genes.file   = genes_file,
        annots.file  = annots_file,
        annots.names = annots_names,
        format       = format,
        verbose      = verbose
    )

    gdb.init(path, rescan = TRUE)

    if ("rmsk" %in% annotations && !is.null(recipe$rmsk)) {
        tryCatch(
            {
                tmp <- tempfile(fileext = ".gz")
                on.exit(unlink(tmp), add = TRUE)
                .download_to(recipe$rmsk, tmp, verbose = verbose)
                unz <- .gunzip_to_file(tmp, sub("\\.gz$", "", tmp))
                .save_post_build_intervals("rmsk", .parse_ucsc_rmsk(unz), verbose = verbose)
            },
            error = function(e) warning(sprintf("rmsk load failed: %s", conditionMessage(e)), call. = FALSE)
        )
    }
    if ("cpgIsland" %in% annotations && !is.null(recipe$cpgIsland)) {
        tryCatch(
            {
                tmp <- tempfile(fileext = ".gz")
                on.exit(unlink(tmp), add = TRUE)
                .download_to(recipe$cpgIsland, tmp, verbose = verbose)
                unz <- .gunzip_to_file(tmp, sub("\\.gz$", "", tmp))
                .save_post_build_intervals("cpgIsland", .parse_ucsc_cpg_island(unz), verbose = verbose)
            },
            error = function(e) warning(sprintf("cpgIsland load failed: %s", conditionMessage(e)), call. = FALSE)
        )
    }
    if ("cytoband" %in% annotations && !is.null(recipe$cytoband)) {
        tryCatch(
            {
                tmp <- tempfile(fileext = ".gz")
                on.exit(unlink(tmp), add = TRUE)
                .download_to(recipe$cytoband, tmp, verbose = verbose)
                unz <- .gunzip_to_file(tmp, sub("\\.gz$", "", tmp))
                .save_post_build_intervals("cytoband", .parse_ucsc_cytoband(unz), verbose = verbose)
            },
            error = function(e) warning(sprintf("cytoband load failed: %s", conditionMessage(e)), call. = FALSE)
        )
    }

    .write_genome_info(path, name, recipe, annotations, list())
    invisible(NULL)
}

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

#' Build a misha genome database from a name
#'
#' Builds a misha genomic database for a named assembly by downloading FASTA
#' and (optionally) gene annotations, repeats, CpG islands, and cytobands from
#' UCSC golden path or NCBI Datasets, then dispatching to \code{\link{gdb.create}}.
#'
#' Compared to \code{\link{gdb.create_genome}} (which downloads a pre-baked
#' tarball from S3), \code{gdb.build_genome} builds from upstream sources at
#' call time and supports any UCSC golden-path assembly or any NCBI assembly
#' accession.
#'
#' \strong{Name resolution.} \code{name} is resolved through a chain:
#' \enumerate{
#'   \item the \code{registry} argument, if given;
#'   \item \code{getOption("misha.genome_registry")};
#'   \item a project-local \code{misha.yaml} (walked up to git root);
#'   \item the package's built-in \code{inst/genomes.yaml};
#'   \item if \code{name} matches \code{GC[FA]_<digits>.<digits>}, it is
#'     synthesized as \code{\{source: ncbi, accession: <name>\}}.
#' }
#'
#' \strong{Backends.}
#' \describe{
#'   \item{\code{ucsc}}{Golden path. Pulls
#'     \code{ncbiRefSeq.txt.gz} (genes, RefSeq), \code{ncbiRefSeqLink.txt.gz}
#'     (gene annotations), \code{rmsk.txt.gz}, \code{cpgIslandExt.txt.gz},
#'     \code{cytoBandIdeo.txt.gz} as requested.}
#'   \item{\code{ncbi}}{NCBI Datasets API. FASTA and GFF3 (RefSeq). Requires
#'     the \code{gff3ToGenePred} binary, downloaded on demand to
#'     \code{tools::R_user_dir("misha", "cache")} after user consent.
#'     \code{rmsk}/\code{cpgIsland}/\code{cytoband} are not loaded for NCBI
#'     assemblies in this version. Recipe field \code{chrom_naming} controls
#'     how the FASTA contigs are named on disk (see "Chromosome naming"
#'     below); the full sequence-report mapping is persisted at
#'     \code{<groot>/chrom_aliases.tsv}.}
#'   \item{\code{s3}}{Alias for \code{\link{gdb.create_genome}} (pre-baked
#'     tarball).}
#'   \item{\code{local}}{No download — \code{\link{gdb.init}} an existing
#'     groot at the registry-specified path.}
#'   \item{\code{manual}}{Use registry-supplied URLs verbatim.}
#' }
#'
#' On success a \code{genome_info.yaml} file is written at the groot root,
#' recording the source, recipe, and download time, and \code{\link{gdb.init}}
#' is called on the new groot.
#'
#' \strong{Chromosome naming (NCBI backend).} The recipe field
#' \code{chrom_naming} controls how the on-disk canonical contig names are
#' chosen. NCBI ships its FASTAs/GFFs keyed by RefSeq accession (e.g.
#' \code{NC_067374.1}); the sequence-report metadata maps these to friendlier
#' names. Allowed values:
#' \describe{
#'   \item{\code{"sequence_name"} (default)}{NCBI's \code{chrName} for
#'     assembled molecules (\code{1}, \code{2}, ..., \code{X}, \code{MT});
#'     RefSeq accession for unplaced/unlocalized scaffolds. Misha's existing
#'     \code{CHROM_ALIAS} mechanism then resolves \code{chr1 <-> 1} and
#'     \code{chrM <-> M <-> MT} automatically.}
#'   \item{\code{"ucsc"}}{UCSC-style names: \code{chr1}, \code{chr2}, ...,
#'     \code{chrX}, \code{chrM}, \code{chrUn_<acc>v<n>} for unplaced.}
#'   \item{\code{"accession"}}{Keep RefSeq accessions verbatim (no rename).}
#' }
#' If the NCBI payload lacks a sequence report, the build falls back to
#' \code{"accession"} with a warning. The full mapping (canonical, RefSeq,
#' GenBank, sequenceName, role) is written to
#' \code{<groot>/chrom_aliases.tsv} for downstream tooling.
#'
#' \strong{Alias loading.} On every \code{\link{gdb.init}}, if a
#' \code{chrom_aliases.tsv} is present at the groot, its
#' \code{refseqAccession} / \code{genbankAccession} / \code{sequenceName} /
#' \code{chrName} columns are added to \code{.misha$CHROM_ALIAS} as aliases
#' pointing at the canonical contig name. This makes
#' \code{gintervals("NC_067374.1", ...)},
#' \code{gintervals("CM028932.1", ...)}, and
#' \code{gintervals("contig_2989", ...)} all resolve to the same canonical
#' chromosome as \code{gintervals("chr1", ...)} or \code{gintervals("1", ...)}.
#'
#' @param name Genome name. Either a registry key (e.g. \code{"hg38"},
#'   \code{"UM_NZW_1.0"}) or an NCBI accession (e.g. \code{"GCF_009806435.1"}).
#' @param path Output directory for the new groot. Must not exist.
#' @param registry Optional path to a YAML registry file. Overrides the
#'   resolution chain.
#' @param annotations Character vector; subset of
#'   \code{c("genes","rmsk","cpgIsland","cytoband")}. Annotations not
#'   available from the chosen source emit a warning and are skipped.
#' @param format Database format: \code{"indexed"} (default) or
#'   \code{"per-chromosome"}. \code{NULL} uses
#'   \code{getOption("gmulticontig.indexed_format", TRUE)}.
#' @param verbose If \code{TRUE}, prints progress messages.
#' @return None; called for side effects (creates groot, then
#'   \code{gdb.init}s it).
#'
#' @seealso \code{\link{gdb.create}}, \code{\link{gdb.create_genome}},
#'   \code{\link{gdb.list_genomes}}, \code{\link{gdb.genome_info}},
#'   \code{\link{gdb.install_gff3_converter}}.
#'
#' @examples
#' \dontrun{
#' # Standard UCSC assembly (whitelisted in inst/genomes.yaml)
#' gdb.build_genome("hg38", path = "~/genomes/hg38")
#'
#' # Any NCBI assembly by accession
#' gdb.build_genome("GCF_009806435.1", path = "~/genomes/UM_NZW_1.0")
#'
#' # Sequence-only build (no gene annotation tracks)
#' gdb.build_genome("dm6", path = "/tmp/dm6", annotations = character(0))
#' }
#'
#' @export
gdb.build_genome <- function(name,
                             path = name,
                             registry = NULL,
                             annotations = c("genes", "rmsk", "cpgIsland", "cytoband"),
                             format = NULL,
                             verbose = TRUE) {
    if (file.exists(path)) {
        stop(sprintf("Output path '%s' already exists; refusing to overwrite. Choose a fresh path.", path), call. = FALSE)
    }
    annotations <- match.arg(annotations,
        choices = c("genes", "rmsk", "cpgIsland", "cytoband"),
        several.ok = TRUE
    )
    res <- .resolve_genome(name, registry = registry)
    recipe <- res$recipe
    if (verbose) {
        message(sprintf("Resolved '%s' from %s -> source=%s", name, res$resolved_from, recipe$source))
    }

    dispatcher <- switch(recipe$source,
        ucsc = .run_ucsc_backend,
        ncbi = .run_ncbi_backend,
        s3 = .run_s3_backend,
        local = .run_local_backend,
        manual = .run_manual_backend,
        stop(sprintf("Internal error: no backend for source '%s'", recipe$source), call. = FALSE)
    )
    dispatcher(name, recipe, path, annotations, format, verbose)
    invisible(NULL)
}


#' List resolvable genome names
#'
#' Returns a data frame describing every genome resolvable from the active
#' registry chain (see \code{\link{gdb.build_genome}} for the chain order).
#'
#' @param registry Optional path to an explicit registry YAML, overriding the
#'   resolution chain.
#' @return A data frame with columns:
#'   \itemize{
#'     \item \code{name} — registry key.
#'     \item \code{source} — backend (\code{ucsc}, \code{ncbi}, \code{s3},
#'       \code{local}, \code{manual}).
#'     \item \code{detail} — assembly / accession / path.
#'     \item \code{resolved_from} — which registry the entry came from.
#'   }
#'
#' @seealso \code{\link{gdb.build_genome}}, \code{\link{gdb.genome_info}}.
#'
#' @examples
#' gdb.list_genomes()
#'
#' @export
gdb.list_genomes <- function(registry = NULL) {
    sources <- list()
    if (!is.null(registry)) {
        sources[[length(sources) + 1]] <- list(label = sprintf("registry arg (%s)", registry), path = registry)
    }
    opt <- getOption("misha.genome_registry")
    if (!is.null(opt)) {
        sources[[length(sources) + 1]] <- list(label = sprintf("getOption (%s)", opt), path = opt)
    }
    proj <- .find_project_misha_yaml()
    if (!is.null(proj)) {
        sources[[length(sources) + 1]] <- list(label = sprintf("project (%s)", proj), path = proj)
    }
    builtin <- .builtin_registry_path()
    if (nzchar(builtin)) {
        sources[[length(sources) + 1]] <- list(label = "built-in", path = builtin)
    }

    rows <- list()
    seen <- character(0)
    for (src in sources) {
        entries <- .parse_genome_registry(src$path)
        for (nm in names(entries)) {
            if (nm %in% seen) next
            seen <- c(seen, nm)
            recipe <- tryCatch(.normalize_recipe(entries[[nm]], nm), error = function(e) NULL)
            if (is.null(recipe)) next
            detail <- recipe$assembly %||% recipe$accession %||% recipe$path %||%
                (if (!is.null(recipe$fasta)) paste(head(recipe$fasta, 1), collapse = ",") else NA_character_)
            rows[[length(rows) + 1]] <- data.frame(
                name = nm,
                source = recipe$source,
                detail = detail,
                resolved_from = src$label,
                stringsAsFactors = FALSE
            )
        }
    }
    if (!length(rows)) {
        return(data.frame(name = character(), source = character(), detail = character(), resolved_from = character(), stringsAsFactors = FALSE))
    }
    do.call(rbind, rows)
}


#' Inspect a resolved genome recipe without building
#'
#' Resolves \code{name} through the registry chain and returns the recipe (a
#' list) along with the source it was resolved from. Useful for previewing
#' what \code{\link{gdb.build_genome}} would do.
#'
#' @param name Genome name.
#' @param registry Optional path to an explicit registry YAML.
#' @return A list with components \code{recipe} (the resolved recipe) and
#'   \code{resolved_from} (the registry source).
#'
#' @seealso \code{\link{gdb.build_genome}}, \code{\link{gdb.list_genomes}}.
#'
#' @examples
#' gdb.genome_info("hg38")
#' gdb.genome_info("GCF_009806435.1")
#'
#' @export
gdb.genome_info <- function(name, registry = NULL) {
    .resolve_genome(name, registry = registry)
}


#' Pre-install UCSC's gff3ToGenePred binary
#'
#' Downloads UCSC's \code{gff3ToGenePred} static binary (~25 MB) into
#' \code{tools::R_user_dir("misha", "cache")/bin/}, verifies its SHA256, and
#' makes it executable. Used by \code{\link{gdb.build_genome}} when the
#' \code{ncbi} backend (or the \code{manual} backend with
#' \code{genes_format: gff3}) is invoked. Calling it directly is useful in CI
#' or in non-interactive scripts where the consent prompt would otherwise
#' fail.
#'
#' Override the binary location by setting environment variable
#' \code{MISHA_GFF3_TO_GENEPRED} to a binary you provide (for example, one
#' installed via \code{conda install -c bioconda ucsc-gff3togenepred}). This
#' is the recommended workaround on systems whose glibc is older than the one
#' UCSC's prebuilt binary requires.
#'
#' @param force If \code{TRUE}, skip the consent prompt and re-download even
#'   if the binary is already cached.
#' @return The cache path of the installed binary (invisibly).
#'
#' @examples
#' \dontrun{
#' gdb.install_gff3_converter()
#' Sys.setenv(MISHA_GFF3_TO_GENEPRED = "/path/to/your/gff3ToGenePred")
#' }
#'
#' @export
gdb.install_gff3_converter <- function(force = FALSE) {
    invisible(.install_gff3_converter(force = force))
}

#' Pre-install UCSC's gtfToGenePred binary
#'
#' Mirrors \code{\link{gdb.install_gff3_converter}}. Required for the
#' \code{ucsc-hub} backend's \code{genes} set (UCSC mammal hubs ship GTFs).
#'
#' Override the binary location by setting environment variable
#' \code{MISHA_GTF_TO_GENEPRED}.
#'
#' @param force If \code{TRUE}, skip consent prompt and re-download even if cached.
#' @return The cache path (invisibly).
#' @examples
#' \dontrun{
#' gdb.install_gtf_converter()
#' Sys.setenv(MISHA_GTF_TO_GENEPRED = "/path/to/your/gtfToGenePred")
#' }
#' @export
gdb.install_gtf_converter <- function(force = FALSE) {
    invisible(.install_gtf_converter(force = force))
}
