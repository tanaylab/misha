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

.MISHA_GENOME_SOURCES <- c("ucsc", "ucsc-hub", "ncbi", "s3", "manual", "local")

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
    } else if (src == "ucsc-hub") {
        if (is.null(recipe$accession)) miss("accession")
        if (!grepl("^GC[FA]_[0-9]+\\.[0-9]+$", recipe$accession)) {
            stop(sprintf(
                "Genome '%s': accession '%s' does not match GC[FA]_<digits>.<digits>",
                name, recipe$accession
            ), call. = FALSE)
        }
        if (!is.null(recipe$chrom_naming) &&
            !recipe$chrom_naming %in% c("ucsc", "accession", "sequence_name")) {
            stop(sprintf(
                "Genome '%s': chrom_naming '%s' invalid for ucsc-hub.",
                name, recipe$chrom_naming
            ), call. = FALSE)
        }
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

    # Pattern fallback for UCSC mammal hub accessions.
    if (grepl("^GC[FA]_[0-9]+\\.[0-9]+$", name)) {
        recipe <- list(source = "ucsc-hub", accession = name)
        return(list(
            recipe = recipe,
            resolved_from = "pattern fallback (UCSC mammal hub accession)"
        ))
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

# Write genome_info.yaml — a record of where this groot came from.
.write_genome_info <- function(groot, name, recipe, sets, files = list()) {
    info <- list(
        name = name,
        source = recipe$source,
        downloaded_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
        misha_version = as.character(utils::packageVersion("misha")),
        sets = as.list(sets),
        recipe = recipe,
        files = files
    )
    yaml::write_yaml(info, file.path(groot, "genome_info.yaml"))
    invisible(NULL)
}


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

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

#' Build a misha genome database from a name
#'
#' Builds a misha genomic database for a named assembly. Resolves the name
#' through the registry chain (or pattern-fallback for \code{GC[FA]_*}
#' accessions), downloads the FASTA, calls \code{\link{gdb.create}} to build
#' the seq-only groot, then dispatches to \code{\link{gdb.install_intervals}}
#' for the requested sets.
#'
#' For details on resolution, sources, sets, and chromosome-alias handling,
#' see \code{\link{gdb.install_intervals}}.
#'
#' @param name Genome name (registry key, alias, or \code{GC[FA]_*} accession).
#' @param path Output directory; must not exist.
#' @param registry Optional path to an explicit registry YAML.
#' @param sets Subset of \code{c("genes", "rmsk", "cgi", "cytoband")}.
#'   Empty vector \code{character(0)} = sequence-only build.
#' @param prefix Character scalar prepended to set names (see
#'   \code{\link{gdb.install_intervals}}).
#' @param gene_sets Named character vector mapping the four
#'   \code{gintervals.import_genes()} roles to on-disk set names; \code{NA} skips
#'   a role.
#' @param gtf_priority Character vector ordering GTF source preference.
#' @param format \code{"indexed"} or \code{"per-chromosome"}; \code{NULL} =>
#'   \code{getOption("gmulticontig.indexed_format", TRUE)}.
#' @param verbose If \code{TRUE}, prints progress.
#' @return None (invisible \code{NULL}).
#'
#' @seealso \code{\link{gdb.install_intervals}}, \code{\link{gdb.create}},
#'   \code{\link{gdb.list_genomes}}, \code{\link{gdb.genome_info}}.
#'
#' @examples
#' \dontrun{
#' gdb.build_genome("hg38", path = "~/genomes/hg38")
#' gdb.build_genome("GCA_004023825.1",
#'     path   = "~/genomes/arctic_fox",
#'     prefix = "intervs.global."
#' )
#' }
#'
#' @export
gdb.build_genome <- function(name,
                             path = name,
                             registry = NULL,
                             sets = c("genes", "rmsk", "cgi", "cytoband"),
                             prefix = "",
                             gene_sets = c(
                                 tss = "tss", exons = "exons",
                                 utr3 = "utr3", utr5 = "utr5"
                             ),
                             gtf_priority = c(
                                 "ncbiRefSeq", "bestRefSeq",
                                 "ensGene", "augustus", "xenoRefGene"
                             ),
                             format = NULL,
                             verbose = TRUE) {
    if (file.exists(path)) {
        stop(sprintf(
            "Output path '%s' already exists; refusing to overwrite. Choose a fresh path.",
            path
        ), call. = FALSE)
    }
    if (length(sets)) {
        sets <- match.arg(sets,
            choices = c("genes", "rmsk", "cgi", "cytoband"),
            several.ok = TRUE
        )
    }

    res <- .resolve_genome(name, registry = registry)
    recipe <- res$recipe
    if (verbose) {
        message(sprintf(
            "Resolved '%s' from %s -> source=%s",
            name, res$resolved_from, recipe$source
        ))
    }

    seq_info <- .build_seq(recipe, path, format = format, verbose = verbose)
    gdb.init(path, rescan = TRUE)

    if (length(sets)) {
        gdb.install_intervals(
            groot        = path,
            source       = recipe,
            sets         = sets,
            prefix       = prefix,
            gene_sets    = gene_sets,
            gtf_priority = gtf_priority,
            overwrite    = FALSE,
            registry     = NULL,
            verbose      = verbose
        )
    }

    .write_genome_info(path, name, recipe, sets, files = seq_info$files_record)
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

#' Install interval sets onto an existing groot
#'
#' Given an existing groot and a source recipe (or registry name, or accession),
#' fetches the relevant annotation files and installs interval sets — one or
#' more of \code{genes / rmsk / cgi / cytoband}.
#'
#' Decoupled from \code{\link{gdb.build_genome}} so that:
#' \itemize{
#'   \item users with a private FASTA build can layer canonical annotations onto it;
#'   \item failed installs can be resumed without re-fetching the FASTA;
#'   \item the same groot can host annotations from multiple sources under
#'         different prefixes (e.g. \code{intervs.global.}, \code{intervs.repeats.}).
#' }
#'
#' @param groot Path to a misha groot. \code{NULL} uses the active groot.
#' @param source Either a registry name, a recipe \code{list}, or a bare
#'   \code{GC[FA]_<digits>.<digits>} accession.
#' @param sets Subset of \code{c("genes", "rmsk", "cgi", "cytoband")}.
#' @param prefix Character scalar prepended verbatim to each set name. Include
#'   the trailing dot if you want one (e.g. \code{"intervs.global."}).
#' @param gene_sets Named character vector mapping
#'   \code{c("tss", "exons", "utr3", "utr5")} to the on-disk set name. \code{NA}
#'   value skips that role.
#' @param gtf_priority Character vector ordering GTF source preference for
#'   sources that ship multiple GTFs (currently \code{ucsc-hub}). First found wins.
#' @param overwrite If \code{FALSE} (default), error on existing target sets.
#'   If \code{TRUE}, remove existing sets before saving.
#' @param registry Optional path to a registry YAML; overrides the resolution chain.
#' @param verbose If \code{TRUE}, prints progress.
#' @return Invisible \code{NULL}. Side effects: writes \code{.interv} files under
#'   \code{<groot>/tracks/}, extends \code{<groot>/chrom_aliases.tsv}, appends to
#'   \code{<groot>/genome_info.yaml}, and re-initializes the active groot.
#'
#' @seealso \code{\link{gdb.build_genome}}, \code{\link{gdb.install_gtf_converter}}.
#'
#' @examples
#' \dontrun{
#' # Standalone install on an existing groot.
#' gdb.install_intervals(
#'     groot  = "/genomes/arctic_fox",
#'     source = "GCA_004023825.1",
#'     prefix = "intervs.global."
#' )
#'
#' # Layered: private FASTA groot + intervals from a UCSC hub assembly.
#' gdb.install_intervals(
#'     groot  = "/genomes/my_private",
#'     source = list(source = "ucsc-hub", accession = "GCF_009806435.1"),
#'     sets   = c("genes", "rmsk")
#' )
#' }
#' @export
gdb.install_intervals <- function(groot,
                                  source,
                                  sets = c("genes", "rmsk", "cgi", "cytoband"),
                                  prefix = "",
                                  gene_sets = c(
                                      tss = "tss", exons = "exons",
                                      utr3 = "utr3", utr5 = "utr5"
                                  ),
                                  gtf_priority = c(
                                      "ncbiRefSeq", "bestRefSeq",
                                      "ensGene", "augustus", "xenoRefGene"
                                  ),
                                  overwrite = FALSE,
                                  registry = NULL,
                                  verbose = TRUE) {
    sets <- match.arg(sets,
        choices = c("genes", "rmsk", "cgi", "cytoband"),
        several.ok = TRUE
    )
    if (!is.null(groot)) {
        if (!dir.exists(groot) ||
            !file.exists(file.path(groot, "chrom_sizes.txt"))) {
            stop(sprintf(
                "'%s' is not a misha groot (no chrom_sizes.txt). ",
                groot
            ), call. = FALSE)
        }
        gdb.init(groot, rescan = TRUE)
    } else if (!exists("GROOT", envir = .misha)) {
        stop("No active groot and no `groot` argument supplied.", call. = FALSE)
    }
    groot <- get("GROOT", envir = .misha)

    # Resolve source: list -> use as recipe; string -> registry/pattern.
    recipe <- if (is.list(source)) {
        .normalize_recipe(source, "<arg>")
        .validate_recipe(source, "<arg>")
        source
    } else if (is.character(source) && length(source) == 1L) {
        res <- .resolve_genome(source, registry = registry)
        if (verbose) {
            message(sprintf(
                "Resolved source '%s' from %s -> source=%s",
                source, res$resolved_from, res$recipe$source
            ))
        }
        res$recipe
    } else {
        stop("`source` must be a length-1 character or a recipe list.", call. = FALSE)
    }

    if (recipe$source %in% c("local", "s3")) {
        stop(sprintf("source '%s' has no fetchable intervals.", recipe$source),
            call. = FALSE
        )
    }

    workdir <- tempfile("misha_install_intervals_")
    dir.create(workdir, recursive = TRUE)
    on.exit(unlink(workdir, recursive = TRUE), add = TRUE)

    fetcher <- switch(recipe$source,
        ucsc = .ucsc_fetch_assets,
        `ucsc-hub` = function(r, s, w, v) .hub_fetch_assets(r, s, w, gtf_priority, v),
        ncbi = .ncbi_fetch_assets,
        manual = .manual_fetch_assets,
        stop(sprintf("No fetcher for source '%s'", recipe$source), call. = FALSE)
    )
    assets <- fetcher(recipe, sets, workdir, verbose)

    # chromAlias: detect groot column and source columns; build translator closure.
    groot_chroms <- as.character(get("ALLGENOME", envir = .misha)[[1]]$chrom)
    alias_df <- assets$chrom_alias$df
    groot_col <- if (!is.null(alias_df)) .detect_alias_column(alias_df, groot_chroms) else NA_character_
    if (!is.null(alias_df) && is.na(groot_col)) {
        scores <- attr(groot_col, "scores")
        stop(sprintf(
            "chromAlias has no column with 100%% coverage of groot chroms.\nPer-column overlap counts: %s\nFirst 5 unmapped groot chroms: %s",
            paste(sprintf("%s=%d/%d", names(scores), scores, length(groot_chroms)), collapse = ", "),
            paste(utils::head(setdiff(groot_chroms, unlist(alias_df, use.names = FALSE)), 5L),
                collapse = ", "
            )
        ), call. = FALSE)
    }
    if (!is.null(alias_df)) {
        .merge_chrom_aliases_tsv(groot, alias_df, groot_col)
        gdb.init(groot, rescan = TRUE) # reload CHROM_ALIAS
    }

    # Helper: per-asset translator. Detects the asset's own source column, asserts 100%.
    make_translator <- function(asset_chroms, asset_label) {
        if (is.null(alias_df)) {
            return(NULL)
        }
        src_col <- .detect_alias_column(alias_df, unique(asset_chroms))
        if (is.na(src_col)) {
            scores <- attr(src_col, "scores")
            stop(sprintf(
                "chromAlias has no column with 100%% coverage of distinct chroms in %s.\nPer-column overlap counts: %s\nFirst 5 unmapped chroms: %s",
                asset_label,
                paste(sprintf("%s=%d/%d", names(scores), scores, length(unique(asset_chroms))),
                    collapse = ", "
                ),
                paste(
                    utils::head(setdiff(
                        unique(asset_chroms),
                        unlist(alias_df, use.names = FALSE)
                    ), 5L),
                    collapse = ", "
                )
            ), call. = FALSE)
        }
        function(rows, chrom_col) {
            .translate_chroms(rows, chrom_col, alias_df, src_col, groot_col)
        }
    }

    # Genes.
    if ("genes" %in% sets && !is.null(assets$genes)) {
        # Sample first ~100 chroms from genePred/GTF/GFF for translator detection.
        chroms <- .sample_chroms_from_file(assets$genes$file, assets$genes$format)
        translator <- make_translator(chroms, "genes file")
        asset <- c(assets$genes, list(translate = translator))
        .install_genes_set(asset,
            prefix = prefix, gene_sets = gene_sets,
            overwrite = overwrite, verbose = verbose
        )
    }
    # rmsk.
    if ("rmsk" %in% sets && !is.null(assets$rmsk)) {
        df <- if (assets$rmsk$format == "rmsk-out") {
            .parse_rm_out(assets$rmsk$file, verbose)
        } else {
            .parse_ucsc_rmsk(assets$rmsk$file)
        }
        if (!is.null(alias_df)) {
            translator <- make_translator(unique(df$chrom), "rmsk")
            df <- translator(df, "chrom")
        }
        .install_rmsk_set(df, prefix = prefix, overwrite = overwrite, verbose = verbose)
    }
    # cgi.
    if ("cgi" %in% sets && !is.null(assets$cgi)) {
        df <- .parse_ucsc_cpg_island(assets$cgi$file)
        if (!is.null(alias_df)) {
            translator <- make_translator(unique(df$chrom), "cgi")
            df <- translator(df, "chrom")
        }
        .install_cgi_set(df, prefix = prefix, overwrite = overwrite, verbose = verbose)
    }
    # cytoband.
    if ("cytoband" %in% sets && !is.null(assets$cytoband)) {
        df <- .parse_ucsc_cytoband(assets$cytoband$file)
        if (!is.null(alias_df)) {
            translator <- make_translator(unique(df$chrom), "cytoband")
            df <- translator(df, "chrom")
        }
        .install_cytoband_set(df, prefix = prefix, overwrite = overwrite, verbose = verbose)
    }

    # Provenance: append to genome_info.yaml.
    .append_tracks_to_genome_info(groot, recipe, sets, prefix)

    # Final reload + summary.
    gdb.init(groot, rescan = TRUE)
    if (verbose) .install_intervals_summary(groot, recipe, sets, prefix)
    invisible(NULL)
}

# Sample distinct chroms from a genePred/GTF/GFF for translator detection.
# Reads up to 50,000 lines.
.sample_chroms_from_file <- function(file, format) {
    con <- if (grepl("\\.gz$", file)) gzfile(file, "rt") else file(file, "rt")
    on.exit(close(con), add = TRUE)
    chroms <- character(0)
    chunk_size <- 50000L
    repeat {
        lines <- readLines(con, n = chunk_size, warn = FALSE)
        if (!length(lines)) break
        # genePred: chrom in column 2; GTF/GFF: column 1.
        col_idx <- if (format == "genepred") 2L else 1L
        # Skip comment lines.
        lines <- lines[!startsWith(lines, "#") & nzchar(lines)]
        if (!length(lines)) next
        f <- strsplit(lines, "\t", fixed = TRUE)
        chroms <- unique(c(
            chroms,
            vapply(
                f, function(x) if (length(x) >= col_idx) x[[col_idx]] else NA_character_,
                character(1)
            )
        ))
        if (length(chroms) > 5000L) break # plenty for detection
    }
    chroms[!is.na(chroms) & nzchar(chroms)]
}

.install_intervals_summary <- function(groot, recipe, sets, prefix) {
    message("\ngdb.install_intervals: completed")
    message(sprintf("  groot:   %s", groot))
    message(sprintf(
        "  source:  %s%s", recipe$source,
        if (!is.null(recipe$accession)) {
            sprintf(" (%s)", recipe$accession)
        } else if (!is.null(recipe$assembly)) {
            sprintf(" (%s)", recipe$assembly)
        } else {
            ""
        }
    ))
    message(sprintf("  prefix:  %s", if (nzchar(prefix)) sprintf("\"%s\"", prefix) else "(none)"))
    message(sprintf("  sets:    %s", paste(sets, collapse = ", ")))
}

.append_tracks_to_genome_info <- function(groot, recipe, sets, prefix) {
    info_path <- file.path(groot, "genome_info.yaml")
    info <- if (file.exists(info_path)) yaml::read_yaml(info_path) else list()
    if (is.null(info$tracks)) info$tracks <- list()
    ts <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    for (s in sets) {
        info$tracks[[length(info$tracks) + 1L]] <- list(
            set          = paste0(prefix, s),
            source       = recipe$source,
            installed_at = ts
        )
    }
    yaml::write_yaml(info, info_path)
    invisible(NULL)
}
