# R/genome-build-fetchers.R
# Per-source asset fetchers. Return a uniform list consumed by per-set installers.

# Internal: parse an Apache directory-index HTML body via regex.
# Exposed for testing without network.
.hub_list_dir_parse <- function(body) {
    hits <- regmatches(body, gregexpr('href="([^"?/][^"]*)"', body, perl = TRUE))[[1L]]
    hits <- sub('^href="', "", hits)
    hits <- sub('"$', "", hits)
    hits[!hits %in% c("/", "../") & !startsWith(hits, "?")]
}

# Fetch a UCSC hub-style directory listing. Returns a character vector of
# filenames (or NULL if URL is unreachable or 404). Uses curl (already in
# Imports) — no xml2/rvest.
.hub_list_dir <- function(url, verbose = TRUE) {
    if (verbose) message(sprintf("Listing %s ...", url))
    h <- curl::new_handle()
    resp <- tryCatch(curl::curl_fetch_memory(url, handle = h), error = function(e) NULL)
    if (is.null(resp) || resp$status_code >= 400) {
        return(NULL)
    }
    body <- rawToChar(resp$content)
    .hub_list_dir_parse(body)
}

# URL for a UCSC mammal-hub directory by GenBank/RefSeq accession.
.hub_url_for <- function(accession) {
    m <- regmatches(accession, regexec(
        "^(GC[FA])_([0-9]{3})([0-9]{3})([0-9]{3})\\.[0-9]+$",
        accession
    ))[[1L]]
    if (length(m) != 5L) {
        stop(sprintf("Invalid accession '%s'; expected GC[FA]_<9 digits>.<n>", accession),
            call. = FALSE
        )
    }
    sprintf(
        "https://hgdownload.soe.ucsc.edu/hubs/%s/%s/%s/%s/%s/",
        m[[2L]], m[[3L]], m[[4L]], m[[5L]], accession
    )
}

# From a list of filenames in a hub's genes/ subdir, pick the first matching
# entry of `priority`. Returns list(file = "<path>", source = "<which>") or NULL.
.pick_gtf <- function(files, priority) {
    for (src in priority) {
        pat <- sprintf("\\.%s\\.gtf(\\.gz)?$", src)
        hit <- files[grepl(pat, files)]
        if (length(hit)) {
            return(list(file = hit[[1L]], source = src))
        }
    }
    NULL
}

# Fetch every asset needed for the requested sets from a UCSC mammal hub.
# Downloads files into a workdir (caller manages cleanup via on.exit).
# Returns a list with named entries per asset, each:
#   list(file = <local path>, format = <tag>, ...)
# plus chrom_alias = list(file=..., df=<parsed data.frame>).
.hub_fetch_assets <- function(recipe, sets, workdir, gtf_priority,
                              verbose = TRUE) {
    base <- .hub_url_for(recipe$accession)
    files <- .hub_list_dir(base, verbose = verbose)
    if (is.null(files)) {
        stop(sprintf(
            "UCSC mammal hub directory not found: %s\nTry source: ncbi for this accession.",
            base
        ), call. = FALSE)
    }
    out <- list()

    # chromAlias is always fetched (we need it even for non-translation cases).
    alias_file <- files[grepl("\\.chromAlias\\.txt$", files)]
    if (!length(alias_file)) {
        out$chrom_alias <- NULL
        if (verbose) message("  No chromAlias.txt at hub; alias-driven translation disabled.")
    } else {
        local <- file.path(workdir, basename(alias_file[[1L]]))
        .download_to(paste0(base, alias_file[[1L]]), local, verbose = verbose)
        out$chrom_alias <- list(file = local, df = .parse_ucsc_chromalias(local))
    }

    if ("rmsk" %in% sets) {
        rmsk_file <- files[grepl("\\.repeatMasker\\.out\\.gz$", files)]
        if (length(rmsk_file)) {
            local <- file.path(workdir, basename(rmsk_file[[1L]]))
            .download_to(paste0(base, rmsk_file[[1L]]), local, verbose = verbose)
            out$rmsk <- list(file = local, format = "rmsk-out")
        } else {
            warning(sprintf(
                "'rmsk' requested but no repeatMasker.out.gz at %s; skipping.",
                base
            ), call. = FALSE)
        }
    }

    if ("genes" %in% sets) {
        genes_files <- .hub_list_dir(paste0(base, "genes/"), verbose = verbose)
        if (is.null(genes_files)) {
            warning(sprintf("'genes' requested but no genes/ subdir at %s; skipping.", base),
                call. = FALSE
            )
        } else {
            picked <- .pick_gtf(genes_files, priority = gtf_priority)
            if (is.null(picked)) {
                warning(sprintf(
                    "'genes' requested but none of %s found at %s/genes/; skipping.",
                    paste(gtf_priority, collapse = ", "), base
                ), call. = FALSE)
            } else {
                local <- file.path(workdir, basename(picked$file))
                .download_to(paste0(base, "genes/", picked$file), local, verbose = verbose)
                out$genes <- list(file = local, format = "gtf", gtf_source = picked$source)
            }
        }
    }

    if ("cgi" %in% sets) {
        cgi_file <- files[grepl("\\.cpgIslandExt\\.txt\\.gz$", files)]
        if (length(cgi_file)) {
            local <- file.path(workdir, basename(cgi_file[[1L]]))
            .download_to(paste0(base, cgi_file[[1L]]), local, verbose = verbose)
            out$cgi <- list(file = local, format = "ucsc-cpg-11col")
        } else {
            warning(sprintf(
                "'cgi' requested but no cpgIslandExt.txt.gz at %s; skipping.",
                base
            ), call. = FALSE)
        }
    }

    if ("cytoband" %in% sets) {
        warning("'cytoband' is not available from UCSC mammal hubs; skipping.", call. = FALSE)
    }

    out
}
