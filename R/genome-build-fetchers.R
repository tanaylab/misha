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
