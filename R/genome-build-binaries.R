# R/genome-build-binaries.R
# UCSC converter binaries (gff3ToGenePred, gtfToGenePred) — install on demand,
# SHA-pinned, cached under tools::R_user_dir("misha", "cache").

# UCSC gff3ToGenePred binary URLs (no version tag — UCSC overwrites in place).
# SHAs pinned at implementation time (2026-04-30); a UCSC rebuild will trigger
# an integrity error pointing users at MISHA_GFF3_TO_GENEPRED.
.GFF3_TO_GENEPRED_URLS <- list(
    linux_x86_64 = "https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gff3ToGenePred",
    macos_x86_64 = "https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/gff3ToGenePred",
    macos_arm64  = "https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/gff3ToGenePred"
)

.GFF3_TO_GENEPRED_SHA256 <- list(
    linux_x86_64 = "8a62fde0c395078ae173f36762bd554804ec1e86b9d737292dd504b05570e678",
    macos_x86_64 = "db1a34d01f00fa21cafc995d16b065b860a89650b4a9804a6bdf6585c0224cb4",
    macos_arm64  = "9896c1566f06b9779ff7031ae00e6abd18f00111eafc2774510e13dece25c9df"
)

.GFF3_TO_GENEPRED_ENV <- "MISHA_GFF3_TO_GENEPRED"

# ---------------------------------------------------------------------------
# Architecture / platform detection
# ---------------------------------------------------------------------------

.detect_arch <- function() {
    sys <- tolower(Sys.info()[["sysname"]])
    machine <- tolower(Sys.info()[["machine"]])

    if (sys == "linux" && machine %in% c("x86_64", "amd64")) {
        return("linux_x86_64")
    }
    if (sys == "darwin" && machine %in% c("x86_64", "amd64")) {
        return("macos_x86_64")
    }
    if (sys == "darwin" && machine %in% c("arm64", "aarch64")) {
        return("macos_arm64")
    }
    stop(sprintf(
        "No prebuilt gff3ToGenePred binary for %s/%s. Set %s to a binary you provide (e.g. via conda: ucsc-gff3togenepred).",
        sys, machine, .GFF3_TO_GENEPRED_ENV
    ), call. = FALSE)
}

# ---------------------------------------------------------------------------
# gff3ToGenePred binary management
# ---------------------------------------------------------------------------

.gff3_to_genepred_cache_dir <- function() {
    file.path(tools::R_user_dir("misha", "cache"), "bin")
}

.gff3_to_genepred_cache_path <- function() {
    file.path(.gff3_to_genepred_cache_dir(), "gff3ToGenePred")
}

# Returns absolute path to a usable gff3ToGenePred binary, or NULL if none
# is available. Never downloads — caller decides whether to install.
.gff3_to_genepred_path <- function() {
    env <- Sys.getenv(.GFF3_TO_GENEPRED_ENV, unset = "")
    if (nzchar(env)) {
        if (!file.exists(env)) {
            stop(sprintf("%s points at non-existent file: %s", .GFF3_TO_GENEPRED_ENV, env), call. = FALSE)
        }
        return(normalizePath(env, mustWork = TRUE))
    }
    cached <- .gff3_to_genepred_cache_path()
    if (!file.exists(cached)) {
        return(NULL)
    }
    arch <- .detect_arch()
    expected_sha <- .GFF3_TO_GENEPRED_SHA256[[arch]]
    actual_sha <- digest::digest(file = cached, algo = "sha256")
    if (!identical(actual_sha, expected_sha)) {
        return(NULL) # Caller will re-prompt.
    }
    cached
}

.prompt_yes_no <- function(message, default_yes = TRUE) {
    if (!interactive()) {
        return(FALSE)
    }
    suffix <- if (default_yes) "[Y/n]" else "[y/N]"
    ans <- readline(prompt = paste0(message, " ", suffix, ": "))
    ans <- tolower(trimws(ans))
    if (!nzchar(ans)) {
        return(default_yes)
    }
    ans %in% c("y", "yes")
}

# Install (download + verify + cache) the binary. Returns the cache path.
# Errors if user declines or session is non-interactive.
.install_gff3_converter <- function(force = FALSE) {
    arch <- .detect_arch()
    url <- .GFF3_TO_GENEPRED_URLS[[arch]]
    expected_sha <- .GFF3_TO_GENEPRED_SHA256[[arch]]
    cache_dir <- .gff3_to_genepred_cache_dir()
    cache_path <- .gff3_to_genepred_cache_path()

    if (!force) {
        message(sprintf(
            "misha needs UCSC's gff3ToGenePred binary (%s) to build genomes from NCBI sources.",
            arch
        ))
        message(sprintf("  source: %s", url))
        message(sprintf("  cache:  %s  (~25 MB)", cache_path))
        if (!.prompt_yes_no("Download now?")) {
            stop(sprintf(
                "Cannot proceed without gff3ToGenePred. Set %s to a binary you provide, or call gdb.install_gff3_converter(force = TRUE).",
                .GFF3_TO_GENEPRED_ENV
            ), call. = FALSE)
        }
    }

    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE, mode = "0755")
    tmp <- tempfile(tmpdir = cache_dir, fileext = ".part")
    on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)

    utils::download.file(url, tmp, mode = "wb", quiet = FALSE)

    actual_sha <- digest::digest(file = tmp, algo = "sha256")
    if (!identical(actual_sha, expected_sha)) {
        stop(sprintf(
            "gff3ToGenePred binary integrity check failed.\n  expected: %s\n  got:      %s\nUCSC may have published a new build. Set %s to override, or update misha.",
            expected_sha, actual_sha, .GFF3_TO_GENEPRED_ENV
        ), call. = FALSE)
    }

    file.rename(tmp, cache_path)
    Sys.chmod(cache_path, mode = "0755")
    message(sprintf("Installed gff3ToGenePred to %s", cache_path))
    cache_path
}

# Resolve a path, installing if needed (with consent).
.gff3_to_genepred_resolve_or_install <- function(force_install = FALSE) {
    p <- .gff3_to_genepred_path()
    if (!is.null(p) && !force_install) {
        return(p)
    }
    .install_gff3_converter(force = force_install)
}

# ---------------------------------------------------------------------------
# gtfToGenePred binary management
# ---------------------------------------------------------------------------

# UCSC gtfToGenePred binary URLs (no version tag — UCSC overwrites in place).
# SHAs pinned at implementation time (2026-05-03); a UCSC rebuild will trigger
# an integrity error pointing users at MISHA_GTF_TO_GENEPRED.
.GTF_TO_GENEPRED_URLS <- list(
    linux_x86_64 = "https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred",
    macos_x86_64 = "https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/gtfToGenePred",
    macos_arm64  = "https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/gtfToGenePred"
)

.GTF_TO_GENEPRED_SHA256 <- list(
    linux_x86_64 = "17932c3b0deb755e961a3f2cfbfef995113f770b334ff90f6359e860e2d3f949",
    macos_x86_64 = "eff2fb6e0bae40f0453ce79b1d751e32fc43cfa3632c15e89ca2735a2788ccb3",
    macos_arm64  = "d0da218eef522c8de3c7f05f98ed0946e3eb864413e20c9415d8ad139bb1a478"
)

.GTF_TO_GENEPRED_ENV <- "MISHA_GTF_TO_GENEPRED"

.gtf_to_genepred_cache_path <- function() {
    file.path(.gff3_to_genepred_cache_dir(), "gtfToGenePred")
}

# Returns absolute path to a usable gtfToGenePred binary, or NULL if none
# is available. Never downloads — caller decides whether to install.
.gtf_to_genepred_path <- function() {
    env <- Sys.getenv(.GTF_TO_GENEPRED_ENV, unset = "")
    if (nzchar(env)) {
        if (!file.exists(env)) {
            stop(sprintf("%s points at non-existent file: %s", .GTF_TO_GENEPRED_ENV, env), call. = FALSE)
        }
        return(normalizePath(env, mustWork = TRUE))
    }
    cached <- .gtf_to_genepred_cache_path()
    if (!file.exists(cached)) {
        return(NULL)
    }
    arch <- .detect_arch()
    expected_sha <- .GTF_TO_GENEPRED_SHA256[[arch]]
    actual_sha <- digest::digest(file = cached, algo = "sha256")
    if (!identical(actual_sha, expected_sha)) {
        return(NULL) # Caller will re-prompt.
    }
    cached
}

# Install (download + verify + cache) the binary. Returns the cache path.
# Errors if user declines or session is non-interactive.
.install_gtf_converter <- function(force = FALSE) {
    arch <- .detect_arch()
    url <- .GTF_TO_GENEPRED_URLS[[arch]]
    expected_sha <- .GTF_TO_GENEPRED_SHA256[[arch]]
    cache_dir <- .gff3_to_genepred_cache_dir()
    cache_path <- .gtf_to_genepred_cache_path()

    if (!force) {
        message(sprintf("misha needs UCSC's gtfToGenePred binary (%s) to import hub GTFs.", arch))
        message(sprintf("  source: %s", url))
        message(sprintf("  cache:  %s", cache_path))
        if (!.prompt_yes_no("Download now?")) {
            stop(sprintf(
                "Cannot proceed without gtfToGenePred. Set %s to a binary you provide, or call gdb.install_gtf_converter(force = TRUE).",
                .GTF_TO_GENEPRED_ENV
            ), call. = FALSE)
        }
    }

    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE, mode = "0755")
    tmp <- tempfile(tmpdir = cache_dir, fileext = ".part")
    on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
    utils::download.file(url, tmp, mode = "wb", quiet = FALSE)

    actual_sha <- digest::digest(file = tmp, algo = "sha256")
    if (!identical(actual_sha, expected_sha)) {
        stop(sprintf(
            "gtfToGenePred binary integrity check failed.\n  expected: %s\n  got:      %s\nUCSC may have published a new build. Set %s to override, or update misha.",
            expected_sha, actual_sha, .GTF_TO_GENEPRED_ENV
        ), call. = FALSE)
    }

    file.rename(tmp, cache_path)
    Sys.chmod(cache_path, mode = "0755")
    message(sprintf("Installed gtfToGenePred to %s", cache_path))
    cache_path
}

# Resolve a path, installing if needed (with consent).
.gtf_to_genepred_resolve_or_install <- function(force_install = FALSE) {
    p <- .gtf_to_genepred_path()
    if (!is.null(p) && !force_install) {
        return(p)
    }
    .install_gtf_converter(force = force_install)
}
