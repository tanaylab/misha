write_chain_entry <- function(con, srcName, srcSize, srcStrand, srcStart, srcEnd,
                              tgtName, tgtSize, tgtStrand, tgtStart, tgtEnd, id, score = 1000, append = NULL) {
    # "chain 1000 source1 100 + 10 30 chr1 44 + 0 20 1"
    # where "source*" is the source (query) side and "chr*" is the target side.
    # Use %.0f for large integers to avoid %d overflow

    # Auto-detect append mode: if file exists and has content, append
    if (is.null(append)) {
        append <- file.exists(con) && file.info(con)$size > 0
    }

    cat(
        sprintf(
            "chain %.0f %s %.0f %s %.0f %.0f %s %.0f %s %.0f %.0f %d\n",
            as.numeric(score),
            srcName, as.numeric(srcSize), srcStrand, as.numeric(srcStart), as.numeric(srcEnd),
            tgtName, as.numeric(tgtSize), tgtStrand, as.numeric(tgtStart), as.numeric(tgtEnd), id
        ),
        file = con, append = append
    )
    # single block size line (srcEnd - srcStart), then blank line
    cat(sprintf("%.0f\n\n", as.numeric(srcEnd - srcStart)), file = con, append = TRUE)
}

setup_db <- function(chrom_defs) {
    # chrom_defs = list(">chr1\nAAAA...", ... as strings to write)
    # OR list(">chr1\n", "AAAA...", "\n", ">chr2\n", "CCCC...", "\n") - header and sequence can be separate
    # .gseq.import() extracts chromosome name from filename: chr(\\w+)
    # We need separate FASTA files for each chromosome with names like chr1.fasta, chr2.fasta

    # Combine all elements into a single string first
    combined <- paste(chrom_defs, collapse = "")

    # Split by FASTA headers to get individual chromosomes
    # Handle any chromosome name (not just "chr" prefixed)
    lines <- strsplit(combined, "\n")[[1]]
    chrom_entries <- list()
    current_entry <- NULL

    for (line in lines) {
        if (grepl("^>", line)) {
            # Start a new chromosome entry (any header starting with ">")
            if (!is.null(current_entry)) {
                chrom_entries[[length(chrom_entries) + 1]] <- paste(current_entry, collapse = "\n")
            }
            current_entry <- c(line)
        } else if (!is.null(current_entry)) {
            # Add to current entry
            current_entry <- c(current_entry, line)
        }
    }
    # Add the last entry
    if (!is.null(current_entry)) {
        chrom_entries[[length(chrom_entries) + 1]] <- paste(current_entry, collapse = "\n")
    }

    if (length(chrom_entries) == 0) {
        stop("No valid FASTA entries found")
    }

    # Extract chromosome names and create files
    target_fastas <- character(length(chrom_entries))

    for (i in seq_along(chrom_entries)) {
        entry <- chrom_entries[[i]]
        # Extract header line (first line)
        entry_lines <- strsplit(entry, "\n")[[1]]
        header_line <- entry_lines[1]
        # Extract chromosome name from header (e.g., "chr1" from ">chr1" or "AncRef" from ">AncRef")
        chrom_name <- gsub("^>(\\S+).*$", "\\1", header_line, perl = TRUE)
        if (chrom_name == header_line || chrom_name == "") {
            stop(sprintf("Invalid FASTA header in entry %d: %s", i, header_line))
        }

        # .gseq.import() requires filename to start with "chr" and match pattern ^chr(\\w+)
        # The database chromosome name is extracted from filename: chr(\\w+) -> adds "chr" prefix
        # So if filename is "chrX.fasta", it extracts "X" and database has "chrX"
        # If filename is "chrAncRef.fasta", it extracts "AncRef" and database has "chrAncRef"
        # We need to update the FASTA header to match what will be in the database
        if (grepl("^chr", chrom_name)) {
            # Already starts with "chr", database will have same name
            filename_chrom <- chrom_name
            db_chrom_name <- chrom_name
        } else {
            # Doesn't start with "chr", filename will be "chr<name>.fasta"
            # Database will extract "<name>" and add "chr" -> "chr<name>"
            filename_chrom <- paste0("chr", chrom_name)
            db_chrom_name <- paste0("chr", chrom_name)
        }
        target_fastas[i] <- file.path(tempdir(), paste0(filename_chrom, ".fasta"))

        # Update FASTA header to match database chromosome name
        entry_lines[1] <- paste0(">", db_chrom_name)
        entry <- paste(entry_lines, collapse = "\n")

        # Ensure entry ends with newline for proper FASTA format
        if (!grepl("\n$", entry)) {
            entry <- paste0(entry, "\n")
        }
        cat(entry, file = target_fastas[i])
    }

    target_db <- tempfile()
    gdb.create(groot = target_db, fasta = target_fastas)
    gdb.init(target_db)
    withr::defer(
        {
            unlink(target_db, recursive = TRUE)
            unlink(target_fastas)
        },
        testthat::teardown_env()
    )
    target_fastas[1] # Return first file path for compatibility
}

new_chain_file <- function() {
    f <- tempfile(fileext = ".chain")
    withr::defer(unlink(f), testthat::teardown_env())
    f
}

has_liftover_binary <- function() {
    # Helper function to check if liftOver binary is available
    result <- tryCatch(
        {
            system2("which", "liftOver", stdout = TRUE, stderr = FALSE)
            TRUE
        },
        error = function(e) FALSE,
        warning = function(w) FALSE
    )
    if (!isTRUE(result)) {
        result <- tryCatch(
            {
                system2("liftOver", stdout = FALSE, stderr = FALSE)
                TRUE
            },
            error = function(e) FALSE,
            warning = function(w) FALSE
        )
    }
    return(isTRUE(result))
}

run_kent_liftover <- function(src_bed, chain_file, multiple = TRUE, minMatch = 1e-7, ret_cmd = FALSE) {
    withr::local_options(list(scipen = 999999))
    # Returns a list with 'mapped' (data.frame) and 'unmapped' (data.frame) results

    kent_result <- tempfile(fileext = ".bed")
    kent_unmapped <- tempfile(fileext = ".unmapped")

    # Run liftOver with -multiple flag if multiple is TRUE
    system2("liftOver", c(if (multiple) "-multiple" else NULL, if (!is.null(minMatch)) paste0("-minMatch=", minMatch) else NULL, src_bed, chain_file, kent_result, kent_unmapped),
        stdout = FALSE, stderr = FALSE
    )
    if (ret_cmd) {
        return(paste("liftOver ", if (multiple) "-multiple" else NULL, if (!is.null(minMatch)) paste0("-minMatch=", minMatch) else NULL, src_bed, chain_file, kent_result, kent_unmapped))
    }
    # Read mapped results (with name and value columns from -multiple)
    kent_data <- if (file.exists(kent_result) && file.info(kent_result)$size > 0) {
        if (multiple) {
            read.table(kent_result,
                header = FALSE,
                col.names = c("chrom", "start", "end", "name", "value"),
                stringsAsFactors = FALSE
            )
        } else {
            read.table(kent_result,
                header = FALSE,
                col.names = c("chrom", "start", "end"),
                stringsAsFactors = FALSE
            )
        }
    } else {
        data.frame(
            chrom = character(), start = integer(), end = integer(),
            name = character(), value = integer(), stringsAsFactors = FALSE
        )
    }

    # Read unmapped results
    kent_unmapped_data <- if (file.exists(kent_unmapped) && file.info(kent_unmapped)$size > 0) {
        if (multiple) {
            read.table(kent_unmapped,
                header = FALSE,
                col.names = c("chrom", "start", "end"),
                stringsAsFactors = FALSE
            )
        } else {
            read.table(kent_unmapped,
                header = FALSE,
                col.names = c("chrom", "start", "end"),
                stringsAsFactors = FALSE
            )
        }
    } else {
        data.frame(
            chrom = character(), start = integer(), end = integer(),
            stringsAsFactors = FALSE
        )
    }

    # Clean up temp files
    unlink(c(kent_result, kent_unmapped))

    return(list(mapped = kent_data, unmapped = kent_unmapped_data))
}
