# Tests for .gsynth_process_parallel helper function
# This helper extracts the common parallel processing pattern from gsynth functions

# Re-roots into a database of its own; hand this process's overlay back to the next
# file in the parallel worker (helper-test_db.R).
restore_groot_on_exit()


test_that(".gsynth_process_parallel returns NULL for small genomes", {
    gdb.init_examples()

    # Small intervals should not trigger parallel processing
    small_intervals <- gintervals(1, 0, 1000)

    # Simple chunk processor that returns a temp file
    process_fn <- function(chunk_interval, chunk_index, ...) {
        temp_file <- tempfile(fileext = ".txt")
        writeLines(sprintf("chunk_%d", chunk_index), temp_file)
        temp_file
    }

    result <- .gsynth_process_parallel(
        intervals = small_intervals,
        output_format = "fasta",
        output_path = tempfile(fileext = ".fa"),
        process_chunk_fn = process_fn,
        max_chunk_size = 1e9
    )

    # Should return NULL indicating no parallel processing was needed
    expect_null(result)
})

test_that(".gsynth_process_parallel processes file output correctly", {
    gdb.init_examples()

    # Create intervals large enough to trigger parallel processing
    # by using a small max_chunk_size
    intervals <- gintervals(c(1, 1, 1), c(0, 1000, 2000), c(1000, 2000, 3000))

    output_file <- tempfile(fileext = ".txt")

    # Chunk processor that writes content to a temp file
    process_fn <- function(chunk_interval, chunk_index, ...) {
        temp_file <- tempfile(fileext = ".txt")
        content <- sprintf(
            "chunk_%d:chrom=%s,start=%d,end=%d",
            chunk_index,
            chunk_interval$chrom,
            chunk_interval$start,
            chunk_interval$end
        )
        writeLines(content, temp_file)
        temp_file
    }

    # Use a small max_chunk_size to force parallel processing
    result <- .gsynth_process_parallel(
        intervals = intervals,
        output_format = "fasta", # text mode
        output_path = output_file,
        process_chunk_fn = process_fn,
        max_chunk_size = 100 # Force chunking
    )

    # Result should be invisible(NULL) for file output
    expect_null(result)

    # Output file should exist and contain all chunks
    expect_true(file.exists(output_file))
    content <- readLines(output_file)
    expect_equal(length(content), 3)
    expect_true(grepl("chunk_1", content[1]))
    expect_true(grepl("chunk_2", content[2]))
    expect_true(grepl("chunk_3", content[3]))

    unlink(output_file)
})

test_that(".gsynth_process_parallel processes vector output correctly", {
    gdb.init_examples()

    # Create intervals
    intervals <- gintervals(c(1, 1, 1), c(0, 1000, 2000), c(1000, 2000, 3000))

    # Chunk processor that returns a vector
    process_fn <- function(chunk_interval, chunk_index, ...) {
        c(sprintf("seq_%d_a", chunk_index), sprintf("seq_%d_b", chunk_index))
    }

    # Use a small max_chunk_size to force parallel processing
    result <- .gsynth_process_parallel(
        intervals = intervals,
        output_format = "vector",
        output_path = NULL,
        process_chunk_fn = process_fn,
        max_chunk_size = 100 # Force chunking
    )

    # Result should be unlisted vector
    expect_true(is.character(result))
    expect_equal(length(result), 6) # 3 chunks * 2 sequences each
    expect_equal(result[1], "seq_1_a")
    expect_equal(result[2], "seq_1_b")
    expect_equal(result[3], "seq_2_a")
})

test_that(".gsynth_process_parallel handles binary file output", {
    gdb.init_examples()

    intervals <- gintervals(c(1, 1), c(0, 1000), c(1000, 2000))
    output_file <- tempfile(fileext = ".bin")

    # Chunk processor that writes binary content
    process_fn <- function(chunk_interval, chunk_index, ...) {
        temp_file <- tempfile(fileext = ".bin")
        writeBin(as.raw(chunk_index:(chunk_index + 4)), temp_file)
        temp_file
    }

    result <- .gsynth_process_parallel(
        intervals = intervals,
        output_format = "misha", # binary mode
        output_path = output_file,
        process_chunk_fn = process_fn,
        max_chunk_size = 100
    )

    expect_null(result)
    expect_true(file.exists(output_file))

    # Read back and verify binary content was concatenated
    content <- readBin(output_file, "raw", n = 100)
    expect_equal(length(content), 10) # 5 bytes from each chunk

    unlink(output_file)
})

test_that(".gsynth_process_parallel generates reproducible seeds", {
    gdb.init_examples()

    intervals <- gintervals(c(1, 1), c(0, 1000), c(1000, 2000))

    # Chunk processor that uses the seed to generate a random value
    # This tests that seeds are passed and produce reproducible results
    process_fn <- function(chunk_interval, chunk_index, chunk_seed = NULL, ...) {
        if (!is.null(chunk_seed)) {
            set.seed(chunk_seed)
            random_val <- runif(1)
        } else {
            random_val <- NA
        }
        sprintf("chunk_%d_seed_%.6f", chunk_index, random_val)
    }

    # Run twice with same seed - results should be identical
    result1 <- .gsynth_process_parallel(
        intervals = intervals,
        output_format = "vector",
        output_path = NULL,
        process_chunk_fn = process_fn,
        max_chunk_size = 100,
        seed = 42
    )

    result2 <- .gsynth_process_parallel(
        intervals = intervals,
        output_format = "vector",
        output_path = NULL,
        process_chunk_fn = process_fn,
        max_chunk_size = 100,
        seed = 42
    )

    # Results should be deterministic given the same input seed
    expect_equal(result1, result2)
    expect_equal(length(result1), 2) # One result per chunk

    # Verify seeds were actually used (not NA)
    expect_false(any(grepl("NA", result1)))

    # Run with different seed - results should differ
    result3 <- .gsynth_process_parallel(
        intervals = intervals,
        output_format = "vector",
        output_path = NULL,
        process_chunk_fn = process_fn,
        max_chunk_size = 100,
        seed = 99
    )

    expect_false(identical(result1, result3))
})

test_that(".gsynth_process_parallel passes extra arguments to process function", {
    gdb.init_examples()

    intervals <- gintervals(c(1, 1), c(0, 1000), c(1000, 2000))

    # Process function that includes extra arguments in output
    # This tests that extra args are actually passed through
    process_fn <- function(chunk_interval, chunk_index, chunk_seed = NULL,
                           extra_arg1, extra_arg2, ...) {
        sprintf("chunk_%d_arg1=%s_arg2=%d", chunk_index, extra_arg1, extra_arg2)
    }

    result <- .gsynth_process_parallel(
        intervals = intervals,
        output_format = "vector",
        output_path = NULL,
        process_chunk_fn = process_fn,
        max_chunk_size = 100,
        extra_arg1 = "test_value",
        extra_arg2 = 123
    )

    expect_equal(length(result), 2)
    # Verify extra args were passed by checking the output contains them
    expect_true(grepl("arg1=test_value", result[1]))
    expect_true(grepl("arg2=123", result[1]))
    expect_true(grepl("arg1=test_value", result[2]))
    expect_true(grepl("arg2=123", result[2]))
})

test_that(".gsynth_process_parallel handles failed chunks gracefully", {
    gdb.init_examples()

    intervals <- gintervals(c(1, 1), c(0, 1000), c(1000, 2000))
    output_file <- tempfile(fileext = ".txt")

    # Chunk processor that fails on second chunk
    process_fn <- function(chunk_interval, chunk_index, ...) {
        if (chunk_index == 2) {
            return("/nonexistent/file/path.txt") # Return non-existent file
        }
        temp_file <- tempfile(fileext = ".txt")
        writeLines(sprintf("chunk_%d", chunk_index), temp_file)
        temp_file
    }

    expect_error(
        .gsynth_process_parallel(
            intervals = intervals,
            output_format = "fasta",
            output_path = output_file,
            process_chunk_fn = process_fn,
            max_chunk_size = 100
        ),
        "Failed to generate chunk"
    )

    unlink(output_file)
})

test_that(".gsynth_process_parallel respects gmax.processes option", {
    gdb.init_examples()

    intervals <- gintervals(
        c(1, 1, 1, 1, 1),
        c(0, 1000, 2000, 3000, 4000),
        c(1000, 2000, 3000, 4000, 5000)
    )

    old_option <- getOption("gmax.processes")
    options(gmax.processes = 2)

    # We can't easily test the actual number of cores used,
    # but we can verify the function respects the option without error
    process_fn <- function(chunk_interval, chunk_index, ...) {
        sprintf("chunk_%d", chunk_index)
    }

    result <- .gsynth_process_parallel(
        intervals = intervals,
        output_format = "vector",
        output_path = NULL,
        process_chunk_fn = process_fn,
        max_chunk_size = 100
    )

    expect_equal(length(result), 5)

    options(gmax.processes = old_option)
})

# Lowers the size above which gsynth.random and gsynth.replace_kmer process the intervals in
# parallel chunks (1e9 bp), so that small intervals reach that path.
local_gsynth_chunk_size <- function(size, env = parent.frame()) {
    ns <- asNamespace("misha")
    old <- get(".GSYNTH_MAX_CHUNK_SIZE", envir = ns)
    unlockBinding(".GSYNTH_MAX_CHUNK_SIZE", ns)
    assign(".GSYNTH_MAX_CHUNK_SIZE", size, envir = ns)
    withr::defer(
        {
            assign(".GSYNTH_MAX_CHUNK_SIZE", old, envir = ns)
            lockBinding(".GSYNTH_MAX_CHUNK_SIZE", ns)
        },
        envir = env
    )
}

test_that("a row larger than the chunk size is processed in one call", {
    gdb.init_examples()
    iv <- gintervals(1, 0, 20000)
    expected_random <- gsynth.random(intervals = iv, output_format = "vector", seed = 3)
    expected_replace <- gsynth.replace_kmer("CG", "GC", intervals = iv, output_format = "vector")

    local_gsynth_chunk_size(1e4)
    # each used to call itself on the same row until the C stack overflowed
    expect_identical(gsynth.random(intervals = iv, output_format = "vector", seed = 3), expected_random)
    expect_identical(gsynth.replace_kmer("CG", "GC", intervals = iv, output_format = "vector"), expected_replace)

    # rows larger than the chunk size among others
    ivs <- gintervals(c(1, 1, 2), c(0, 20000, 0), c(20000, 25000, 15000))
    expect_equal(nchar(gsynth.random(intervals = ivs, output_format = "vector", seed = 3)), c(20000, 5000, 15000))
})

test_that("parallel FASTA output gets a .fai and leaves no chunk files behind", {
    gdb.init_examples()
    ivs <- gintervals(c(1, 1, 2), c(0, 9000, 0), c(9000, 18000, 5000))
    serial <- c(random = tempfile(fileext = ".fa"), replace = tempfile(fileext = ".fa"))
    parallel <- c(random = tempfile(fileext = ".fa"), replace = tempfile(fileext = ".fa"))
    withr::defer(unlink(c(serial, parallel, paste0(c(serial, parallel), ".fai"))))
    gsynth.random(intervals = ivs, output_path = serial[["random"]], output_format = "fasta", seed = 3)
    gsynth.replace_kmer("CG", "GC", intervals = ivs, output_path = serial[["replace"]], output_format = "fasta")

    local_gsynth_chunk_size(1e4)
    before <- list.files(tempdir())
    gsynth.random(intervals = ivs, output_path = parallel[["random"]], output_format = "fasta", seed = 3)
    gsynth.replace_kmer("CG", "GC", intervals = ivs, output_path = parallel[["replace"]], output_format = "fasta")
    # gsynth.replace_kmer writes no .fai, in one call or in chunks
    expect_setequal(setdiff(list.files(tempdir()), before), basename(c(parallel, paste0(parallel[["random"]], ".fai"))))

    # the chunks are joined byte for byte
    expect_identical(readBin(parallel[["replace"]], "raw", 1e6), readBin(serial[["replace"]], "raw", 1e6))

    # the chunks hold the sequences the vector output gives with the same seed
    fa <- parallel[["random"]]
    l <- readLines(fa)
    is_header <- startsWith(l, ">")
    records <- vapply(split(l[!is_header], cumsum(is_header)[!is_header]), paste, character(1), collapse = "")
    expect_equal(unname(records), gsynth.random(intervals = ivs, output_format = "vector", seed = 3))

    # the same records, lengths and layout as the serial output: the same index
    expect_equal(readLines(paste0(fa, ".fai")), readLines(paste0(serial[["random"]], ".fai")))
    # each offset points at the first base of its record
    fai <- utils::read.table(paste0(fa, ".fai"), sep = "\t")
    con <- file(fa, "rb")
    first_lines <- vapply(seq_len(nrow(fai)), function(i) {
        seek(con, fai$V3[i])
        readChar(con, fai$V4[i])
    }, character(1))
    close(con)
    expect_equal(first_lines, l[which(is_header) + 1])
})
