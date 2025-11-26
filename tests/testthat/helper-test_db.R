load_test_db <- function() {
    db_path <- if (getOption("gmulticontig.indexed_format", FALSE)) {
        "/net/mraid20/ifs/wisdom/tanay_lab/tgdata/db/tgdb/misha_test_db_indexed/"
    } else {
        "/net/mraid20/export/tgdata/db/tgdb/misha_test_db/"
    }

    # Skip if database doesn't exist
    if (!dir.exists(db_path)) {
        return(invisible(NULL))
    }

    if (getOption("misha.test.verbose", FALSE)) {
        if (getOption("gmulticontig.indexed_format", FALSE)) {
            message("Loading indexed test database")
        } else {
            message("Loading per-chromosome test database")
        }
    }
    gsetroot(db_path)

    # remove temp directory if it exists
    db_dir <- .misha$GROOT
    temp_dir <- file.path(db_dir, "tracks", "temp")
    if (dir.exists(temp_dir)) {
        unlink(temp_dir, recursive = TRUE)
    }
    gdb.reload()
    gdir.create("temp", showWarnings = FALSE)
}

#' Save and restore the current database state
#'
#' This is useful for tests that temporarily switch to a different database.
#' Uses withr-style automatic cleanup.
#'
#' @param env The environment to use for defer (defaults to parent frame)
local_db_state <- function(env = parent.frame()) {
    # Save current state
    original_groot <- if (exists("GROOT", envir = .misha, inherits = FALSE)) {
        get("GROOT", envir = .misha)
    } else {
        NULL
    }

    # Register cleanup
    withr::defer(
        {
            if (!is.null(original_groot)) {
                suppressMessages(gdb.init(original_groot))
            }
        },
        envir = env
    )
}


#' Create an isolated test database for this test file
#'
#' Creates a temporary copy of the test database with:
#' - Symlinked seq/ directory (read-only, large)
#' - Symlinked test fixture tracks (read-only)
#' - Copied small files (chrom_sizes.txt, intervs/, pssms/)
#' - Fresh empty tracks/ directory for test-created tracks
#'
#' This approach provides complete isolation between parallel test processes
#' while minimizing disk space and setup time.
#'
#' @return Path to the isolated test database
create_isolated_test_db <- function() {
    source_db <- if (getOption("gmulticontig.indexed_format", FALSE)) {
        "/net/mraid20/ifs/wisdom/tanay_lab/tgdata/db/tgdb/misha_test_db_indexed/"
    } else {
        "/net/mraid20/export/tgdata/db/tgdb/misha_test_db/"
    }

    # Create unique temp dir for this test file/process
    testdb_dir <- tempfile(pattern = "misha_testdb_", tmpdir = tempdir())
    dir.create(testdb_dir, showWarnings = FALSE)

    # Copy small files that might be read
    file.copy(
        file.path(source_db, "chrom_sizes.txt"),
        testdb_dir
    )

    # Copy interval sets and PSSM directories
    system(sprintf("cp -r %s/intervs %s/", source_db, testdb_dir))
    system(sprintf("cp -r %s/pssms %s/", source_db, testdb_dir))

    # Symlink the large seq directory (read-only)
    system(sprintf("ln -s %s/seq %s/seq", source_db, testdb_dir))

    # Create tracks directory for this test file
    tracks_dir <- file.path(testdb_dir, "tracks")
    dir.create(tracks_dir, showWarnings = FALSE)

    # Symlink test fixture tracks (read-only, used by tests)
    # Get list of non-temp tracks from source
    source_tracks <- list.dirs(file.path(source_db, "tracks"),
        full.names = FALSE, recursive = FALSE
    )
    source_tracks <- source_tracks[!grepl("^temp", source_tracks)]

    # Symlink each fixture track
    for (track in source_tracks) {
        system(sprintf(
            "ln -s %s/tracks/%s %s/tracks/%s",
            source_db, track, testdb_dir, track
        ))
    }

    # Set as root and reload to recognize symlinked tracks
    gsetroot(testdb_dir)
    gdb.reload()

    withr::defer(
        {
            unlink(testdb_dir, recursive = TRUE)
        },
        envir = parent.frame()
    )

    invisible(testdb_dir)
}
