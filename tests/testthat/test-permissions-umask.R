# Permissions of the files and directories misha creates.
#
# misha writes into databases that are shared between the members of a group,
# so its own writes go through .gwith_umask() (R/utils.R), which applies the
# gpermissions.umask option and restores the process umask straight after.
# Before 5.11.21 the umask was instead set once in .onLoad() and never put
# back, which changed the permissions of everything the R session wrote.
#
# These tests set the process umask. Every one of them registers the restore
# with withr::defer() on the line after the change, so a failing expectation
# cannot leak a modified umask into the tests that follow it.

restore_groot_on_exit()

local_umask <- function(value, envir = parent.frame()) {
    old <- Sys.umask(value)
    withr::defer(Sys.umask(old), envir = envir)
    invisible(old)
}

octal_mode <- function(path) {
    if (!file.exists(path)) {
        return(sprintf("<absent: %s>", path))
    }
    format(as.octmode(file.info(path)$mode), width = 4L)
}

# A throwaway two-contig database, built by gdb.create() so the database
# skeleton itself is one of the artefacts under test.
local_scratch_db <- function(envir = parent.frame()) {
    dir <- withr::local_tempdir(.local_envir = envir)
    fasta <- file.path(dir, "genome.fa")
    set.seed(17)
    writeLines(c(
        ">chr1", paste(sample(c("A", "C", "G", "T"), 4000, TRUE), collapse = ""),
        ">chr2", paste(sample(c("A", "C", "G", "T"), 4000, TRUE), collapse = "")
    ), fasta)
    db <- file.path(dir, "db")
    gdb.create(db, fasta)
    gdb.init(db)
    db
}

# Creates one of every artefact misha writes from R and returns a named vector
# of octal modes. Track data itself is written by C++, which sets its own
# umask, so it is included as a control: it must not move.
db_artefact_modes <- function(db) {
    intervs <- data.frame(chrom = "chr1", start = seq(0, 900, by = 100), end = seq(50, 950, by = 100))

    gdir.create("sub", showWarnings = FALSE)
    gtrack.create("dense1", "d", "1", iterator = 50)
    gtrack.var.set("dense1", "tv", 1:3)
    gintervals.save("small_set", intervs)
    gintervals.attr.set("small_set", "a", "b")
    withr::with_options(list(gbig.intervals.size = 3), gintervals.save("big_set", intervs))
    gtrack.copy("dense1", "sub.dense_copy")
    gdb.mark_cache_dirty()
    gdb.reload(rescan = TRUE)

    dataset <- file.path(dirname(db), "dataset")
    gdataset.save(dataset, "a dataset", tracks = "dense1", intervals = "small_set", symlinks = FALSE)

    c(
        "db root (dir)"          = octal_mode(db),
        "db tracks (dir)"        = octal_mode(file.path(db, "tracks")),
        "chrom_sizes.txt"        = octal_mode(file.path(db, "chrom_sizes.txt")),
        ".ro_attributes"         = octal_mode(file.path(db, ".ro_attributes")),
        ".db.cache"              = octal_mode(file.path(db, ".db.cache")),
        "gdir.create (dir)"      = octal_mode(file.path(db, "tracks", "sub")),
        "track (dir)"            = octal_mode(file.path(db, "tracks", "dense1.track")),
        "track .attributes"      = octal_mode(file.path(db, "tracks", "dense1.track", ".attributes")),
        "track vars (dir)"       = octal_mode(file.path(db, "tracks", "dense1.track", "vars")),
        "track var"              = octal_mode(file.path(db, "tracks", "dense1.track", "vars", "tv")),
        "copied track (dir)"     = octal_mode(file.path(db, "tracks", "sub", "dense_copy.track")),
        "copied track file"      = octal_mode(file.path(db, "tracks", "sub", "dense_copy.track", ".attributes")),
        "small interval set"     = octal_mode(file.path(db, "tracks", "small_set.interv")),
        "interval set .iattr"    = octal_mode(file.path(db, "tracks", "small_set.iattr")),
        "big interval set (dir)" = octal_mode(file.path(db, "tracks", "big_set.interv")),
        "big interval set .meta" = octal_mode(file.path(db, "tracks", "big_set.interv", ".meta")),
        "big interval set chunk" = octal_mode(file.path(db, "tracks", "big_set.interv", "chr1")),
        "dataset (dir)"          = octal_mode(dataset),
        "dataset misha.yaml"     = octal_mode(file.path(dataset, "misha.yaml")),
        "dataset chrom_sizes"    = octal_mode(file.path(dataset, "chrom_sizes.txt")),
        "dataset track (dir)"    = octal_mode(file.path(dataset, "tracks", "dense1.track")),
        "dataset interval set"   = octal_mode(file.path(dataset, "tracks", "small_set.interv"))
    )
}

# What misha's default (gpermissions.umask = "0007") must produce, for every
# artefact and whatever the process umask is: group rwx, no world access.
MODES_DEFAULT <- c(
    "db root (dir)"          = "0770",
    "db tracks (dir)"        = "0770",
    "chrom_sizes.txt"        = "0660",
    ".ro_attributes"         = "0660",
    ".db.cache"              = "0660",
    "gdir.create (dir)"      = "0770",
    "track (dir)"            = "0770",
    "track .attributes"      = "0660",
    "track vars (dir)"       = "0770",
    "track var"              = "0660",
    "copied track (dir)"     = "0770",
    "copied track file"      = "0660",
    "small interval set"     = "0660",
    "interval set .iattr"    = "0660",
    "big interval set (dir)" = "0770",
    "big interval set .meta" = "0660",
    "big interval set chunk" = "0660",
    "dataset (dir)"          = "0770",
    "dataset misha.yaml"     = "0660",
    "dataset chrom_sizes"    = "0660",
    "dataset track (dir)"    = "0770",
    "dataset interval set"   = "0660"
)

# What misha produced before 5.11.21, reproducible with
# options(gpermissions.umask = "0002"). The 770/660 entries are the ones the
# C++ layer writes (or copies of them): it sets umask(07) itself and the
# option never applied to them.
MODES_LEGACY_0002 <- c(
    "db root (dir)"          = "0775",
    "db tracks (dir)"        = "0775",
    "chrom_sizes.txt"        = "0664",
    ".ro_attributes"         = "0664",
    ".db.cache"              = "0664",
    "gdir.create (dir)"      = "0775",
    "track (dir)"            = "0770",
    "track .attributes"      = "0660",
    "track vars (dir)"       = "0775",
    "track var"              = "0664",
    "copied track (dir)"     = "0775",
    "copied track file"      = "0660",
    "small interval set"     = "0664",
    "interval set .iattr"    = "0664",
    "big interval set (dir)" = "0775",
    "big interval set .meta" = "0664",
    "big interval set chunk" = "0664",
    "dataset (dir)"          = "0775",
    "dataset misha.yaml"     = "0664",
    "dataset chrom_sizes"    = "0664",
    "dataset track (dir)"    = "0770",
    "dataset interval set"   = "0664"
)

# Note on which of these actually discriminate: with the wrapper removed, the
# artefacts follow the process umask, so only the umask-022 run tells the two
# apart - under umask 007 the ambient umask happens to give the same 660/770.
# Do not drop the 022 variant as a duplicate of the 007 one. The two
# "process umask is unchanged" tests cannot fail against a wrapper that does
# nothing at all; they guard the opposite mistake, a wrapper that sets the
# umask and forgets to put it back.
test_that("misha's own writes are 660/770 whatever the process umask is (umask 022)", {
    local_umask("0022")
    db <- local_scratch_db()
    expect_equal(db_artefact_modes(db), MODES_DEFAULT)
})

test_that("misha's own writes are 660/770 whatever the process umask is (umask 007)", {
    local_umask("0007")
    db <- local_scratch_db()
    expect_equal(db_artefact_modes(db), MODES_DEFAULT)
})

test_that("gpermissions.umask = '0002' reproduces the pre-5.11.21 permissions", {
    local_umask("0022")
    withr::local_options(list(gpermissions.umask = "0002"))
    db <- local_scratch_db()
    expect_equal(db_artefact_modes(db), MODES_LEGACY_0002)
})

test_that("gpermissions.umask = NULL leaves the process umask in force", {
    local_umask("0027")
    withr::local_options(list(gpermissions.umask = NULL))
    expect_null(getOption("gpermissions.umask"))
    db <- local_scratch_db()
    expect_equal(octal_mode(file.path(db, "chrom_sizes.txt")), "0640")
    expect_equal(octal_mode(file.path(db, "tracks")), "0750")

    intervs <- data.frame(chrom = "chr1", start = 0, end = 100)
    gintervals.save("null_set", intervs)
    expect_equal(octal_mode(file.path(db, "tracks", "null_set.interv")), "0640")
    expect_equal(as.character(Sys.umask(NA)), "27")
})

test_that("the process umask is unchanged after a misha call returns", {
    local_umask("0022")
    db <- local_scratch_db()
    before <- as.character(Sys.umask(NA))

    gintervals.save("ok_set", data.frame(chrom = "chr1", start = 0, end = 100))
    expect_equal(as.character(Sys.umask(NA)), before)

    gtrack.create("t1", "d", "1", iterator = 50)
    expect_equal(as.character(Sys.umask(NA)), before)
})

test_that("the process umask is restored when a write errors inside the scope", {
    local_umask("0022")
    db <- local_scratch_db()
    before <- as.character(Sys.umask(NA))

    # The directory exists, so gintervals.save() gets as far as
    # .gintervals.save_file(), whose file() open then fails inside the scope.
    gdir.create("readonly", showWarnings = FALSE)
    ro <- file.path(db, "tracks", "readonly")
    Sys.chmod(ro, "0500")
    withr::defer(Sys.chmod(ro, "0700"))

    expect_error(suppressWarnings(
        gintervals.save("readonly.nope", data.frame(chrom = "chr1", start = 0, end = 100))
    ))
    expect_equal(as.character(Sys.umask(NA)), before)

    expect_error(misha:::.gwith_umask(stop("boom")), "boom")
    expect_equal(as.character(Sys.umask(NA)), before)
})

test_that("the process umask is restored when a wrapped expression is interrupted", {
    local_umask("0022")
    before <- as.character(Sys.umask(NA))
    inside <- NA_character_

    result <- tryCatch(
        misha:::.gwith_umask({
            inside <- as.character(Sys.umask(NA))
            tools::pskill(Sys.getpid(), tools::SIGINT)
            # R only acts on the pending signal at its next check point, which
            # is in this loop - i.e. still inside the umask scope.
            for (i in 1:1e7) {
                x <- sqrt(i)
            }
            "not interrupted"
        }),
        interrupt = function(i) "interrupted"
    )

    expect_equal(result, "interrupted")
    expect_equal(inside, "7") # the umask really was in force when it happened
    expect_equal(as.character(Sys.umask(NA)), before)
})

test_that("a forked child cannot leak its umask back to the parent", {
    local_umask("0022")
    before <- as.character(Sys.umask(NA))

    kids <- parallel::mclapply(1:2, function(i) {
        Sys.umask("0077")
        as.character(Sys.umask(NA))
    }, mc.cores = 2)

    expect_equal(unlist(kids), c("77", "77"))
    expect_equal(as.character(Sys.umask(NA)), before)

    db <- local_scratch_db()
    gtrack.create("mt", "d", "1", iterator = 10)
    withr::local_options(list(gmultitasking = TRUE, gmax.processes = 4, gmin.scope4process = 1))
    invisible(gextract("mt", gintervals(1), iterator = 10))
    expect_equal(as.character(Sys.umask(NA)), before)
})

test_that("an unusable gpermissions.umask value is reported, not silently applied", {
    local_umask("0022")
    before <- as.character(Sys.umask(NA))
    withr::local_options(list(gpermissions.umask = "not-a-umask"))
    expect_error(misha:::.gwith_umask(1), "Invalid gpermissions.umask option")
    expect_equal(as.character(Sys.umask(NA)), before)
})

test_that(".onLoad does not touch the process umask", {
    # The regression this file exists for: misha must never set the umask
    # anywhere that outlives the call.
    expect_false(any(grepl("Sys.umask", deparse(misha:::.onLoad), fixed = TRUE)))
    expect_false(any(grepl("Sys.umask", deparse(misha:::.onAttach), fixed = TRUE)))
})
