.onLoad <- function(lib, pkg) {
    Sys.umask("0002")

    options(.ginteractive = FALSE)
    options(.gautocompletion = FALSE)
    options(gmax.data.size = 10000000)
    options(gmax.mem.usage = 10000000) # in KB
    options(gbig.intervals.size = 1000000)
    options(gbuf.size = 1000)
    options(gmax.processes = 16)
    options(gmax.processes2core = 2)
    options(gmin.scope4process = 10000)
    options(gmultitasking = TRUE)

    options(gquantile.edge.data.size = 100000)
    options(gpv.middle.size = 0.96)
    options(gpv.middle.precision = 10^(-4))
    options(gpv.edge.precision = 10^(-9))

    options(gtrack.chunk.size = 100000)
    options(gtrack.num.chunks = 0)

    # set the groot to samples dir
    if (!exists("GROOT", envir = .GlobalEnv)) {
        gdb.init_examples()
    }
}

.onAttach <- function(lib, pkg) {
    assign(".GFUNCS", getNamespaceExports("misha"), envir = .GlobalEnv)
    assign("GITERATOR.INTERVALS", NULL, envir = .GlobalEnv)

    assign(".GLIBDIR", path.package("misha"), envir = .GlobalEnv)
}

.onDetach <- function(lib) {

}

.onUnload <- function(lib) {
    if (exists(".GFUNCS", envir = .GlobalEnv)) {
        remove(".GFUNCS", envir = .GlobalEnv)
    }

    if (exists(".GLIBDIR", envir = .GlobalEnv)) {
        remove(".GLIBDIR", envir = .GlobalEnv)
    }

    if (exists("ALLGENOME", envir = .GlobalEnv)) {
        remove("ALLGENOME", envir = .GlobalEnv)
    }

    if (exists("GINTERVID", envir = .GlobalEnv)) {
        remove("GINTERVID", envir = .GlobalEnv)
    }

    if (exists("GITERATOR.INTERVALS", envir = .GlobalEnv)) {
        remove("GITERATOR.INTERVALS", envir = .GlobalEnv)
    }

    if (exists("GROOT", envir = .GlobalEnv)) {
        remove("GROOT", envir = .GlobalEnv)
    }

    if (exists("GWD", envir = .GlobalEnv)) {
        remove("GWD", envir = .GlobalEnv)
    }

    if (exists("GTRACKS", envir = .GlobalEnv)) {
        remove(list = get("GTRACKS"), envir = .GlobalEnv)
        remove("GTRACKS", envir = .GlobalEnv)
    }

    if (exists("GINTERVS", envir = .GlobalEnv)) {
        remove(list = get("GINTERVS"), envir = .GlobalEnv)
        remove("GINTERVS", envir = .GlobalEnv)
    }
}
