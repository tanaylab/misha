.onLoad <- function(lib, pkg) {
}

.onAttach <- function(lib, pkg) {
	Sys.umask("0002")

	assign(".GFUNCS", getNamespaceExports("misha"), envir = .GlobalEnv)
	assign("GITERATOR.INTERVALS", NULL, envir = .GlobalEnv)

	if (R.Version()$major >= 3)
		assign(".GLIBDIR", path.package("misha"), envir = .GlobalEnv)
	else
		assign(".GLIBDIR", .path.package("misha"), envir = .GlobalEnv)
	
	options(.ginteractive = FALSE)
	options(.gautocompletion = FALSE)
	options(gmax.data.size = 10000000)
	options(gmax.mem.usage = 10000000)    # in KB
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

	# if grnd.seed==0, the seed is determined internally in a random manner
	options(grnd.seed = 1)

	options(gtrack.chunk.size = 100000)
	options(gtrack.num.chunks = 0)

	# set the groot to samples dir
	gdb.init_examples()
}

.onDetach <- function(lib) {
	if (exists(".GFUNCS", envir = .GlobalEnv))
		remove(".GFUNCS", envir = .GlobalEnv)

	if (exists(".GLIBDIR", envir = .GlobalEnv))
		remove(".GLIBDIR", envir = .GlobalEnv)

	if (exists("ALLGENOME", envir = .GlobalEnv))
		remove("ALLGENOME", envir = .GlobalEnv)

	if (exists("GINTERVID", envir = .GlobalEnv))
		remove("GINTERVID", envir = .GlobalEnv)

	if (exists("GITERATOR.INTERVALS", envir = .GlobalEnv))
		remove("GITERATOR.INTERVALS", envir = .GlobalEnv)

	if (exists("GROOT", envir = .GlobalEnv))
		remove("GROOT", envir = .GlobalEnv)

	if (exists("GWD", envir = .GlobalEnv))
		remove("GWD", envir = .GlobalEnv)

	if (exists("GTRACKS", envir = .GlobalEnv)) {
		remove(list = get("GTRACKS"), envir = .GlobalEnv)
		remove("GTRACKS", envir = .GlobalEnv)
	}

	if (exists("GINTERVS", envir = .GlobalEnv)) {
		remove(list = get("GINTERVS"), envir = .GlobalEnv)
		remove("GINTERVS", envir = .GlobalEnv)
	}
}

.gcall <- function(...) {
	tryCatch({ res <- .Call(...) },
			 interrupt = function(interrupt){ stop("Command interrupted!", call. = FALSE); } )
	res
}

.gcall_noninteractive <- function(FUN, ...) {
	.ginteractive = .ggetOption(".ginteractive")
	tryCatch({
		options(.ginteractive = F)
		do.call(FUN, list(...))
	},
	finally = {
		options(.ginteractive = .ginteractive)
	})
}

.gcheckroot <- function() {
	if (!exists("GROOT", envir = .GlobalEnv) || !exists("ALLGENOME", envir = .GlobalEnv) || is.null(get("GROOT")) || is.null(get("ALLGENOME")))
		stop("Database root directory is not set. Please call gdb.init().", call. = F)
}

.gdir.cd <- function(dir, rescan) {
	oldwd <- getwd()
	setwd(get("GWD"))
	tryCatch({
		t <- .gfindtrackinpath(dir)
		if (!is.null(t))
			stop(sprintf("Directory %s belongs to track %s", dir, t), call. = F)

		setwd(dir)
		newwd <- getwd()

		if (.ggetOption(".gautocompletion", FALSE))
			.gundefine_autocompletion_vars()
		assign("GWD", newwd, envir = .GlobalEnv)
		setwd(oldwd)
		gdb.reload(rescan)
	},
	interrupt = function(interrupt){ setwd(oldwd) },
	finally = {
		setwd(oldwd)
	})
}

.ggetOption <- function(x, default = NULL) {
    if (missing(default)) 
        return(options(x)[[1L]])
    if (x %in% names(options())) 
        options(x)[[1L]]
    else default
}

.gexpr2str <- function(x) {
	if (.ggetOption(".ginteractive", FALSE)) {
		if (is.null(substitute(x)))
			NULL
		else {
			str <- deparse(substitute(x), width.cutoff = 500)[1]
			gsub("^\"(.*)\"$", "\\1", str, perl = T)
		}
	} else
		eval.parent(x)
}

.giterator <- function(iterator) {
	if (typeof(iterator) == "integer" || typeof(iterator) == "double")
		return(iterator)
	
	iterator.str <- do.call(.gexpr2str, list(substitute(iterator)), envir = parent.frame())

	if (typeof(iterator.str) == "character") {
		if (!is.na(match(iterator.str, get("GTRACKS"))) || !is.na(match(iterator.str, get("GINTERVS"))))
			return(iterator.str)
	}
	
	iterator
}

.grbind <- function(...) {
	objs <- list(...)

	zerolines <- lapply(objs, function(obj) { obj[0,] })

	diffs <- sapply(zerolines, FUN = attr.all.equal, zerolines[[1]])
	if (!all(sapply(diffs, FUN = is.null)))
	   stop("Cannot rbind objects: columns differ", call. = F)

	.gcall("grbind", objs, new.env(parent = parent.frame()))
}

.gverify_max_data_size <- function(size, data_name = "Result", arguments = NULL) {
	max.data.size <- .ggetOption("gmax.data.size", 10000000)

	if (size > max.data.size) {
		if (is.null(arguments))
			stop(sprintf(
				paste("%s size exceeded the maximal allowed (%d).",
					"Note: the maximum data size is controlled via gmax.data.size option (see options, getOptions).", sep = "\n"),
				data_name, max.data.size), call. = F)
		else
			stop(sprintf(
				paste("%s size exceeded the maximal allowed (%d).",
					"Consider saving the result in a file (use %s argument).",
					"Note: the maximum data size is controlled via gmax.data.size option (see options, getOptions).", sep = "\n"),
				data_name, max.data.size, paste(arguments, collapse = " or ")), call. = F)
	}
}

.gvtrack <- function(vtrack) {
	vtrackstr <- do.call(.gexpr2str, list(substitute(vtrack)), envir = parent.frame())
	if (!is.character(vtrackstr) || length(vtrackstr) != 1)
		stop(sprintf("Virtual track must be specified as a character string"), call. = F)

    if (is.na(match(vtrackstr, gvtrack.ls())))
        stop(sprintf("Virtual track %s does not exist", vtrackstr), call. = F)
	vtrackstr
}

.gvtrack.get <- function(vtrackstr) {
    if (is.na(match(vtrackstr, gvtrack.ls())))
        stop(sprintf("Virtual track %s does not exist", vtrackstr), call. = F)

    gwd <- get("GWD")
	get("GVTRACKS")[[gwd]][[vtrackstr]]
}

.gvtrack.set <- function(vtrackstr, var) {
    if (exists("GVTRACKS", envir = .GlobalEnv))
        gvtracks <- get("GVTRACKS")
    else
        gvtracks <- list()

    gwds <- names(gvtracks)
	if (!is.list(gvtracks) || (length(gvtracks) && !is.character(gwds)) || length(gvtracks) != length(gwds))
		stop("Invalid format of GVTRACKS variable.\nTo continue working with virtual tracks please remove this variable from the environment.", call. = F)

    gwd <- get("GWD")
    idx1 <- match(gwd, gwds)
    if (is.na(idx1)) {
        gwds <- c(gwds, gwd)
        idx1 <- length(gwds)
        gvtracks[[idx1]] <- list()
    }

    vtracks <- gvtracks[[idx1]]
    names(gvtracks) <- gwds

    vtracknames <- names(vtracks)
	if (!is.list(vtracks) || (length(vtracks) && !is.character(vtracknames)) || length(vtracks) != length(vtracknames))
		stop("Invalid format of GVTRACKS variable.\nTo continue working with virtual tracks please remove this variable from the environment.", call. = F)

    idx2 <- match(vtrackstr, vtracknames)
    if (is.na(idx2)) {
		if (!is.na(match(vtrackstr, get("GTRACKS"))))
			stop(sprintf("Track %s already exists", vtrackstr), call. = F)

		if (!is.na(match(vtrackstr, get("GINTERVS"))))
			stop(sprintf("Interval %s already exists", vtrackstr), call. = F)

 		if (.ggetOption(".gautocompletion", FALSE) && exists(vtrackstr, envir = .GlobalEnv))
			stop(sprintf("Variable \"%s\" shadows the name of identically named virtual track.\nPlease remove this variable from the environment or switch off autocompletion mode.", vtrackstr), call. = F)

        vtracknames <- c(vtracknames, vtrackstr)
        idx2 <- length(vtracknames)
    }

    gvtracks[[idx1]][[idx2]] <- var

	if (!is.null(var)) {
		names(gvtracks[[idx1]]) <- vtracknames
		
		envir <- new.env(parent = parent.frame())
		assign("GVTRACKS", gvtracks, envir)
		.gcall("gcheck_vtrack", vtrackstr, envir)
	}

	success <- F
	old.gvtracks <- NULL
	if (exists("GVTRACKS", envir = .GlobalEnv))
		old.gvtracks <- get("GVTRACKS")

	success <- F
	tryCatch({
		if (.ggetOption(".gautocompletion", FALSE))
			.gundefine_autocompletion_vars()
			assign("GVTRACKS", gvtracks, envir = .GlobalEnv)
		if (.ggetOption(".gautocompletion", FALSE))
			.gdefine_autocompletion_vars(get("GTRACKS"), get("GINTERVS"), gvtrack.ls(), .ggetOption(".ginteractive", FALSE))
		success <- T
	},
	finally = {
		if (!success && .ggetOption(".gautocompletion", FALSE))
			.gdefine_autocompletion_vars(get("GTRACKS"), get("GINTERVS"), gvtrack.ls(), .ggetOption(".ginteractive", FALSE))
	})
}

.gchroms <- function(chroms) {
	if (!is.character(chroms))
		chroms <- as.character(chroms)

	idx <- substr(chroms, 1, 3) != "chr"
	chroms[idx] <- paste("chr", chroms[idx], sep = "")
	
	indices <- match(chroms, get("ALLGENOME")[[1]]$chrom)
	err.chroms <- chroms[is.na(indices)]
	if (length(err.chroms) > 0)
		stop(sprintf("Chromosome %s does not exist in the database", err.chroms[1]))
	get("ALLGENOME")[[1]]$chrom[indices]   # return factor
}

.gintervals <- function(chroms, starts, ends, strands) {
	if (is.null(strands))
    	intervals <- data.frame(chrom = .gchroms(chroms), start = starts, end = ends)
	else
	    intervals <- data.frame(chrom = .gchroms(chroms), start = starts, end = ends, strand = strands)

	numintervals <- nrow(intervals)

	maxends <- get("ALLGENOME")[[1]]$end[match(as.character(intervals$chrom), get("ALLGENOME")[[1]]$chrom)]
	maxidx <- intervals$end == -1
	intervals$end[maxidx] <- maxends[maxidx]
	
	err.intervs <- intervals[intervals$start < 0, ]
	if (nrow(err.intervs) > 0)
		stop(sprintf("Invalid interval (%s, %g, %g): start coordinate is out of range", err.intervs$chrom[1], err.intervs$start[1], err.intervs$end[1]), call. = F)
	
	err.intervs <- intervals[intervals$end > maxends, ]
	if (nrow(err.intervs) > 0)
		stop(sprintf("Invalid interval (%s, %g, %g): end coordinate exceeds chromosome boundaries", err.intervs$chrom[1], err.intervs$start[1], err.intervs$end[1]), call. = F)
	
	err.intervs <- intervals[intervals$start >= intervals$end, ]
	if (nrow(err.intervs) > 0)
		stop(sprintf("Invalid interval (%s, %g, %g): start coordinate exceeds or equals to end coordinate", err.intervs$chrom[1], err.intervs$start[1], err.intervs$end[1]), call. = F)

	if (!is.null(strands)) {
	    if (!is.numeric(intervals$strand))
		    stop("Invalid strand values", call. = F)
			
	    err.intervs <- intervals[intervals$strand != as.integer(intervals$strand) | intervals$strand < -1 | intervals$strand > 1, ]
		if (nrow(err.intervs) > 0)
		    stop(sprintf("Invalid strand value %g of interval (%s, %g, %g)", err.intervs$strand[1], err.intervs$chrom[1], err.intervs$start[1], err.intervs$end[1]))
	}

	intervals
}

.gintervals.apply <- function(chroms, intervals, intervals.set.out, FUN, ...) {
	if (!is.null(intervals.set.out))
		fullpath <- .gintervals.check_new_set(intervals.set.out)

	if (is.data.frame(intervals)) 
		intervals <- list(intervals)

	# sort chroms
	chroms$size <- NULL
	if ("chrom" %in% colnames(chroms))
		chroms <- data.frame(chrom = chroms[with(chroms, order(chrom)), ])
	else
		chroms <- chroms[with(chroms, order(chrom1, chrom2)), ]


	# let's assume that if any of the input intervals sets is big then intervals.set.out should be big as well
	if (any(unlist(lapply(intervals, function(intervals) { .gintervals.is_bigset(intervals) || .gintervals.needs_bigset(intervals) })))) {
		stats <- NULL
		zeroline <- NULL
		success <- FALSE
		t <- Sys.time()
		progress.percentage <- -1
		tryCatch({
			# if any of the source intervals sets is big then create the output intervals set big too
			if (!is.null(intervals.set.out))
				dir.create(fullpath, recursive = T, mode = "0777")

			if (.gintervals.is1d(intervals[[1]])) {
				mapply(function(chrom) {
						loaded_intervals <- lapply(intervals, function(intervals) { .gintervals.load_ext(intervals, chrom = chrom) })
						res <- do.call(FUN, list(loaded_intervals, ...))
						if (!is.null(intervals.set.out) && !is.null(res) && nrow(res) > 0) {
							zeroline <<- res[0, ]
							.gintervals.big.save(fullpath, res, chrom = chrom)
							stat <- .gcall("gintervals_stats", res, new.env(parent = parent.frame()))
							stats <<- rbind(stats, data.frame(chrom = chrom, stat))
						}
						if (as.integer(difftime(Sys.time(), t, units="secs")) > 3) {
							t <<- Sys.time()
							percentage <- as.integer(100 * match(chrom, chroms$chrom) / nrow(chroms))
							if (percentage < 100 && progress.percentage != percentage) {
								cat(sprintf("%d%%...", percentage))
								progress.percentage <<- percentage
							}
						}
					}, chroms$chrom)
			} else {
				mapply(function(chrom1, chrom2) {
						loaded_intervals <- lapply(intervals, function(intervals) { .gintervals.load_ext(intervals, chrom1 = chrom1, chrom2 = chrom2) })
						res <- do.call(FUN, list(loaded_intervals, ...))
						if (!is.null(intervals.set.out) && !is.null(res) && nrow(res) > 0) {
							zeroline <<- res[0, ]
							.gintervals.big.save(fullpath, res, chrom1 = chrom1, chrom2 = chrom2)
							stat <- .gcall("gintervals_stats", res, new.env(parent = parent.frame()))
							stats <<- rbind(stats, data.frame(chrom1 = chrom1, chrom2 = chrom2, stat))
						}
						if (as.integer(difftime(Sys.time(), t, units="secs")) > 3) {
							t <<- Sys.time()
							percentage <- as.integer(100 * which(chroms$chrom1 == chrom1 & chroms$chrom2 == chrom2) / nrow(chroms))
							if (percentage < 100 && progress.percentage != percentage) {
								cat(sprintf("%d%%...", percentage))
								progress.percentage <<- percentage
							}
						}
					}, chroms$chrom1, chroms$chrom2)
			}

			if (!is.null(intervals.set.out)) {
				if (is.null(stats))
					return(retv <- NULL)
				.gintervals.big.save_meta(fullpath, stats, zeroline)
			}

			if (progress.percentage >= 0)
				cat("100%\n")

			success <- TRUE

			# check whether the output intervals set needs to remain in big format
			if (!is.null(intervals.set.out) && !.gintervals.needs_bigset(intervals.set.out))
				.gintervals.big2small(intervals.set.out)
		},
		finally = {
			if (!success && !is.null(intervals.set.out))
				unlink(fullpath, recursive = TRUE)
		})
	} else {
		loaded_intervals <- lapply(intervals, .gintervals.load_ext)
		res <- do.call(FUN, list(loaded_intervals, ...))
		if (!is.null(intervals.set.out) && !is.null(res) && nrow(res) > 0) {
			if (.gintervals.is1d(res))
				res <- res[order(res$chrom), ]
			else
				res <- res[order(res$chrom1, res$chrom2), ]
			if (.gintervals.needs_bigset(res))
				.gintervals.small2big(intervals.set.out, res)
			else
				.gintervals.save_file(fullpath, res)
		} else
			return(NULL)
	}

	# refresh the list of GINTERVS, etc.
	if (!is.null(intervals.set.out))
		.gdb.add_intervals.set(intervals.set.out)
}

.gintervals.check_new_set <- function(intervals.set) {
	if (!is.na(match(intervals.set, get("GINTERVS"))))
		stop(sprintf("Intervals set %s already exists", intervals.set), call. = F)
	
	if (!length(grep("^[A-Za-z][\\w.]*$", intervals.set, perl=T)))
		stop("Invalid interval name %s. Only alphanumeric characters and _ are allowed in the name.")

	path <- gsub(".", "/", intervals.set, fixed = T)
	path <- paste(path, ".interv", sep = "")
	fullpath <- paste(get("GWD"), path, sep = "/")
	dir <- dirname(path)
	fulldir <- paste(get("GWD"), dir, sep = "/")

	if (!file.exists(fulldir))
		stop(sprintf("Directory %s does not exist", dir), call. = F)

	if (file.exists(fullpath))
		stop(sprintf("File %s already exists", path), call. = F)

	if (!is.na(match(intervals.set, get("GTRACKS"))))
		stop(sprintf("Track %s already exists", intervals.set), call. = F)

	if (!is.na(match(intervals.set, gvtrack.ls())))
		stop(sprintf("Virtual track %s already exists", intervals.set), call. = F)

	if (.ggetOption(".gautocompletion", FALSE) && exists(intervals.set, envir = .GlobalEnv))
		stop(sprintf("Variable \"%s\" shadows the name of the new intervals set.\nPlease remove this variable from the environment or switch off autocompletion mode.", intervals.set), call. = F)

	fullpath
}

.gintervals.load_ext <- function(intervals.set = NULL, chrom = NULL, chrom1 = NULL, chrom2 = NULL, progress = F) {
	if (is.null(intervals.set))
		stop("Usage: gintervals.load(intervals.set, chrom = NULL, chrom1 = NULL, chrom2 = NULL)", call. = F)
	.gcheckroot()

	if (is.character(intervals.set) && length(intervals.set) == 1 && is.na(match(intervals.set, get("GINTERVS"))) && is.na(match(intervals.set, get("GTRACKS"))))
		stop(sprintf("Intervals set %s does not exist", intervals.set), call. = F)

	.gintervals.load(intervals.set, chrom, chrom1, chrom2, progress)
}

.gintervals.load <- function(intervals.set = NULL, chrom = NULL, chrom1 = NULL, chrom2 = NULL, progress = F) {
	if (!is.null(chrom)) {
		chrom <- .gchroms(chrom)
		if (length(chrom) > 1)
			stop("chrom parameter should mark only one chromosome")
	}

	if (!is.null(chrom1)) {
		chrom1 <- .gchroms(chrom1)
		if (length(chrom1) > 1)
			stop("chrom1 parameter should mark only one chromosome")
	}

	if (!is.null(chrom2)) {
		chrom2 <- .gchroms(chrom2)
		if (length(chrom2) > 1)
			stop("chrom2 parameter should mark only one chromosome")
	}

	if (!is.null(chrom) && !is.null(chrom1))
		stop("Cannot use chrom and chrom1 parameters in the same call", call. = F)

	if (!is.null(chrom) && !is.null(chrom2))
		stop("Cannot use chrom and chrom2 parameters in the same call", call. = F)

	if (is.character(intervals.set) && length(intervals.set) != 1 || !is.character(intervals.set) && !.gintervals.is1d(intervals.set) && !.gintervals.is2d(intervals.set))
		stop("Invalid format of intervals", call. = F)

	res <- NULL
	if (.gintervals.is_bigset(intervals.set)) {
		meta <- .gintervals.big.meta(intervals.set)
		zeroline <- meta$zeroline
		t <- Sys.time()
		progress.percentage <- -1

		if (.gintervals.big.is1d(intervals.set)) {
			if (!is.null(chrom1) || !is.null(chrom2))
				stop(sprintf("%s is a 1D big intervals set.\nchrom1 or chrom2 parameters can be applied only to 2D intervals.", intervals.set), call. = F)

			if (!is.null(chrom))
				meta$stats <- meta$stats[meta$stats$chrom == chrom,]

			if (!.gintervals.loadable(intervals.set, chrom = chrom)) {
				if (is.null(chrom))
					stop(sprintf("Cannot load a big intervals set %s: its size (%d) exceeds the limit (%d) controlled by gmax.data.size option.\nFor big intervals sets only one chromosome pair can be loaded at a time.",
						intervals.set, sum(meta$stats$size), .ggetOption("gmax.data.size", 10000000)), call. = F)
				else
					stop(sprintf("Cannot load chromosome %s of an intervals set %s: its size (%d) exceeds the limit (%d) controlled by gmax.data.size option.",
						chrom, intervals.set, sum(meta$stats$size), .ggetOption("gmax.data.size", 10000000)), call. = F)
			}

			if (nrow(meta$stats) > 1) {
				res <- list(zeroline)
				lapply(meta$stats$chrom,
					function(chrom) {
						loaded_intervs <- .gintervals.load_file(intervals.set, chrom = chrom)
						if (!identical(sapply(loaded_intervs, "class"), sapply(zeroline, "class")))
							stop(sprintf("Intervals set %s, chrom %s: invalid columns definition", intervals.set, chrom), call. = F)
						res <<- c(res, list(loaded_intervs))
						if (as.integer(difftime(Sys.time(), t, units="secs")) > 3) {
							t <<- Sys.time()
							percentage <- as.integer(100 * match(chrom, meta$stats$chrom) / nrow(meta$stats$chrom))
							if (progress && percentage < 100 && progress.percentage != percentage) {
								cat(sprintf("%d%%...", percentage))
								progress.percentage <<- percentage
							}
						}
					}
				)
				res <- do.call(.grbind, res)   # much faster than calling rbind incrementally in mapply
			} else if (nrow(meta$stats) == 1) {
				res <- .gintervals.load_file(intervals.set, chrom = meta$stat$chrom[1])
				if (!identical(sapply(res, "class"), sapply(zeroline, "class")))
					stop(sprintf("Intervals set %s, chrom %s: invalid columns definition", intervals.set, meta$stat$chrom[1]), call. = F)
			} else
				res <- meta$zeroline
		} else {
			if (!is.null(chrom))
				stop(sprintf("%s is a 2D big intervals set.\nchrom parameter can be applied only to 1D intervals.", intervals.set), call. = F)

			if (!is.null(chrom1))
				meta$stats <- meta$stats[meta$stats$chrom1 == chrom1,]
			if (!is.null(chrom2))
				meta$stats <- meta$stats[meta$stats$chrom2 == chrom2,]

			if (!.gintervals.loadable(intervals.set, chrom1 = chrom1, chrom2 = chrom2)) {
				if (!is.null(chrom1) && !is.null(chrom2))
					stop(sprintf("Cannot load chromosome pair (%s, %s) of an intervals set %s: its size (%d) exceeds the limit (%d) controlled by gmax.data.size option.",
						chrom1, chrom2, intervals.set, sum(meta$stats$size), .ggetOption("gmax.data.size", 10000000)), call. = F)
				else if (!is.null(chrom1))
					stop(sprintf("Cannot load chromosome %s of an intervals set %s: its size (%d) exceeds the limit (%d) controlled by gmax.data.size option.",
						chrom1, intervals.set, sum(meta$stats$size), .ggetOption("gmax.data.size", 10000000)), call. = F)
				else if (!is.null(chrom2))
					stop(sprintf("Cannot load chromosome %s of an intervals set %s: its size (%d) exceeds the limit (%d) controlled by gmax.data.size option.",
						chrom2, intervals.set, sum(meta$stats$size), .ggetOption("gmax.data.size", 10000000)), call. = F)
				else
					stop(sprintf("Cannot load a big intervals set %s: its size (%d) exceeds the limit (%d) controlled by gmax.data.size option.\nFor big intervals sets only one chromosome pair can be loaded at a time.",
						intervals.set, sum(meta$stats$size), .ggetOption("gmax.data.size", 10000000)), call. = F)
			}

			if (nrow(meta$stats) > 1) {
				res <- list(zeroline)
				mapply(function(chrom1, chrom2) {
						loaded_intervs <- .gintervals.load_file(intervals.set, chrom1 = chrom1, chrom2 = chrom2)
						if (!identical(sapply(loaded_intervs, "class"), sapply(zeroline, "class")))
							stop(sprintf("Interval set %s, chrom1 %s, chrom2 %s: invalid columns definition", intervals.set, chrom1, chrom2), call. = F)
						res <<- c(res, list(loaded_intervs))
						if (as.integer(difftime(Sys.time(), t, units="secs")) > 3) {
							t <<- Sys.time()
							percentage <- as.integer(100 * which(meta$stats$chrom1 == chrom1 & meta$stats$chrom2 == chrom2) / nrow(meta$stats))
							if (progress && percentage < 100 && progress.percentage != percentage) {
								cat(sprintf("%d%%...", percentage))
								progress.percentage <<- percentage
							}
						}
					},
					meta$stats$chrom1, meta$stats$chrom2)
				res <- do.call(.grbind, res)   # much faster than calling rbind incrementally in mapply
			} else if (nrow(meta$stats) == 1) {
				res <- .gintervals.load_file(intervals.set, chrom1 = meta$stat$chrom1[1], chrom2 = meta$stat$chrom2[1])
				if (!identical(sapply(res, "class"), sapply(zeroline, "class")))
					stop(sprintf("Interval set %s, chrom1 %s, chrom2 %s: invalid columns definition", intervals.set, meta$stat$chrom1[1], meta$stat$chrom2[1]), call. = F)
			} else
				res <- meta$zeroline
		}

		if (progress.percentage >= 0)
			cat("100%\n")
	} else {
		if (is.character(intervals.set) && length(intervals.set) == 1)
			res <- .gintervals.load_file(intervals.set)
		else
			res <- intervals.set
		if (!is.null(res)) {
			if (!.gintervals.is1d(res) && !is.null(chrom))
				stop("chrom parameter can be applied only to 1D intervals", call. = F)

			if (!.gintervals.is2d(res) && (!is.null(chrom1) || !is.null(chrom2)))
				stop("chrom1 or chrom2 parameters can be applied only to 2D intervals", call. = F)

			if (nrow(res) > 0) {
				if (!is.null(chrom)) {
					res <- res[res$chrom == chrom, ]
					if (nrow(res))
						rownames(res) <- 1 : nrow(res)
				}
				if (!is.null(chrom1)) {
					res <- res[res$chrom1 == chrom1, ]
					if (nrow(res))
						rownames(res) <- 1 : nrow(res)
				}
				if (!is.null(chrom2)) {
					res <- res[res$chrom2 == chrom2, ]
					if (nrow(res))
						rownames(res) <- 1 : nrow(res)
				}
			}
		}
	}
	res
}

.gintervals.load_file <- function(intervals.set, chrom = NULL, chrom1 = NULL, chrom2 = NULL) {
	intervfname <- sprintf("%s.interv", paste(get("GWD"), gsub("\\.", "/", intervals.set), sep = "/"))
	if (!is.null(chrom)) {
		chrom <- .gchroms(chrom)
		intervfname <- sprintf("%s/%s", intervfname, chrom)
	} else if (!is.null(chrom1) && !is.null(chrom2)) {
		chrom1 <- .gchroms(chrom1)
		chrom2 <- .gchroms(chrom2)
		intervfname <- sprintf("%s/%s-%s", intervfname, chrom1, chrom2)
	}

	if (intervals.set %in% get("GTRACKS"))
		.gcall("gtrack_intervals_load", intervals.set, chrom, chrom1, chrom2, new.env(parent = parent.frame()))
	else {
		if (file.exists(intervfname) || (is.null(chrom) && is.null(chrom1) && is.null(chrom2))) {
			f <- file(intervfname, "rb")
			intervals.set <- unserialize(f)
			close(f)
			intervals.set
		} else {
			if (.gintervals.is_bigset(intervals.set))
				.gintervals.big.meta(intervals.set)$zeroline
			else
				stop(sprintf("File %s does not exist", intervfname), call. = F)
		}
	}
}

.gintervals.save_file <- function(filename, intervs) {
	if (nrow(intervs)) {
		if (.gintervals.is1d(intervs)) {
			intervs$chrom <- as.factor(intervs$chrom)
			point.intervs <- intervs$start == intervs$end
			intervs[point.intervs, ]$end = intervs[point.intervs, ]$end + 1
		} else {
			intervs$chrom1 <- as.factor(intervs$chrom1)
			intervs$chrom2<- as.factor(intervs$chrom2)
			point.intervs <- intervs$start1 == intervs$end1
			intervs[point.intervs, ]$end1 = intervs[point.intervs, ]$end1 + 1
			point.intervs <- intervs$start2 == intervs$end2
			intervs[point.intervs, ]$end2 = intervs[point.intervs, ]$end2 + 1
		}
	}

	f <- file(filename, "wb")
	serialize(intervs, f)
	close(f)
	nrow(intervs)
}

.gintervals.is_bigset <- function(intervals.set, err_if_non_exist = T) {
	if (is.character(intervals.set) & length(intervals.set) == 1) {
		if (intervals.set %in% get("GTRACKS"))
			intervfname <- sprintf("%s.track", paste(get("GWD"), gsub("\\.", "/", intervals.set), sep = "/"))
		else
			intervfname <- sprintf("%s.interv", paste(get("GWD"), gsub("\\.", "/", intervals.set), sep = "/"))
		if (file.exists(intervfname)) {
			if (file.info(intervfname)$isdir)
				return(TRUE)
		} else if (err_if_non_exist)
			stop(sprintf("Intervals set %s does not exist", intervals.set), call. = F)
	}
	FALSE
}

.gintervals.needs_bigset <- function(intervals = NULL, size = NULL) {
	if (!is.null(intervals)) {
		chromsizes <- gintervals.chrom_sizes(intervals)
		nrow(chromsizes) && sum(chromsizes$size) > min(.ggetOption("gmax.data.size", 10000000), .ggetOption("gbig.intervals.size", 1000000))
	} else
		size > min(.ggetOption("gmax.data.size", 10000000), .ggetOption("gbig.intervals.size", 1000000))
}

.gintervals.loadable <- function(intervals = NULL, size = NULL, chrom = NULL, chrom1 = NULL, chrom2 = NULL) {
	if (!is.null(intervals)) {
		chromsizes <- gintervals.chrom_sizes(intervals)
		if (nrow(chromsizes)) {
			if (!is.null(chrom))
				chromsizes <- chromsizes[chromsizes$chrom == chrom,]
			if (!is.null(chrom1))
				chromsizes <- chromsizes[chromsizes$chrom1 == chrom1,]
			if (!is.null(chrom2))
				chromsizes <- chromsizes[chromsizes$chrom2 == chrom2,]
		}
		!nrow(chromsizes) || sum(chromsizes$size) <= .ggetOption("gmax.data.size", 10000000)
	} else
		size <= .ggetOption("gmax.data.size", 10000000)
}

.gintervals.big.is1d <- function(intervals.set) {
	"chrom" %in% colnames(.gintervals.big.meta(intervals.set)$stats)
}

.gintervals.big.is2d <- function(intervals.set) {
	"chrom1" %in% colnames(.gintervals.big.meta(intervals.set)$stats)
}

.gintervals.big2small <- function(intervals.set) {
	# We assume that writing the intervals might be a lengthy process.
	# During this time the process might get interrupted leaving the intervals set in incomplete state.
	# Even though it's not fully transaction-safe, we prefer to create a temporary file and then move it hoping it's fast enough.
	path <- gsub(".", "/", intervals.set, fixed = T)
	path <- paste(get("GWD"), path, sep = "/")
	path <- paste(path, ".interv", sep = "")

	intervals <- .gintervals.load(intervals.set)
	tmpfilename <- tempfile(".", dirname(path)) # tmpdir = the parent directory of intervals set, otherwise rename might not work
	                                            # (tmpdir might be at another file system)
	file.rename(path, tmpfilename)
	success <- FALSE
	tryCatch({
		.gintervals.save_file(path, intervals)
		success <- TRUE
	},
	finally = {
		if (!success) {
			unlink(path, recursive = TRUE)
			file.rename(tmpfilename, path)
		}
		unlink(tmpfilename, recursive = TRUE)
	})
}

.gintervals.small2big <- function(intervals.set, intervals = NULL) {
	# We assume that writing the intervals might be a lengthy process.
	# During this time the process might get interrupted leaving the intervals set in incomplete state.
	# Even though it's not fully transaction-safe, we prefer to create a temporary file and then move it hoping it's fast enough.
	if (is.null(intervals))
		intervals <- .gintervals.load(intervals.set)

	path <- gsub(".", "/", intervals.set, fixed = T)
	path <- paste(get("GWD"), path, sep = "/")
	path <- paste(path, ".interv", sep = "")

	tmpfilename <- tempfile(".", dirname(path)) # tmpdir = the parent directory of intervals set, otherwise rename might not work
	                                            # (tmpdir might be at another file system)
	file.rename(path, tmpfilename)
	gintervs <- get("GINTERVS")
	assign("GINTERVS", gintervs[gintervs != intervals.set], envir = .GlobalEnv)
	success <- FALSE
	tryCatch({
		.gcall_noninteractive(gintervals.save, intervals.set, intervals)
		success <- TRUE
	},
	finally = {
		if (!success) {
			unlink(path, recursive = TRUE)
			file.rename(tmpfilename, path)
			gintervs <- c(gintervs, intervals.set)
			gintervs <- unique(gintervs)
			sort(gintervs)
			assign("GINTERVS", gintervs, envir = .GlobalEnv)
		}
		unlink(tmpfilename, recursive = TRUE)
	})
}

.gintervals.big.meta <- function(intervals.set) {
	metafname <- ""
	if (intervals.set %in% get("GTRACKS")) {
		metafname <- sprintf("%s.track/.meta", paste(get("GWD"), gsub("\\.", "/", intervals.set), sep = "/"))
		if (!file.exists(metafname))
			.gcall("gtrack_create_meta", intervals.set, new.env(parent = parent.frame()))
	} else
		metafname <- sprintf("%s.interv/.meta", paste(get("GWD"), gsub("\\.", "/", intervals.set), sep = "/"))
	f <- file(metafname, "rb")
	res <- unserialize(f)
	close(f)
	res
}

.gintervals.big.save_meta <- function(path, stats, zeroline) {
	f <- file(sprintf("%s/.meta", path), "wb")
	serialize(list(stats = stats, zeroline = zeroline), f)
	close(f)
}

.gintervals.big.save <- function(path, intervs, chrom = NULL, chrom1 = NULL, chrom2 = NULL) {
	if (!is.null(chrom))
		filename <- sprintf("%s/%s", path, chrom)
	else
		filename <- sprintf("%s/%s-%s", path, chrom1, chrom2)

	if (is.null(intervs) || nrow(intervs) == 0)
		unlink(filename, recursive = TRUE)
	else {
		if (!is.null(chrom))
			intervs <- intervs[intervs$chrom == chrom, ]
		else
			intervs <- intervs[intervs$chrom1 == chrom1 & intervs$chrom2 == chrom2, ]
		.gintervals.save_file(filename, intervs)
	}
}

#' @export
.gintervals.is1d <- function(intervals) {
	if (is.character(intervals)) {
		if (.gintervals.is_bigset(intervals))
			return(.gintervals.big.is1d(intervals))
		intervals <- .gintervals.load(intervals)
	}
	all(colnames(intervals)[1 : 3] == c("chrom", "start", "end"))
}

#' @export
.gintervals.is2d <- function(intervals) {
	if (is.character(intervals)) {
		if (.gintervals.is_bigset(intervals))
			return(.gintervals.big.is2d(intervals))
		intervals <- .gintervals.load(intervals)
	}
	all(colnames(intervals)[1 : 6] == c("chrom1", "start1", "end1", "chrom2", "start2", "end2"))
}

# tries to locate track name in a path, return the track name
# examples:
#   .gfindtrackinpath("aaa/bbb/ccc.track/ddd/eee")       returns "aaa.bbb.ccc"
#   .gfindtrackinpath("aaa/bbb/ccc.track/ddd.track/eee") returns "aaa.bbb.ccc"
#   .gfindtrackinpath("aaa/bbb/ccc/ddd/eee")             returns NULL
.gfindtrackinpath <- function(path) {
	dirs <- unlist(strsplit(path, split = "/"))
	r <- grep("\\.track$", unlist(dirs))
	if (length(r) > 0) {
		idx <- r[1]
		dirs[idx] <- paste(substr(dirs[idx], 0, nchar(dirs[idx]) - nchar(".tracks") + 1))
		return(paste(dirs[1 : idx], collapse = "."))
	}
	NULL				   
}

.gtrack.array.get_colnames <- function(trackstr) {
	.gcheckroot()

	if (is.na(match(trackstr, get("GTRACKS"))))
		stop(sprintf("Track %s does not exist", trackstr), call. = F)

	if (.gcall_noninteractive(gtrack.info, trackstr)$type != "array")
		stop("gtrack.array.get_colnames can only be applied to array tracks", call. = F)

	trackdir <- sprintf("%s.track", paste(get("GWD"), gsub("\\.", "/", trackstr), sep = "/"))
	filename <- paste(trackdir, ".colnames", sep = "/")

	if (!file.exists(filename))
		stop(sprintf("File %s does not exist", filename))

	f <- file(filename, "rb")
	colnames <- unserialize(f)
	close(f)
	colnames
}

.gtrack.array.set_colnames <- function(trackstr, names, check_num_cols) {
	.gcheckroot()

	if (is.na(match(trackstr, get("GTRACKS"))))
		stop(sprintf("Track %s does not exist", trackstr), call. = F)

	if (typeof(names) != "character")
		stop(sprintf("names parameter must be a character vector", trackstr), call. = F)

	if (.gcall_noninteractive(gtrack.info, trackstr)$type != "array")
		stop("gtrack.array.set_colnames can only be applied to array tracks", call. = F)

	if ("" %in% names)
		stop(sprintf("Column names cannot be empty", duplicated[1]), call. = F)

	duplicated <- names[duplicated(names)]
	if (length(duplicated)) 
		stop(sprintf("Column %s appears more than once", duplicated[1]), call. = F)

	if (check_num_cols) {
		oldnames <- .gtrack.array.get_colnames(trackstr)
		if (length(oldnames) != length(names))
			stop(sprintf("The number of columns in the track (%d) does not match the number of column names (%d)",
				length(oldnames), length(names)), call. = F)
	}

	colnames <- as.integer(1 : length(names))
	names(colnames) <- names

	trackdir <- sprintf("%s.track", paste(get("GWD"), gsub("\\.", "/", trackstr), sep = "/"))
	filename <- paste(trackdir, ".colnames", sep = "/")
	f <- file(filename, "wb")
	serialize(colnames, f)
	close(f)
}

.gtrack.attr.set <- function(trackstr, attr, value, overwrite.readonly) {
	table <- data.frame(value)
	colnames(table) <- attr
	rownames(table) <- trackstr

	.gtrack.attr.import(table, FALSE, overwrite.readonly)
}

.gtrack.attr.import <- function(table, remove.others, overwrite.readonly) {
	table[, 1 : ncol(table)] <- sapply(table[, 1 : ncol(table)], as.character)

	readonly_attrs <- NULL
	if (!overwrite.readonly)
	   readonly_attrs <- gdb.get_readonly_attrs()
	
	.gcall("gset_tracks_attrs", table, remove.others, readonly_attrs, new.env(parent = parent.frame()))
}

#' @export
.gtrack.prepare.pvals <- function(track) {
	.gcheckroot()
	
	trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())

	if (is.na(match(trackstr, get("GTRACKS"))))
		stop(sprintf("Track %s does not exist", trackstr), call. = F)

	quantile.edge.data.size <- .ggetOption("gquantile.edge.data.size", 100000)
	middle.size = .ggetOption("gpv.middle.size", 0.98)
	edge.size = (1 - middle.size) / 2
	middle.precision <- .ggetOption("gpv.middle.precision", 10^(-4))
	edge.precision <- .ggetOption("gpv.edge.precision", 10^(-9))

	# In the middle section the percentiles increase by a constant step = middle.precision.
	# At the edges the step is exponential. For instance, step i at the edge close to 0 follows the following function:
	#   k+i
	#  b     , where b and k are unknown. We also don't know n - the number of total steps required to cover the edge.
	#
	# b, k, n be calculated from the following:
	#
	# 1. b^(k+1) - b^k = edge.precision
	# 2. b^(k+n) - b^(k+n-1) = middle.precision
	# 3. b^(k+n) - b^k = edge.size

	b <- (edge.precision - edge.size) / (middle.precision - edge.size)
	k <- log(edge.precision / (b - 1), b)
	n <- ceiling(log(edge.size + b^k, b) - k)

	percentiles <- 0
	percentiles <- c(percentiles, b ^ (k + (0 : n)) - b^k)
	percentiles <- c(percentiles, 1 - percentiles)
	num.middle.steps <- middle.size / middle.precision
	percentiles <- c(percentiles, edge.size + (0 : num.middle.steps) * middle.precision)
	percentiles <- sort(percentiles)

	selected.percentiles <- NULL
	
	multitasking <- .ggetOption("gmultitasking")
	tryCatch({
		options(warn = -1) # disable warnings since gquantiles is going to warn about random sampling
		options(gmultitasking = FALSE)
		quantiles <- do.call(gquantiles, list(substitute(track), percentiles = c(0, percentiles)), envir = parent.frame())
		names(quantiles) <- NULL
		minval <- quantiles[1]
		maxval <- quantiles[length(quantiles)]
		quantiles <- quantiles[2 : length(quantiles)]
		
		# for each group of quantiles with identical value choose the maximal one
		selected.percentiles <- sapply(split(percentiles, quantiles), max)
		names(selected.percentiles) <- NULL

		# if all percentiles are equal create an artificial table
		if (length(selected.percentiles) == 1) {
			selected.percentiles <- c(1, 1)
			attr(selected.percentiles, "breaks") <- c(minval, maxval + 1)
		} else {
			indices <- match(selected.percentiles, percentiles)
			selected.quantiles <- quantiles[indices]
			attr(selected.percentiles, "breaks") <- selected.quantiles
		}
		
		attr(selected.percentiles, "minval") <- minval
		attr(selected.percentiles, "maxval") <- maxval
	},
			 finally = {
				options(warn = 0)
				options(gmultitasking = multitasking)
			}
	)

	# save the percentiles
	.gtrack.var.set(trackstr, "pv.percentiles", selected.percentiles)

	retv <- 0
}

.gtrack.var.exists <- function(trackname, varname) {
	trackdir <- sprintf("%s.track", paste(get("GWD"), gsub("\\.", "/", trackname), sep = "/"))
	filename <- paste(trackdir, "vars", varname, sep = "/")
	file.exists(filename)
}

.gtrack.var.get <- function(trackname, varname) {
	if (is.na(match(trackname, get("GTRACKS"))))
		stop(sprintf("Track %s does not exist", trackname), call. = F)

	trackdir <- sprintf("%s.track", paste(get("GWD"), gsub("\\.", "/", trackname), sep = "/"))
	filename <- paste(trackdir, "vars", varname, sep = "/")
	if (!file.exists(filename))
		stop(sprintf("Track variable %s does not exist", varname), call. = F)
	f <- file(filename, "rb")
	val <- unserialize(f)
	close(f)
	val
}

.gtrack.var.set <- function(trackname, varname, value) {
	if (is.na(match(trackname, get("GTRACKS"))))
		stop(sprintf("Track %s does not exist", trackname), call. = F)

	# if vars directory does not exist, create it
	trackdir <- sprintf("%s.track", paste(get("GWD"), gsub("\\.", "/", trackname), sep = "/"))
	dirname <- paste(trackdir, "vars", sep = "/")
	if (!file.exists(dirname))
		dir.create(dirname, mode = "0777")

	# save the variable
	filename <- paste(dirname, varname, sep = "/")
	f <- file(filename, "wb")
	serialize(value, f)
	close(f)
}

.gcluster.running.jobs <- function(jobids) {
	str <- system("qstat | sed 's/^[ ]*//' | cut -f 1 -d\" \"", intern=T)
	if (length(str) > 2)
		intersect(jobids, str)
	else
		c()
}

.gconfirmtrackcreate <- function(track) {
	if (!is.na(match(track, get("GTRACKS"))))
		stop(sprintf("Track %s already exists", track), call. = F)

	path <- gsub(".", "/", track, fixed = T)
	dir <- dirname(path)
	fulldir <- paste(get("GWD"), dir, sep = "/")
	fullpath <- sprintf("%s.track", paste(get("GWD"), path, sep = "/"))

	if (!file.exists(fulldir))
		stop(sprintf("Directory %s does not exist", dir), call. = F)

	if (file.exists(fullpath))
		stop(sprintf("File %s already exists", path), call. = F)

	if (!is.na(match(track, get("GINTERVS"))))
		stop(sprintf("Interval %s already exists", track), call. = F)

	if (!is.na(match(track, gvtrack.ls())))
		stop(sprintf("Virtual track %s already exists", track), call. = F)

	if (.ggetOption(".gautocompletion", FALSE) && exists(track))
		stop(sprintf("Variable \"%s\" shadows the name of the new track.\nPlease remove this variable from the environment or switch off autocompletion mode.", track), call. = F)
}

.gdb.add_track <- function(track) {
    .gcheckroot()
	
    trackdir <- sprintf("%s.track", paste(get("GWD"), gsub("\\.", "/", track), sep = "/"))
    if (file.exists(trackdir)) {
    	tracks <- sort(c(get("GTRACKS"), track))
    	intervals <- sort(get("GINTERVS"))
    	
    	res <- intersect(tracks, intervals)
    	if (length(res) > 0)
    		stop("The following tracks exist also as intervals: ", paste(res, collapse = " "))

    	if (.ggetOption(".gautocompletion", FALSE)) {
            if (exists(track, envir = .GlobalEnv))
                stop(sprintf("Variable \"%s\" shadows the name of identically named track.\nPlease remove this variable from the environment or switch off autocompletion mode.", track), call. = F)

        	if (.ggetOption(".ginteractive", FALSE))  # set track to NULL otherwise evaluation of track expression pmin(track, 2) will produce a string "2"
       			assign(track, NULL, envir = .GlobalEnv)
        	else
       			assign(track, track, envir = .GlobalEnv)
        }
    	
    	assign("GTRACKS", tracks, envir = .GlobalEnv)
    }
}

.gdb.rm_track <- function(track) {
    .gcheckroot()
	
    trackdir <- sprintf("%s.track", paste(get("GWD"), gsub("\\.", "/", track), sep = "/"))
    if (!file.exists(trackdir)) {
    	if (.ggetOption(".gautocompletion", FALSE)) {
            if (exists(track, envir = .GlobalEnv))
    		    remove(list = track, envir = .GlobalEnv)
        }

    	tracks <- get("GTRACKS")
        tracks <- tracks[tracks != track]	
    	assign("GTRACKS", tracks, envir = .GlobalEnv)
    }
}

.gdb.add_intervals.set <- function(intervals.set) {
    .gcheckroot()
	
    fname <- sprintf("%s.interv", paste(get("GWD"), gsub("\\.", "/", intervals.set), sep = "/"))
    if (file.exists(fname)) {
    	tracks <- get("GTRACKS")
    	intervals <- sort(c(get("GINTERVS"), intervals.set))
    	
    	res <- intersect(tracks, intervals)
    	if (length(res) > 0)
    		stop("The following tracks exist also as intervals: ", paste(res, collapse = " "))

    	if (.ggetOption(".gautocompletion", FALSE)) {
            if (exists(intervals.set, envir = .GlobalEnv))
                stop(sprintf("Variable \"%s\" shadows the name of identically named intervals set.\nPlease remove this variable from the environment or switch off autocompletion mode.", intervals.set), call. = F)

    	    assign(intervals.set, intervals.set, envir = .GlobalEnv)
        }
    	
    	assign("GINTERVS", intervals, envir = .GlobalEnv)
    }
}

.gdb.rm_intervals.set <- function(intervals.set) {
    .gcheckroot()
	
    fname <- sprintf("%s.interv", paste(get("GWD"), gsub("\\.", "/", intervals.set), sep = "/"))
    if (!file.exists(fname)) {
    	if (.ggetOption(".gautocompletion", FALSE)) {
            if (exists(intervals.set, envir = .GlobalEnv))
    		    remove(list = intervals.set, envir = .GlobalEnv)
        }

    	intervals <- get("GINTERVS")
        intervals <- intervals[intervals != intervals.set]
    	assign("GINTERVS", intervals, envir = .GlobalEnv)
    }
}

.gseq.import <- function(groot = NULL, path = NULL) {
	chroms <- c()
	files <- c()
	tmp.dirname <- tempfile(pattern = "", tmpdir = paste(groot, "/downloads", sep = ""))
	if (!dir.create(tmp.dirname, recursive = T, mode = "0777"))
		stop(sprintf("Failed to create a directory %s", tmp.dirname), call. = F)

	tryCatch({
		files <- c()

		for (fasta in path) {
        	protocol <- "ftp://"
			.files <- c()
        	if (substr(fasta, 1, nchar(protocol)) == protocol) {
        		cat("Downloading FASTA files...\n")

        	    # ftp
        		.files <- gwget(fasta, tmp.dirname)
        		.files <- paste(tmp.dirname, "/", .files, sep = "")
        	} else
            	# local path
        		.files <- system(paste("/bin/sh -c \"ls -d -A", fasta, "\""), intern = T)

			files <- c(files, .files)
		}

    	cat("Building Seq files...\n")
    	fastas <- files[grep("^chr.+$", basename(files), perl = T)]
    	for (fasta in fastas) {
    		chrom <- gsub("^chr(\\w+)(\\..*)*$", "\\1", basename(fasta), perl = T)
			if (!is.na(match(chrom, chroms)))
				next

    		fasta.original <- fasta
    				
    		if (length(grep("^.+\\.gz$", fasta, perl = T))) {
    			.chrom <- basename(gsub("^(.+)\\.gz$", "\\1", fasta, perl = T))
    			fasta.unzipped <- paste(tmp.dirname, "/", .chrom, sep = "")
    			cmd <- paste("/bin/sh -c \"gunzip -q -c", fasta, ">", fasta.unzipped, "\"")
    			if (system(cmd))
    				stop(sprintf("Command failed: %s", cmd), call. = F)
    			fasta <- fasta.unzipped
    		}

    		seq <- sprintf("%s/seq/chr%s.seq", groot, chrom)

    		cat(sprintf("chr%s\n", chrom))
			.gcall("gseqimport", fasta, seq, new.env(parent = parent.frame()), silent = TRUE)

    		chroms <- c(chroms, chrom)
    	}
	}, finally = {
		unlink(tmp.dirname, recursive = TRUE) }
	)
	chroms
}

.gslice <- function(trackstr, slice) {
    res <- list()
	colnames <- .gtrack.array.get_colnames(trackstr)

    if (is.null(slice)) {
        res$slice <- NULL
        res$colnames <- names(colnames)
    } else if (typeof(slice) == "character") {
        slice <- unique(slice)
		res$slice <- colnames[slice]
        res$colnames <- names(res$slice)

		idx <- match(NA, res$slice)
		if (!is.na(idx))
		    stop(sprintf("%s does not appear among the column names of track %s", slice[idx], trackstr), call. = F)
    } else if (is.numeric(slice) || is.integer(slice)) {
        if (TRUE %in% (as.integer(slice) != slice))
            stop("Invalid type of slice parameter", call. = F)

        slice <- unique(slice)
		slice <- as.integer(slice)

        outofrange <- slice < 1 | slice > length(colnames)
        if (TRUE %in% outofrange)
		    stop(sprintf("Slice index %d is out of range", slice[match(TRUE, outofrange)], trackstr), call. = F)

        res$slice <- colnames[slice]
        res$colnames <- names(res$slice)
    } else
        stop("Invalid type of slice parameter", call. = F)

    res
}

.gtrack.create_test_arrays <- function(track, minsize, maxsize, intervals = get("ALLGENOME"), iterator = NULL) {
	if (is.null(substitute(track)))
		stop("Usage: .gtrack.create_test_arrays(track, expr, iterator = NULL, band = NULL)", call. = F)
	.gcheckroot()

	trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
	.iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())
	trackdir <- sprintf("%s.track", paste(get("GWD"), gsub("\\.", "/", trackstr), sep = "/"))

	direxisted <- file.exists(trackdir)

	if (!is.na(match(trackstr, get("GTRACKS"))))
		stop(sprintf("Track %s already exists", trackstr), call. = F)
	
	.gconfirmtrackcreate(trackstr)
	success <- FALSE
	tryCatch({
		colnames <- .gcall("_gcreate_arrays_track", trackstr, minsize, maxsize, "1", intervals, .iterator, new.env(parent = parent.frame()), silent = TRUE)
        .gdb.add_track(trackstr)
		.gtrack.array.set_colnames(trackstr, colnames, FALSE)
		.gtrack.attr.set(trackstr, "created.by", ".gtrack.create_test_arrays", T)
		.gtrack.attr.set(trackstr, "created.date", date(), T)
		success <- TRUE
	},
			 finally = {
				 if (!success && !direxisted) {
					 unlink(trackdir, recursive = TRUE)
					 .gdb.rm_track(trackstr)
				 }
			 }
	)
	retv <- 0 # suppress return value
}

.gtrack.create_test_computer2d <- function(track = NULL, prob.skip.chrom = 0.5, max.rect = 100000, max.rect.size = 10000) {
	if (is.null(substitute(track)))
		stop("Usage: gtrack.create_test_computer2d(trackname, prob.skip.chrom, max.rect, max.rect.size)", call. = F)
	.gcheckroot()

	trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
	trackdir <- sprintf("%s.track", paste(get("GWD"), gsub("\\.", "/", trackstr), sep = "/"))

	direxisted <- file.exists(trackdir)

	if (!is.na(match(trackstr, get("GTRACKS"))))
		stop(sprintf("Track %s already exists", trackstr), call. = F)
	
	.gconfirmtrackcreate(trackstr)
	success <- FALSE
	tryCatch({
		.gcall("gcreate_test_computer2d_track", trackstr, prob.skip.chrom, max.rect, max.rect.size, new.env(parent = parent.frame()), silent = TRUE)
        .gdb.add_track(trackstr)
		.gtrack.attr.set(trackstr, "created.by", sprintf(".gtrack.create_test_computer2d(%s, %g, %g, %g)", trackstr, prob.skip.chrom, max.rect, max.rect.size), T)
		.gtrack.attr.set(trackstr, "created.date", date(), T)
		success <- TRUE
	},
			 finally = {
				 if (!success && !direxisted) {
					 unlink(trackdir, recursive = TRUE)
					 .gdb.rm_track(trackstr)
				 }
			 }
	)
	retv <- 0 # suppress return value
}

.gdefine_autocompletion_vars <- function(tracks, intervals, vtracks, interactive) {
	for (track in tracks) {
		if (exists(track))
			stop(sprintf("Variable \"%s\" shadows the name of identically named track.\nPlease remove this variable from the environment or switch off autocompletion mode.", track), call. = F)
	}

	for (interval in intervals) {
		if (exists(interval))
			stop(sprintf("Variable \"%s\" shadows the name of identically named intervals set.\nPlease remove this variable from the environment or switch off autocompletion mode.", interval), call. = F)
	}
	
	for (vtrack in vtracks) {
		if (exists(vtrack))
			stop(sprintf("Variable \"%s\" shadows the name of identically named virtual track.\nPlease remove this variable from the environment or switch off autocompletion mode.", vtrack), call. = F)
	}

	if (interactive) {
		# set track to NULL otherwise evaluation of track expression pmin(track, 2) will produce a string "2"
		for (track in tracks)
			assign(track, NULL, envir = .GlobalEnv)

		for (vtrack in vtracks)
			assign(vtrack, NULL, envir = .GlobalEnv)
	} else {
		for (track in tracks)
			assign(track, track, envir = .GlobalEnv)

		for (vtrack in vtracks)
			assign(vtrack, vtrack, envir = .GlobalEnv)
	}

	for (interval in intervals)
		assign(interval, interval, envir = .GlobalEnv)
}

.gundefine_autocompletion_vars <- function() {
	options(warn = -1) # disable warnings: some variables might be removed already by the user

	if (exists("GTRACKS"))
		remove(list = get("GTRACKS"), envir = .GlobalEnv)

	if (exists("GINTERVS"))
		remove(list = get("GINTERVS"), envir = .GlobalEnv)


	vtracks <- gvtrack.ls()
	if (!is.null(vtracks))
		remove(list = vtracks, envir = .GlobalEnv)

	options(warn = -1) # restore warnings
}



#' Computes auto-correlation between the strands for a file of mapped sequences
#' 
#' Calculates auto-correlation between plus and minus strands for the given
#' chromosome in a file of mapped sequences.
#' 
#' This function calculates auto-correlation between plus and minus strands for
#' the given chromosome in a file of mapped sequences. Each line in the file
#' describes one read. Each column is separated by a TAB character.
#' 
#' The following columns must be presented in the file: sequence, chromosome,
#' coordinate and strand. The position of these columns are controlled by
#' 'cols.order' argument accordingly. The default value of 'cols.order' is a
#' vector (9,11,13,14) meaning that sequence is expected to be found at column
#' number 9, chromosome - at column 11, coordinate - at column 13 and strand -
#' at column 14. The first column should be referenced by 1 and not by 0.
#' 
#' Coordinates that are not in [min.coord, max.coord] range are ignored.
#' 
#' gcompute_strands_autocorr outputs the total statistics and the
#' auto-correlation given by bins. The size of the bin is indicated by
#' 'binsize' parameter. Statistics is calculated for bins in the range of
#' [-maxread, maxread].
#' 
#' @param file the name of the file containing mapped sequences
#' @param chrom chromosome for which the auto-correlation is computed
#' @param binsize calculate the auto-correlation for bins in the range of
#' [-maxread, maxread]
#' @param maxread maximal length of the sequence used for statistics
#' @param cols.order order of sequence, chromosome, coordinate and strand
#' columns in file
#' @param min.coord minimal coordinate used for statistics
#' @param max.coord maximal coordinate used for statistics
#' @return Statistics for each strand and auto-correlation by given bins.
#' @keywords ~gcompute_strands_autocorr ~auto-correlation ~autocorrelation
#' ~correlation
#' @examples
#' 
#' gdb.init_examples()
#' gcompute_strands_autocorr(paste(GROOT, "reads", sep = "/"),
#'                           "chr1", 50, maxread = 300)
#' 
#' @export gcompute_strands_autocorr
gcompute_strands_autocorr <- function(file = NULL, chrom = NULL, binsize = NULL, maxread = 400, cols.order = c(9, 11, 13, 14), min.coord = 0, max.coord = 3e+8) {
	if (is.null(file) || is.null(chrom) || is.null(binsize))
		stop("Usage: gcompute_strands_autocorr(file, chrom, binsize, maxread = 400, cols.order = c(9, 11, 13, 14), min.coord = 0, max.coord = 3e+8)", call. = F)
	.gcheckroot()

	if (substr(chrom, 1, 3) != "chr")
		chrom <- paste("chr", chrom, sep = "")
	
	res <- .gcall("gcompute_strands_autocorr", file, chrom, binsize, maxread, cols.order, min.coord, max.coord, new.env(parent = parent.frame()))
	res
}



#' Creates a new Genomic Database
#' 
#' Creates a new Genomic Database.
#' 
#' This function creates a new Genomic Database at the location specified by
#' 'groot'. FASTA files are converted to 'Seq' format and appropriate
#' 'chrom_sizes.txt' file is generated (see "User Manual" for more details).
#' 
#' If 'genes.file' is not 'NULL' four sets of intervals are created in the
#' database: \code{tss}, \code{exons}, \code{utr3} and \code{utr5}. See
#' \link{gintervals.import_genes} for more details about importing genes
#' intervals.
#' 
#' 'fasta', 'genes.file' and 'annots.file' can be either a file path or URL in
#' a form of 'ftp://[address]/[file]'. 'fasta' can also contain wildcards to
#' indicate multiple files. Files that these arguments point to can be zipped
#' or unzipped.
#' 
#' @param groot path to newly created database
#' @param fasta an array of names or URLs of FASTA files. Can contain wildcards
#' for multiple files
#' @param genes.file name or URL of file that contains genes. If 'NULL' no
#' genes are imported
#' @param annots.file name of URL file that contains annotations. If 'NULL' no
#' annotations are imported
#' @param annots.names annotations names
#' @return None.
#' @seealso \code{\link{gdb.init}}, \code{\link{gdb.reload}},
#' \code{\link{gintervals.import_genes}}
#' @keywords ~database ~create ~genes
#' @examples
#' 
#' \dontrun{
#' ftp = "ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19"
#' gdb.create("mydb",
#'            c(paste(ftp, "chromosomes/chr2*", sep = "/"),
#'              paste(ftp, "chromosomes/chr7.*", sep = "/")),
#'            paste(ftp, "database/knownGene.txt.gz", sep = "/"),
#'            paste(ftp, "database/kgXref.txt.gz", sep = "/"),
#'            c("kgID", "mRNA", "spID", "spDisplayID", "geneSymbol",
#'              "refseq", "protAcc", "description", "rfamAcc",
#'              "tRnaName"))
#' gdb.init("mydb")
#' gintervals.ls()
#' gintervals.all()
#' }
#' 
#' @export gdb.create
gdb.create <- function(groot = NULL, fasta = NULL, genes.file = NULL, annots.file = NULL, annots.names = NULL) {
	if (is.null(groot) || is.null(fasta))
		stop("Usage: gdb.create(groot, fasta, genes.file = NULL, annots.file = NULL, annots.names = NULL)", call. = F);

	if (file.exists(groot))
		stop(sprintf("Directory %s already exists", groot), call. = F)

	success <- FALSE
	allgenome.old <- get("ALLGENOME")

	tryCatch({
		dir.create(groot, showWarnings = F, recursive = TRUE, mode = "0777")
		dir.create(paste(groot, "pssms", sep = "/"), showWarnings = F, recursive = TRUE, mode = "0777")
		dir.create(paste(groot, "seq", sep = "/"), showWarnings = F, recursive = TRUE, mode = "0777")
		dir.create(paste(groot, "tracks", sep = "/"), showWarnings = F, recursive = TRUE, mode = "0777")

		chroms = .gseq.import(groot, fasta)

		if (!length(chroms))
			stop("No FASTA files were imported", call. = F)

		seq.files <- paste("chr", chroms, ".seq", sep = "")
		seq.files <- paste(paste(groot, "seq", sep = "/"), seq.files, sep = "/")
		chrom.sizes <- data.frame(chrom = chroms, size = file.info(seq.files)$size)
		write.table(chrom.sizes, paste(groot, "chrom_sizes.txt", sep = "/"), quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)

		# before calling gintervals.import_genes new ALLGENOME must be set
    	intervals <- data.frame(chrom = as.factor(paste("chr", as.character(chrom.sizes$chrom), sep = "")),
    							start = 0, end = as.numeric(chrom.sizes$size))
    	intervals <- intervals[order(intervals$chrom), ]
    	rownames(intervals) <- 1 : nrow(intervals)

    	cartesian <- expand.grid(1 : nrow(intervals), 1 : nrow(intervals))
    	intervals2d <- cbind(intervals[cartesian[, 2], ], intervals[cartesian[, 1], ])
    	names(intervals2d) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2")
    	rownames(intervals2d) <- 1 : nrow(intervals2d)

    	assign("ALLGENOME", list(intervals, intervals2d), envir = .GlobalEnv)

		if (!is.null(genes.file)) {
			intervs <- gintervals.import_genes(genes.file, annots.file, annots.names)
			if (!is.null(intervs)) {
				for (i in 1 : length(intervs)) {
					if (!is.null(intervs$tss))
						.gcall_noninteractive(.gintervals.save_file, sprintf("%s/tracks/%s.interv", groot, names(intervs)[i]), intervs[[i]])
				}
			}
		}

		# write read-only attributes
		f <- file(paste(groot, ".ro_attributes", sep = "/"), "wb")
		serialize(c("created.by", "created.date"), f)
		close(f)

		cat("Database was successfully created\n")
		success <- TRUE
	}, finally = {
		assign("ALLGENOME", allgenome.old, envir = .GlobalEnv)
		if (!success)
			unlink(groot, recursive = TRUE)
	})
	retv <- 0 # suppress return value
}



#' Returns a list of read-only track attributes
#' 
#' Returns a list of read-only track attributes.
#' 
#' This function returns a list of read-only track attributes. These attributes
#' are not allowed to be modified or deleted.
#' 
#' If no attributes are marked as read-only a 'NULL' is returned.
#' 
#' @return A list of read-only track attributes.
#' @seealso \code{\link{gdb.set_readonly_attrs}},
#' \code{\link{gtrack.attr.get}}, \code{\link{gtrack.attr.set}}
#' @keywords ~attr ~attribute
#' @export gdb.get_readonly_attrs
gdb.get_readonly_attrs <- function() {
	.gcheckroot()

	filename <- paste(get("GROOT"), ".ro_attributes", sep = "/")
	attrs <- NULL
	if (file.exists(filename)) {
		f <- file(filename, "rb")
		attrs <- unserialize(f)
		close(f)
		if (!is.character(attrs))
		   stop(sprintf("Invalid format of read-only atrributes file %s", filename), call. = F)

		attrs <- unique(attrs)
		attrs <- attrs[attrs != ""]
	}
	attrs
}



#' Initializes connection with Genomic Database
#' 
#' Initializes connection with Genomic Database: loads the list of tracks,
#' intervals, etc.
#' 
#' 'gdb.init' initializes the connection with the Genomic Database. It is
#' typically called first prior to any other function. When the package is
#' attached it internally calls to 'gdb.init.examples' which opens the
#' connection with the database located at 'PKGDIR/trackdb/test' directory,
#' where 'PKGDIR' is the directory where the package is installed.
#' 
#' The current working directory inside the Genomic Database is set to 'dir'.
#' If 'dir' is 'NULL', the current working directory is set to 'GROOT/tracks'.
#' 
#' If 'rescan' is 'TRUE', the list of tracks and intervals is achieved by
#' rescanning directory structure under the current current working directory.
#' Otherwise 'gdb.init' attempts to use the cached list that resides in
#' 'groot/.db.cache' file.
#' 
#' Upon completion the connection is established with the database. If
#' auto-completion mode is switched on (see 'gset_input_method') the list of
#' tracks and intervals sets is loaded and added as variables to the global
#' environment allowing auto-completion of object names with <TAB> key. Also a
#' few global variables are defined. These variables should not be modified by
#' user.
#' 
#' \tabular{ll}{ GROOT \tab Root directory of Genomic Database\cr GWD \tab
#' Current working directory inside Genomic Database\cr GTRACKS \tab List of
#' all available tracks\cr GINTERVS \tab List of all available intervals\cr
#' GVTRACKS \tab List of all available virtual tracks\cr ALLGENOME \tab List of
#' all chromosomes and their sizes\cr GITERATOR.INTERVALS \tab A set of
#' iterator intervals for which the track expression is evaluated\cr }
#' 
#' @aliases gdb.init gdb.init.examples
#' @param groot the root directory of the Genomic Database
#' @param dir the current working directory inside the Genomic Database
#' @param rescan indicates whether the file structure should be rescanned
#' @return None.
#' @seealso \code{\link{gdb.reload}}, \code{\link{gdb.create}},
#' \code{\link{gdir.cd}}, \code{\link{gtrack.ls}}, \code{\link{gintervals.ls}},
#' \code{\link{gvtrack.ls}}, \code{\link{gset_input_mode}}
#' @keywords ~db ~data ~database
#' @export gdb.init
gdb.init <- function(groot = NULL, dir = NULL, rescan = FALSE) {
	gsetroot(groot, dir, rescan)
}

gdb.init_examples <- function() {
	gsetroot(paste(.GLIBDIR, "trackdb/test", sep = "/"))
}



#' Reloads database from the disk
#' 
#' Reloads database from disk: list of tracks, intervals, etc.
#' 
#' Reloads Genomic Database from disk: list of tracks, intervals, etc. Use this
#' function if you manually add tracks or if for any reason the database
#' becomes corrupted. If 'rescan' is 'TRUE', the list of tracks and intervals
#' is achieved by rescanning directory structure under the current current
#' working directory. Otherwise 'gdb.reload' attempts to use the cached list
#' that resides in 'GROOT/.db.cache' file.
#' 
#' @param rescan indicates whether the file structure should be rescanned
#' @seealso \code{\link{gdb.init}}, \code{\link{gdb.create}},
#' \code{\link{gdir.cd}}, \code{\link{gset_input_mode}}
#' @keywords ~db
#' @export gdb.reload
gdb.reload <- function(rescan = TRUE) {
	if (!exists("GROOT"))
		stop("gdb.init() must be called beforehand.", call. = F)
	
	if (.ggetOption(".gautocompletion", FALSE))
		.gundefine_autocompletion_vars()

	assign("GTRACKS", NULL, envir = .GlobalEnv)
	assign("GINTERVS", NULL, envir = .GlobalEnv)

	dir <- get("GWD")

	res <- ""

	if (get("GWD") != paste(get("GROOT"), "tracks", sep = "/"))
		rescan = TRUE

	db.filename <- paste(get("GROOT"), ".db.cache", sep = "/")
	
	options(warn = -1) # disable warnings since dir() on non dir or non existing dir produces warnings
	if (!rescan) {
		retv <- try({
			f <- file(db.filename, "rb")
			res <- unserialize(f)
			close(f)
		}, silent = T)

		if (inherits(retv, "try-error"))
			rescan = TRUE
	}

	if (rescan) {
		res <- .gcall("gfind_tracks_n_intervals", dir, new.env(parent = parent.frame()), silent = TRUE)
		if (get("GWD") == paste(get("GROOT"), "tracks", sep = "/"))
			try({
				f <- file(db.filename, "wb")
				serialize(res, f)
				close(f)
			}, silent = T)
		else
			unlink(db.filename, recursive = TRUE)
	}
	options(warn = 0) # restore the warning behavior
	
	tracks <- res[[1]]
	intervals <- res[[2]]
	
	tracks <- sort(tracks)
	intervals <- sort(intervals)
	
	res <- intersect(tracks, intervals)
	if (length(res) > 0)
		stop("The following tracks exist also as intervals: ", paste(res, collapse = " "))

	if (.ggetOption(".gautocompletion", FALSE))
		.gdefine_autocompletion_vars(tracks, intervals, gvtrack.ls(), .ggetOption(".ginteractive", FALSE))
	
	assign("GTRACKS", tracks, envir = .GlobalEnv)
	assign("GINTERVS", intervals, envir = .GlobalEnv)
}



#' Sets read-only track attributes
#' 
#' Sets read-only track attributes.
#' 
#' This function sets the list of read-only track attributes. The specified
#' attributes may or may not already exist in the tracks.
#' 
#' If 'attrs' is 'NULL' the list of read-only attributes is emptied.
#' 
#' @param attrs a vector of read-only attributes names or 'NULL'
#' @return None.
#' @seealso \code{\link{gdb.get_readonly_attrs}},
#' \code{\link{gtrack.attr.get}}, \code{\link{gtrack.attr.set}}
#' @keywords ~attr ~attribute
#' @export gdb.set_readonly_attrs
gdb.set_readonly_attrs <- function(attrs) {
	.gcheckroot()

	filename <- paste(get("GROOT"), ".ro_attributes", sep = "/")
	
	if (is.null(attrs))
	   unlink(filename)
	else {
	    attrs <- as.character(attrs)

		idx <- which(duplicated(attrs))[1]
		if (!is.na(idx))
	       stop(sprintf("Attribute %s appears more than once", attrs[idx]), call. = F)
		
		idx <- which(attrs == "")[1]
		if (!is.na(idx))
	       stop("Attribute name cannot be an empty string", call. = F)
		
		f <- file(filename, "wb")
		serialize(attrs, f)
		close(f)
	}
	retv <- 0 # suppress return value
}



#' Calculates distribution of contact distances
#' 
#' Calculates distribution of contact distances.
#' 
#' A 2D iterator interval '(chrom1, start1, end1, chrom2, start2, end2)' is
#' said to represent a contact between two 1D intervals I1 and I2: '(chrom1,
#' start1, end1)' and '(chrom2, start2, end2)'.
#' 
#' For contacts where 'chrom1' equals to 'chrom2' and I1 is within source
#' intervals the function calculates the distribution of distances between I1
#' and I2. The distribution is calculated separately for intra-domain and
#' inter-domain contacts.
#' 
#' An interval is within source intervals if the unification of all source
#' intervals fully overlaps it. 'src' intervals are allowed to contain
#' overlapping intervals.
#' 
#' Two intervals I1 and I2 are within the same domain (intra-domain contact) if
#' among the domain intervals exists an interval that fully overlaps both I1
#' and I2. Otherwise the contact is considered to be inter-domain. 'domain'
#' must contain only non-overlapping intervals.
#' 
#' The distance between I1 and I2 is the absolute distance between the centers
#' of these intervals, i.e.: '|(start1 + end1 - start2 - end2) / 2|'.
#' 
#' The range of distances for which the distribution is calculated is defined
#' by 'breaks' argument. For example: 'breaks=c(x1, x2, x3, x4)' represents
#' three different intervals (bins): (x1, x2], (x2, x3], (x3, x4].
#' 
#' If 'include.lowest' is 'TRUE' the the lowest value will be included in the
#' first interval, i.e. in [x1, x2]
#' 
#' @param expr track expression
#' @param breaks breaks that determine the bin
#' @param src source intervals
#' @param domain domain intervals
#' @param intervals genomic scope for which the function is applied
#' @param include.lowest if 'TRUE', the lowest value of the range determined by
#' breaks is included
#' @param iterator 2D track expression iterator. If 'NULL' iterator is
#' determined implicitly based on track expressions.
#' @param band track expression band. If 'NULL' no band is used.
#' @return 2-dimensional vector representing the distribution of contact
#' distances for inter and intra domains.
#' @seealso \code{\link{gdist}}, \code{\link{gtrack.2d.import_contacts}}
#' @keywords ~contacts
#' @examples
#' 
#' gdb.init_examples()
#' 
#' src <- rbind(gintervals(1, 10, 100),
#'              gintervals(1, 200, 300),
#'              gintervals(1, 400, 500),
#'              gintervals(1, 600, 700),
#'              gintervals(1, 7000, 9100),
#'              gintervals(1, 9000, 18000),
#'              gintervals(1, 30000, 31000),
#'              gintervals(2, 1130, 15000))
#'              
#' domain <- rbind(gintervals(1, 0, 483000),
#'                 gintervals(2, 0, 300000))
#' 
#' gcis_decay("rects_track", 50000 * (1 : 10), src, domain)
#' 
#' @export gcis_decay
gcis_decay <- function(expr = NULL, breaks = NULL, src = NULL, domain = NULL, intervals = NULL, include.lowest = FALSE, iterator = NULL, band = NULL) {
	if (is.null(substitute(expr)) || is.null(breaks) || is.null(src) || is.null(domain))
		stop("Usage: gcis_decay(expr, breaks, src, domain, intervals = ALLGENOME, include.lowest = FALSE, iterator = NULL, band = NULL)", call. = F)
	.gcheckroot()

	if (is.null(intervals))
		intervals <- get("ALLGENOME")
	
	exprstr <- do.call(.gexpr2str, list(substitute(expr)), envir = parent.frame())
	.iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())

	res <- .gcall("gcis_decay", exprstr, breaks, src, domain, intervals, include.lowest, .iterator, band, new.env(parent = parent.frame()))
	attr(res, "breaks") <- breaks
	res
}



#' Calculates distribution of track expressions
#' 
#' Calculates distribution of track expressions' values over the given set of
#' bins.
#' 
#' This function calculates the distribution of values of the numeric track
#' expressions over the given set of bins.
#' 
#' The range of bins is determined by 'breaks' argument. For example:
#' 'breaks=c(x1, x2, x3, x4)' represents three different intervals (bins): (x1,
#' x2], (x2, x3], (x3, x4].
#' 
#' If 'include.lowest' is 'TRUE' the the lowest value will be included in the
#' first interval, i.e. in [x1, x2]
#' 
#' 'gdist' can work with any number of dimensions. If more than one
#' 'expr'-'breaks' pair is passed, the result is a multidimensional vector, and
#' an individual value can be accessed by [i1,i2,...,iN] notation, where 'i1'
#' is the first track and 'iN' is the last track expression.
#' 
#' @param expr track expression
#' @param breaks breaks that determine the bin
#' @param intervals genomic scope for which the function is applied
#' @param include.lowest if 'TRUE', the lowest value of the range determined by
#' breaks is included
#' @param iterator track expression iterator. If 'NULL' iterator is determined
#' implicitly based on track expressions.
#' @param band track expression band. If 'NULL' no band is used.
#' @return N-dimensional vector where N is the number of 'expr'-'breaks' pairs.
#' @seealso \code{\link{gextract}}
#' @keywords ~distribution
#' @examples
#' 
#' gdb.init_examples()
#' 
#' ## calculate the distribution of dense_track for bins:
#' ## (0, 0.2], (0.2, 0.5] and (0.5, 1]
#' gdist("dense_track", c(0, 0.2, 0.5, 1))
#' 
#' ## calculate two-dimensional distribution:
#' ## dense_track vs. sparse_track
#' gdist("dense_track", seq(0, 1, by = 0.1), "sparse_track",
#'       seq(0, 2, by = 0.2), iterator = 100)
#' 
#' @export gdist
gdist <- function(..., intervals = NULL, include.lowest = FALSE, iterator = NULL, band = NULL) {
	args <- as.list(substitute(list(...)))[-1L]
	if (length(args) < 2 || (length(args) %% 2 != 0 && (length(args) - 1) %% 2 != 0))
		stop("Usage: gdist([expr, breaks]+, intervals = ALLGENOME, include.lowest = FALSE, iterator = NULL, band = NULL)", call. = F)
	.gcheckroot()

	if (length(args) %% 2 != 0)
		intervals <- eval.parent(args[[length(args)]])
	else if (is.null(intervals))
		intervals <- get("ALLGENOME")

	exprs <- c()
	breaks <- list()
	
	for (i in (0 : (length(args) / 2 - 1))) {
		exprs <- append(exprs, do.call(.gexpr2str, list(args[[i * 2 + 1]]), envir = parent.frame()))
		breaks[[length(breaks) + 1]] <- eval.parent(args[[i * 2 + 2]])
	}

	.iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())

	if (.ggetOption("gmultitasking"))
		res <- .gcall("gtrackdist_multitask", intervals, exprs, breaks, include.lowest, .iterator, band, new.env(parent = parent.frame()))
	else
		res <- .gcall("gtrackdist", intervals, exprs, breaks, include.lowest, .iterator, band, new.env(parent = parent.frame()))
	attr(res, "breaks") <- breaks
	res
}



#' Returns evaluated track expression
#' 
#' Returns the result of track expressions evaluation for each of the iterator
#' intervals.
#' 
#' This function returns the result of track expressions evaluation for each of
#' the iterator intervals. The returned value is a set of intervals with an
#' additional column for each of the track expressions. This value can be used
#' as an input for any other function that accepts intervals. If the intervals
#' inside 'intervals' argument overlap gextract returns the overlapped
#' coordinate more than once.
#' 
#' The order inside the result might not be the same as the order of intervals.
#' An additional column 'intervalID' is added to the return value. Use this
#' column to refer to the index of the original interval from the supplied
#' 'intervals'.
#' 
#' If 'file' parameter is not 'NULL' the result is outputed to a tab-delimited
#' text file (without 'intervalID' column) rather than returned to the user.
#' This can be especially useful when the result is too big to fit into the
#' physical memory.  The resulted file can be used as an input for
#' 'gtrack.import' or 'gtrack.array.import' functions.
#' 
#' If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
#' set. Similarly to 'file' parameter 'intervals.set.out' can be useful to
#' overcome the limits of the physical memory.
#' 
#' 'colnames' parameter controls the names of the columns that contain the
#' evaluated expressions. By default the column names match the track
#' expressions.
#' 
#' @param expr track expression
#' @param intervals genomic scope for which the function is applied
#' @param colnames sets the columns names in the returned value. If 'NULL'
#' names are set to track expression.
#' @param iterator track expression iterator. If 'NULL' iterator is determined
#' implicitly based on track expressions.
#' @param band track expression band. If 'NULL' no band is used.
#' @param file file name where the function result is optionally outputed in
#' tab-delimited format
#' @param intervals.set.out intervals set name where the function result is
#' optionally outputed
#' @return If 'file' and 'intervals.set.out' are 'NULL' a set of intervals with
#' an additional column for each of the track expressions and 'columnID'
#' column.
#' @seealso \code{\link{gtrack.array.extract}}, \code{\link{gsample}},
#' \code{\link{gtrack.import}}, \code{\link{gtrack.array.import}},
#' \code{\link{glookup}}, \code{\link{gpartition}}, \code{\link{gdist}}
#' @keywords ~extract
#' @examples
#' 
#' gdb.init_examples()
#' 
#' ## get values of 'dense_track' for [0, 500), chrom 1
#' gextract("dense_track", gintervals(1, 0, 500))
#' 
#' ## get values of 'rects_track' (a 2D track) for a 2D interval
#' gextract("rects_track",
#'          gintervals.2d("chr1", 0, 4000, "chr2", 2000, 5000))
#' 
#' ## get values of two track expressions 'dense_track' and
#' ## 'array_track * 2' running over '100' iterator
#' gextract("dense_track", "array_track * 2", gintervals(1, 0, 500),
#'          iterator = 100, colnames = c("expr1", "expr2"))
#' 
#' @export gextract
gextract <- function(..., intervals = NULL, colnames = NULL, iterator = NULL, band = NULL, file = NULL, intervals.set.out = NULL) {
	args <- as.list(substitute(list(...)))[-1L]
	if (is.null(intervals) && length(args) < 2 || !is.null(intervals) && length(args) < 1)
		stop("Usage: gextract([expr]+, intervals, colnames = NULL, iterator = NULL, band = NULL, file = NULL, intervals.set.out = NULL)", call. = F)
	.gcheckroot()

	if (is.null(intervals)) {
		intervals <- eval.parent(args[[length(args)]])
		args <- args[1 : (length(args) - 1)]
	}

	tracks <- c()
	for (track in args)
		tracks <- c(tracks, do.call(.gexpr2str, list(track), envir = parent.frame()))

	.iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())
	intervals.set.out <- do.call(.gexpr2str, list(substitute(intervals.set.out)), envir = parent.frame())

	if (!is.null(intervals.set.out))
		fullpath <- .gintervals.check_new_set(intervals.set.out)

	# intervals can be NULL if gextract is piped with gscreen and the latter returns NULL
	success <- FALSE
	res <- NULL
	tryCatch({
		if (!is.null(intervals)) {
			if (.ggetOption("gmultitasking"))
				res <- .gcall("gextract_multitask", intervals, tracks, colnames, .iterator, band, file, intervals.set.out, new.env(parent = parent.frame()))
			else
				res <- .gcall("gextract", intervals, tracks, colnames, .iterator, band, file, intervals.set.out, new.env(parent = parent.frame()))

			if (!is.null(intervals.set.out) && .gintervals.is_bigset(intervals.set.out, F) && !.gintervals.needs_bigset(intervals.set.out))
				.gintervals.big2small(intervals.set.out)
		}

		success <- TRUE
	},
	finally = {
		if (!success && !is.null(intervals.set.out))
			unlink(fullpath, recursive = TRUE)
	})

	# refresh the list of GINTERVS, etc.
	if (!is.null(intervals.set.out)) {
		.gdb.add_intervals.set(intervals.set.out)
		retv <- 0 # suppress return value
	} else if (!is.null(file))
		retv <- 0 # suppress return value
	else
		res
}



#' Partitions the values of track expression
#' 
#' Converts the values of track expression to intervals that match
#' corresponding bin.
#' 
#' This function converts first the values of track expression into 1-based
#' bin's index according 'breaks' argument. It returns then the intervals with
#' the corresponding bin's index.
#' 
#' The range of bins is determined by 'breaks' argument. For example:
#' 'breaks=c(x1, x2, x3, x4)' represents three different intervals (bins): (x1,
#' x2], (x2, x3], (x3, x4].
#' 
#' If 'include.lowest' is 'TRUE' the the lowest value will be included in the
#' first interval, i.e. in [x1, x2].
#' 
#' If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
#' set. Use this parameter if the result size exceeds the limits of the
#' physical memory.
#' 
#' @param expr track expression
#' @param breaks breaks that determine the bin
#' @param intervals genomic scope for which the function is applied
#' @param include.lowest if 'TRUE', the lowest value of the range determined by
#' breaks is included
#' @param iterator track expression iterator. If 'NULL' iterator is determined
#' implicitly based on track expression.
#' @param band track expression band. If 'NULL' no band is used.
#' @param intervals.set.out intervals set name where the function result is
#' optionally outputed
#' @return If 'intervals.set.out' is 'NULL' a set of intervals with an
#' additional column that indicates the corresponding bin index.
#' @seealso \code{\link{gscreen}}, \code{\link{gextract}},
#' \code{\link{glookup}}, \code{\link{gdist}}
#' @keywords ~partition
#' @examples
#' 
#' gdb.init_examples()
#' breaks <- seq(0, 0.2, by = 0.05)
#' gpartition("dense_track", breaks, gintervals(1, 0, 5000))
#' 
#' @export gpartition
gpartition <- function(expr = NULL, breaks = NULL, intervals = NULL, include.lowest = FALSE, iterator = NULL, band = NULL, intervals.set.out = NULL) {
	if (is.null(substitute(expr)) || is.null(breaks) || is.null(intervals))
		stop("Usage: gpartition(expr, breaks, intervals, include.lowest = FALSE, iterator = NULL, band = NULL, intervals.set.out = NULL)", call. = F);
	.gcheckroot()

	exprstr <- do.call(.gexpr2str, list(substitute(expr)), envir = parent.frame())
	.iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())
	intervals.set.out <- do.call(.gexpr2str, list(substitute(intervals.set.out)), envir = parent.frame())
	if (!is.null(intervals.set.out))
		fullpath <- .gintervals.check_new_set(intervals.set.out)

	# intervals can be NULL if piped with gscreen and the latter returns NULL
	success <- FALSE
	res <- NULL
	tryCatch({
		if (!is.null(intervals)) {
			res <- .gcall("gpartition", intervals, exprstr, breaks, include.lowest, .iterator, band, intervals.set.out, new.env(parent = parent.frame()))
			if (!is.null(intervals.set.out) && .gintervals.is_bigset(intervals.set.out, F) && !.gintervals.needs_bigset(intervals.set.out))
				.gintervals.big2small(intervals.set.out)
		}
		success <- TRUE
	},
	finally = {
		if (!success && !is.null(intervals.set.out))
			unlink(fullpath, recursive = TRUE)
	})

	# refresh the list of GINTERVS, etc.
	if (!is.null(intervals.set.out)) {
		.gdb.add_intervals.set(intervals.set.out)
		retv <- 0 # suppress return value
	} else
		res
}



#' Calculates quantiles of a track expression
#' 
#' Calculates the quantiles of a track expression for the given percentiles.
#' 
#' This function calculates the quantiles for the given percentiles.
#' 
#' If data size exceeds the limit (see: 'getOption(gmax.data.size)'), the data
#' is randomly sampled to fit the limit. A warning message is generated. The
#' seed of the pseudo-random generator can be controled through 'grnd.seed'
#' option.
#' 
#' Note: this function is capable to run in multitasking mode. Sampling may
#' vary according to the extent of multitasking. Since multitasking depends on
#' the number of available CPU cores, running the function on two different
#' machines might give different results. Please switch off multitasking if you
#' want to achieve identical results on any machine. For more information
#' regarding multitasking please refer "User Manual".
#' 
#' @param expr track expression
#' @param percentiles an array of percentiles of quantiles in [0, 1] range
#' @param intervals genomic scope for which the function is applied
#' @param iterator track expression iterator. If 'NULL' iterator is determined
#' implicitly based on track expression.
#' @param band track expression band. If 'NULL' no band is used.
#' @return An array that represent quantiles.
#' @seealso \code{\link{gbins.quantiles}}, \code{\link{gintervals.quantiles}},
#' \code{\link{gdist}}
#' @keywords ~quantiles ~percentiles
#' @examples
#' 
#' gdb.init_examples()
#' gquantiles("dense_track", c(0.1, 0.6, 0.8), gintervals(c(1, 2)))
#' 
#' @export gquantiles
gquantiles <- function(expr = NULL, percentiles = 0.5, intervals = get("ALLGENOME"), iterator = NULL, band = NULL) {
	if (is.null(substitute(expr)))
		stop("Usage: gquantiles(expr, percentiles = 0.5, intervals = ALLGENOME, iterator = NULL, band = NULL)", call. = F)
	.gcheckroot()

	exprstr <- do.call(.gexpr2str, list(substitute(expr)), envir = parent.frame())
	.iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())
	
	if (.ggetOption("gmultitasking"))
		res <- .gcall("gquantiles_multitask", intervals, exprstr, percentiles, .iterator, band, new.env(parent = parent.frame()))
	else
		res <- .gcall("gquantiles", intervals, exprstr, percentiles, .iterator, band, new.env(parent = parent.frame()))
	res
}



#' Returns values from a lookup table based on track expression
#' 
#' Evaluates track expression and translates the values into bin indices that
#' are used in turn to retrieve and return values from a lookup table.
#' 
#' This function evaluates the track expression for all iterator intervals and
#' translates this value into an index based on the breaks. This index is then
#' used to address the lookup table and return the according value. More than
#' one 'expr'-'breaks' pair can be used. In that case 'lookup_table' is
#' addressed in a multidimensional manner, i.e. 'lookup_table[i1, i2, ...]'.
#' 
#' The range of bins is determined by 'breaks' argument. For example: 'breaks =
#' c(x1, x2, x3, x4)' represents three different intervals (bins): (x1, x2],
#' (x2, x3], (x3, x4].
#' 
#' If 'include.lowest' is 'TRUE' then the lowest value is included in the first
#' interval, i.e. in [x1, x2].
#' 
#' 'force.binning' parameter controls what should be done when the value of
#' 'expr' exceeds the range determined by 'breaks'. If 'force.binning' is
#' 'TRUE' then values smaller than the minimal break will be translated to
#' index 1, and the values exceeding the maximal break will be translated to
#' index 'M-1' where 'M' is the number of breaks. If 'force.binning' is 'FALSE'
#' the out-of-range values will produce 'NaN' values.
#' 
#' Regardless of 'force.binning' value if the value of 'expr' is 'NaN' then
#' result is 'NaN' too.
#' 
#' The order inside the result might not be the same as the order of intervals.
#' Use 'intervalID' column to refer to the index of the original interval from
#' the supplied 'intervals'.
#' 
#' If 'intervals.set.out' is not 'NULL' the result (without 'columnID' column)
#' is saved as an intervals set. Use this parameter if the result size exceeds
#' the limits of the physical memory.
#' 
#' @param lookup_table a multi-dimensional array containing the values that are
#' returned by the function
#' @param expr track expression
#' @param breaks breaks that determine the bin
#' @param intervals genomic scope for which the function is applied
#' @param include.lowest if 'TRUE', the lowest value of the range determined by
#' breaks is included
#' @param force.binning if 'TRUE', the values smaller than the minimal break
#' will be translated to index 1, and the values that exceed the maximal break
#' will be translated to index N-1 where N is the number of breaks. If 'FALSE'
#' the out-of-range values will produce NaN values.
#' @param iterator track expression iterator. If 'NULL' iterator is determined
#' implicitly based on track expressions.
#' @param band track expression band. If 'NULL' no band is used.
#' @param intervals.set.out intervals set name where the function result is
#' optionally outputed
#' @return If 'intervals.set.out' is 'NULL' a set of intervals with additional
#' 'value' and 'columnID' columns.
#' @seealso \code{\link{gtrack.lookup}}, \code{\link{gextract}},
#' \code{\link{gpartition}}, \code{\link{gdist}}
#' @keywords ~lookup ~extract
#' @examples
#' 
#' gdb.init_examples()
#' 
#' ## one-dimensional lookup table
#' breaks1 = seq(0.1, 0.2, length.out = 6)
#' glookup(1:5, "dense_track", breaks1, gintervals(1, 0, 200))
#' 
#' ## two-dimensional lookup table
#' t <- array(1:15, dim = c(5, 3))
#' breaks2 = seq(0.31, 0.37, length.out = 4)
#' glookup(t, "dense_track", breaks1, "2 * dense_track", breaks2,
#'         gintervals(1, 0, 200))
#' 
#' @export glookup
glookup <- function(lookup_table = NULL, ..., intervals = NULL, include.lowest = FALSE, force.binning = TRUE, iterator = NULL, band = NULL, intervals.set.out = NULL) {
	args <- as.list(substitute(list(...)))[-1L]
	if (is.null(lookup_table) || length(args) < 2 || (!is.null(intervals) && length(args) %% 2 != 0) || (is.null(intervals) && length(args) %% 2 == 0))
		stop("Usage: glookup(lookup_table, [expr, breaks]+, intervals, include.lowest = FALSE, force.binning = TRUE, iterator = NULL, band = NULL, intervals.set.out = NULL)", call. = F)
	.gcheckroot()

	if (length(args) %% 2 != 0)
		intervals <- eval.parent(args[[length(args)]])
	
	exprs <- c()
	breaks <- list()
	
	for (i in (0 : (length(args) / 2 - 1))) {
		exprs <- append(exprs, do.call(.gexpr2str, list(args[[i * 2 + 1]]), envir = parent.frame()))
		breaks[[length(breaks) + 1]] <- eval.parent(args[[i * 2 + 2]])
	}

	.iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())
	intervals.set.out <- do.call(.gexpr2str, list(substitute(intervals.set.out)), envir = parent.frame())
	if (!is.null(intervals.set.out))
		fullpath <- .gintervals.check_new_set(intervals.set.out)

	# intervals can be NULL if gextract is piped with gscreen and the latter returns NULL
	success <- FALSE
	res <- NULL
	tryCatch({
		if (!is.null(intervals)) {
			res <- .gcall("gbintransform", intervals, exprs, breaks, include.lowest, force.binning, lookup_table, .iterator, band, intervals.set.out, new.env(parent = parent.frame()))
			if (!is.null(intervals.set.out) && .gintervals.is_bigset(intervals.set.out, F) && !.gintervals.needs_bigset(intervals.set.out))
				.gintervals.big2small(intervals.set.out)
		}
		success <- TRUE
	},
	finally = {
		if (!success && !is.null(intervals.set.out))
			unlink(fullpath, recursive = TRUE)
	})

	# refresh the list of GINTERVS, etc.
	if (!is.null(intervals.set.out)) {
		.gdb.add_intervals.set(intervals.set.out)
		retv <- 0 # suppress return value
	} else
		res
}



#' Returns samples from the values of track expression
#' 
#' Returns a sample of the specified size from the values of track expression.
#' 
#' This function returns a sample of the specified size from the values of
#' track expression. If 'n' is less than the total number of values, the data
#' is randomally sampled. The seed of the pseudo-random generator can be
#' controled through 'grnd.seed' option.
#' 
#' If 'n' is higher than the total number of values, all values are returned
#' (yet reshuffled).
#' 
#' @param expr track expression
#' @param n a number of items to choose
#' @param intervals genomic scope for which the function is applied
#' @param iterator track expression iterator. If 'NULL' iterator is determined
#' implicitly based on track expression.
#' @param band track expression band. If 'NULL' no band is used.
#' @return An array that represent quantiles.
#' @seealso \code{\link{gextract}}
#' @keywords ~sample
#' @examples
#' 
#' gdb.init_examples()
#' gsample("sparse_track", 10)
#' 
#' @export gsample
gsample <- function(expr = NULL, n = NULL, intervals = NULL, iterator = NULL, band = NULL) {
	if (is.null(substitute(expr)))
		stop("Usage: gsample(expr, n, intervals = ALLGENOME, iterator = NULL, band = NULL)", call. = F)
	.gcheckroot()

	if (is.null(intervals))
		intervals <- get("ALLGENOME")
	
	exprstr <- do.call(.gexpr2str, list(substitute(expr)), envir = parent.frame())
	.iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())
	
	.gcall("gsample", exprstr, n, intervals, .iterator, band, new.env(parent = parent.frame()))
}

gsetroot <- function(groot = NULL, dir = NULL, rescan = FALSE) {
	if (is.null(groot))
		stop("Usage: gdb.init(groot, dir = NULL, rescan = FALSE)", call. = F);

	# replace ~ mark by full path
	oldwd <- getwd()
	
	groot <- path.expand(groot)
	setwd(groot)
	groot <- getwd()  # get absolute path
	setwd(oldwd)

	if (exists("GROOT") && exists("ALLGENOME") && !is.null(get("GROOT")) && !is.null(get("ALLGENOME")) && .ggetOption(".gautocompletion", FALSE))
		.gundefine_autocompletion_vars()

	assign("ALLGENOME", NULL, envir = .GlobalEnv)
	assign("GROOT", NULL, envir = .GlobalEnv)
	
	chromsizes <- read.csv(paste(groot, "chrom_sizes.txt", sep = "/"), sep = "\t", header = F)
	colnames(chromsizes) <- c("chrom", "size")
	intervals <- data.frame(chrom = as.factor(paste("chr", as.character(chromsizes$chrom), sep = "")),
							start = 0, end = as.numeric(chromsizes$size))

	if (nrow(intervals) == 0)
		stop("chrom_sizes.txt file does not contain any chromosomes", call. = F)
	
	for (chrom in intervals$chrom) {
		if (length(grep(sprintf("^%s$", chrom), intervals$chrom)) > 1)
			stop(sprintf("Chromosome \"%s\" appears more than once in chrom_sizes.txt", chrom))
	}
	intervals <- intervals[order(intervals$chrom), ]
	rownames(intervals) <- 1 : nrow(intervals)

	cartesian <- expand.grid(1 : nrow(intervals), 1 : nrow(intervals))
	intervals2d <- cbind(intervals[cartesian[, 2], ], intervals[cartesian[, 1], ])
	names(intervals2d) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2")
	rownames(intervals2d) <- 1 : nrow(intervals2d)

	assign("ALLGENOME", list(intervals, intervals2d), envir = .GlobalEnv)
	assign("GROOT", groot, envir = .GlobalEnv)
	assign("GWD", groot, envir = .GlobalEnv)

	success <- F
	tryCatch({
		if (is.null(dir))
			.gdir.cd(paste(groot, "tracks", sep = "/"), rescan)
		else {
			if (nchar(dir) < 1)
				stop("dir argument is an empty string")
			
			c <- substr(dir, 1, 1)
			if (c == "~" || c == "/")
				.gdir.cd(dir, rescan)
			else
				.gdir.cd(paste(groot, dir, sep = "/"), rescan)
		}
		success <- T
	},
	finally = {
		if (!success) {
			assign("ALLGENOME", NULL, envir = .GlobalEnv)
			assign("GROOT", NULL, envir = .GlobalEnv)
			assign("GWD", NULL, envir = .GlobalEnv)
		}
	})
}



#' Finds intervals that match track expression
#' 
#' Finds all intervals where track expression is 'TRUE'.
#' 
#' This function finds all intervals where track expression's value is 'TRUE'.
#' 
#' If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
#' set. Use this parameter if the result size exceeds the limits of the
#' physical memory.
#' 
#' @param expr logical track expression
#' @param intervals genomic scope for which the function is applied
#' @param iterator track expression iterator. If 'NULL' iterator is determined
#' implicitly based on track expression.
#' @param band track expression band. If 'NULL' no band is used.
#' @param intervals.set.out intervals set name where the function result is
#' optionally outputed
#' @return If 'intervals.set.out' is 'NULL' a set of intervals that match track
#' expression.
#' @seealso \code{\link{gsegment}}, \code{\link{gextract}}
#' @keywords ~screen ~interval ~intervals
#' @examples
#' 
#' gdb.init_examples()
#' gscreen("dense_track > 0.2 & sparse_track < 0.4",
#'         iterator = "dense_track")
#' 
#' @export gscreen
gscreen <- function(expr = NULL, intervals = NULL, iterator = NULL, band = NULL, intervals.set.out = NULL) {
	if (is.null(substitute(expr)))
		stop("Usage: gscreen(expr, intervals = ALLGENOME, iterator = NULL, band = NULL, intervals.set.out = NULL)", call. = F);
	.gcheckroot()

	if (is.null(intervals))
		intervals <- get("ALLGENOME")

	exprstr <- do.call(.gexpr2str, list(substitute(expr)), envir = parent.frame())
	.iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())	
	intervals.set.out <- do.call(.gexpr2str, list(substitute(intervals.set.out)), envir = parent.frame())

	if (!is.null(intervals.set.out))
		fullpath <- .gintervals.check_new_set(intervals.set.out)

	success <- FALSE
	res <- NULL
	tryCatch({
		if (.ggetOption("gmultitasking"))
			res <- .gcall("gscreen_multitask", exprstr, intervals, .iterator, band, intervals.set.out, new.env(parent = parent.frame()))
		else
			res <- .gcall("gscreen", exprstr, intervals, .iterator, band, intervals.set.out, new.env(parent = parent.frame()))

		if (!is.null(intervals.set.out) && .gintervals.is_bigset(intervals.set.out, F) && !.gintervals.needs_bigset(intervals.set.out))
			.gintervals.big2small(intervals.set.out)

		success <- TRUE
	},
	finally = {
		if (!success && !is.null(intervals.set.out))
			unlink(fullpath, recursive = TRUE)
	})

	# refresh the list of GINTERVS, etc.
	if (!is.null(intervals.set.out)) {
		.gdb.add_intervals.set(intervals.set.out)
		retv <- 0 # suppress return value
	} else
		res
}



#' Divides track expression into segments
#' 
#' Divides the values of track expression into segments by using Wilcoxon test.
#' 
#' This function divides the values of track expression into segments, where
#' each segment size is at least of 'minsegment' size and the P-value of
#' comparing the segment with the first 'minsegment' values from the next
#' segment is at most 'maxpval'. Comparison is done using Wilcoxon (also known
#' as Mann-Whitney) test.
#' 
#' If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
#' set. Use this parameter if the result size exceeds the limits of the
#' physical memory.
#' 
#' @param expr track expression
#' @param minsegment minimal segment size
#' @param maxpval maximal P-value that separates two adjacent segments
#' @param onetailed if 'TRUE', Wilcoxon test is performed one tailed, otherwise
#' two tailed
#' @param intervals genomic scope for which the function is applied
#' @param iterator track expression iterator of "fixed bin" type. If 'NULL'
#' iterator is determined implicitly based on track expression.
#' @param intervals.set.out intervals set name where the function result is
#' optionally outputed
#' @return If 'intervals.set.out' is 'NULL' a set of intervals where each
#' interval represents a segment.
#' @seealso \code{\link{gscreen}}, \code{\link{gwilcox}}
#' @keywords ~segment ~wilcoxon ~Mann-Whitney
#' @examples
#' 
#' gdb.init_examples()
#' gsegment("dense_track", 5000, 0.0001)
#' 
#' @export gsegment
gsegment <- function(expr = NULL, minsegment = NULL, maxpval = 0.05, onetailed = TRUE, intervals = NULL, iterator = NULL, intervals.set.out = NULL) {
	if (is.null(substitute(expr)) || is.null(minsegment))
		stop("Usage: gsegment(expr, minsegment, maxpval = 0.05, onetailed = TRUE, intervals = ALLGENOME, iterator = NULL, intervals.set.out = NULL)", call. = F)
	.gcheckroot()

	if (is.null(intervals))
		intervals <- get("ALLGENOME")

	exprstr <- do.call(.gexpr2str, list(substitute(expr)), envir = parent.frame())
	.iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())
	intervals.set.out <- do.call(.gexpr2str, list(substitute(intervals.set.out)), envir = parent.frame())
	if (!is.null(intervals.set.out))
		fullpath <- .gintervals.check_new_set(intervals.set.out)
	
	if (!onetailed)
		maxpval <- maxpval / 2
	
	# intervals can be NULL if piped with gscreen and the latter returns NULL
	success <- FALSE
	res <- NULL
	tryCatch({
		if (!is.null(intervals)) {
			res <- .gcall("gsegment", exprstr, intervals, minsegment, qnorm(maxpval), onetailed, .iterator, intervals.set.out, new.env(parent = parent.frame()))
			if (!is.null(intervals.set.out) && .gintervals.is_bigset(intervals.set.out, F) && !.gintervals.needs_bigset(intervals.set.out))
				.gintervals.big2small(intervals.set.out)
		}
		success <- TRUE
	},
	finally = {
		if (!success && !is.null(intervals.set.out))
			unlink(fullpath, recursive = TRUE)
	})

	# refresh the list of GINTERVS, etc.
	if (!is.null(intervals.set.out)) {
		.gdb.add_intervals.set(intervals.set.out)
		retv <- 0 # suppress return value
	} else
		res
}



#' Calculates summary statistics of track expression
#' 
#' Calculates summary statistics of track expression.
#' 
#' This function returns summary statistics of a track expression: total number
#' of bins, total number of bins whose value is NaN, min, max, sum, mean and
#' standard deviation of the values.
#' 
#' @param expr track expression
#' @param intervals genomic scope for which the function is applied
#' @param iterator track expression iterator. If 'NULL' iterator is determined
#' implicitly based on track expression.
#' @param band track expression band. If 'NULL' no band is used.
#' @return An array that represents summary statistics.
#' @seealso \code{\link{gintervals.summary}}, \code{\link{gbins.summary}}
#' @keywords ~summary ~statistics
#' @examples
#' 
#' gdb.init_examples()
#' gsummary("rects_track")
#' 
#' @export gsummary
gsummary <- function(expr = NULL, intervals = NULL, iterator = NULL, band = NULL) {
	if (is.null(substitute(expr)))
		stop("Usage: gsummary(expr, intervals = ALLGENOME, iterator = NULL, band = NULL)", call. = F)
	.gcheckroot()

	if (is.null(intervals))
		intervals <- get("ALLGENOME")
	
	exprstr <- do.call(.gexpr2str, list(substitute(expr)), envir = parent.frame())
	.iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())
	
	if (.ggetOption("gmultitasking"))
		res <- .gcall("gtracksummary_multitask", exprstr, intervals, .iterator, band, new.env(parent = parent.frame()))
	else
		res <- .gcall("gtracksummary", exprstr, intervals, .iterator, band, new.env(parent = parent.frame()))
	res
}



#' Sets input mode and auto-completion for track expressions, track names,
#' virtual tracks and interval sets
#' 
#' Sets input mode and auto-completion for track expressions, track names,
#' virtual tracks and interval sets.
#' 
#' This function enables / disables auto-completion for track names, virtual
#' tracks and interval sets. It also controls whether these objects together
#' with track expressions should be passed as strings or "as is" to the various
#' package functions.
#' 
#' If 'autocompletion' is 'TRUE' all the track names, virtual tracks and
#' intervals sets are defined as R variables (auxiliary variables) which allows
#' them to be auto-completed by TAB key. The values of these variables are
#' meaningless for the user and they should not be altered.
#' 
#' If 'interactive' is 'TRUE' track names, virtual tracks, interval sets and
#' track expressions are passed to the package functions "as is", i.e.
#' unquoted.
#' 
#' 'autocompletion' is required to be switched on if 'interactive' mode is on
#' too.
#' 
#' Please beware of the consequences of using interactive mode as it creates a
#' bunch of new variables in R environment. Though collision with the existing
#' variables is checked at the time of the call to 'gset_input_mode', yet
#' nothing prevents the user to modify the value of the auxiliary variables
#' later. This might cause unexpected behaviour in some of the package
#' functions. Also the auxiliary variables are automatically undefined once the
#' interactive mode is switched off. User who mistakenly uses auxiliary
#' variables to store the data might therefore accidentially loose it.
#' 
#' @return None.
#' @keywords ~autocompletion ~interactive ~quoted
#' @examples
#' 
#' gdb.init_examples()
#' gset_input_mode(interactive = FALSE)
#' gsummary("dense_track + 10")
#' gset_input_mode(autocompletion = TRUE, interactive = TRUE)
#' gsummary(dense_track + 10)
#' 
#' # roll-back to default input mode
#' gset_input_mode()
#' 
#' @export gset_input_mode
gset_input_mode <- function(autocompletion = FALSE, interactive = FALSE) {
    if (!is.logical(autocompletion) || !is.logical(interactive) || length(autocompletion) != 1 || length(interactive) != 1)
		stop("Usage: gset_input_mode(autocompletion = FALSE, interactive = FALSE)", call. = F)

	if (!autocompletion && interactive)
		stop("Autocompletion must be switched on in interactive mode", call. = F)

	if (.ggetOption(".gautocompletion") != autocompletion) {
		if (autocompletion) {
			tracks <- NULL
			intervals <- NULL
			if (exists("GTRACKS"))
				tracks <- get("GTRACKS")
			if (exists("GINTERVS"))
				intervals <- get("GINTERVS")
			.gdefine_autocompletion_vars(tracks, intervals, gvtrack.ls(), interactive)
		} else
			.gundefine_autocompletion_vars()
	}

	options(.gautocompletion = autocompletion)
	options(.ginteractive = interactive)
}



#' Prints call stack of the last uncaught error
#' 
#' Prints call stack of the last uncaught error in a friendly way.
#' 
#' Similarly to 'traceback' this function prints the call stack of the last
#' uncaught error. Yet 'gtraceback' does it in a more friendly way by omitting
#' the calls that occurred inside the library.
#' 
#' @param x see 'traceback'
#' @param max.lines see 'traceback'
#' @return See 'traceback'.
#' @seealso \code{\link{traceback}}
#' @keywords ~trace ~error ~exception
#' @examples
#' 
#' gdb.init_examples()
#' f <- function() { gscreen("blablabla") }
#' f()
#' gtraceback()
#' 
#' @export gtraceback
gtraceback <- function(x = NULL, max.lines = getOption("deparse.max.lines")) {
	x <- NULL

    if (is.null(x) && (exists(".Traceback", envir = baseenv()))) 
        x <- get(".Traceback", envir = baseenv())
	
	if (!is.null(x) && length(x) > 0) {
		# get the call stack and concatenate all complex commands together
		x <- sapply(x, paste, collapse="")

		# extract call stack function names
		fnames <- gsub("^(\\S+)\\s*\\(.*\\)$", "\\1", x, perl = T)

		# get the indices of lib functions
		libindices <- which(fnames %in% get(".GFUNCS"))

		# cut whatever comes after the first lib function
		if (length(libindices) > 0) {
			x <- get(".Traceback")[libindices[length(libindices)] : length(get(".Traceback"))]
		}
	}

	traceback(x, max.lines)
}



#' Downloads files from FTP server
#' 
#' Downloads multiple files from FTP server
#' 
#' This function downloads files from FTP server given by 'url'. The address in
#' 'url' can contain wildcards to download more than one file at once. Files
#' are downloaded to a directory given by 'path' argument.  If 'path' is
#' 'NULL', file are downloaded into 'GROOT/downloads'.
#' 
#' @param url URL of FTP server
#' @param path directory path where the downloaded files are stored
#' @return An array of file names that have been downloaded.
#' @seealso \code{\link{gtrack.import_set}}
#' @keywords ~ftp
#' @examples
#' 
#' gdb.init_examples()
#' gwget("ftp://hgdownload.cse.ucsc.edu/logs/2012/access_log.2012071*")
#' 
#' @export gwget
gwget <- function(url = NULL, path = NULL) {
	if (is.null(url))
		stop("Usage: gwget(url, path = NULL)", call. = F)

	if (is.null(path)) {
		.gcheckroot()
		path <- paste(get("GROOT"), "/downloads", sep = "")
		dir.create(path, showWarnings = F, recursive = T, mode = "0777")
	}

	if (!length(grep("^ftp\\:\\/\\/(\\w+(\\.\\w+)+)\\/(.+)", url, perl=T)))
		url <- paste("ftp://", url, sep="")
	
	if (!length(grep("^ftp\\:\\/\\/(\\w+(\\.\\w+)+)\\/(.+)", url, perl=T)))
		stop("Invalid format of URL", call. = F);
	
	old.files <- dir(path)
	oldwd <- getwd()
	tryCatch({
		setwd(path)
		
		server <- gsub("^ftp\\:\\/\\/(\\w+(\\.\\w+)+)\\/(.+)", "\\1", url, perl=T)
		rpath <- gsub("^ftp\\:\\/\\/(\\w+(\\.\\w+)+)\\/(.+)", "\\3", url, perl=T)

		rdir <- dirname(rpath)
		rfiles <- basename(rpath)
		if (system(paste("/bin/sh -c \"",
					 "ftp -i -n -v",
					 server,
					 "<< EOF\n",
					 "user anonymous anonymous\n",
					 "cd", rdir, "\n",
					 paste("mget \"", rfiles, "\"", sep=""), "\n",
					 "quit\n",
					 "EOF\"")))
			stop("Command failed", call. = F)
	},
			 interrupt = function(interrupt){ setwd(oldwd) },
			 finally = { setwd(oldwd) })

	new.files <- dir(path)
	setdiff(new.files, old.files)
}



#' Calculates Wilcoxon test on sliding windows over track expression
#' 
#' Calculates Wilcoxon test on sliding windows over the values of track
#' expression.
#' 
#' This function runs a Wilcoxon test (also known as a Mann-Whitney test) over
#' the values of track expression in the two sliding windows having an
#' identical center. The sizes of the windows are specified by 'winsize1' and
#' 'winsize2'. 'gwilcox' returns intervals where the smaller window tested
#' against a larger window gives a P-value below 'maxpval'. The test can be one
#' or two tailed.
#' 
#' 'what2find' argument controls what should be searched: peaks, lows or both.
#' 
#' If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
#' set. Use this parameter if the result size exceeds the limits of the
#' physical memory.
#' 
#' @param expr track expression
#' @param winsize1 number of values in the first sliding window
#' @param winsize2 number of values in the second sliding window
#' @param maxpval maximal P-value
#' @param onetailed if 'TRUE', Wilcoxon test is performed one tailed, otherwise
#' two tailed
#' @param what2find if '-1', lows are searched. If '1', peaks are searched. If
#' '0', both peaks and lows are searched
#' @param intervals genomic scope for which the function is applied
#' @param iterator track expression iterator of "fixed bin" type. If 'NULL'
#' iterator is determined implicitly based on track expression.
#' @param intervals.set.out intervals set name where the function result is
#' optionally outputed
#' @return If 'intervals.set.out' is 'NULL' a data frame representing the
#' intervals with an additional 'pval' column where P-value is below 'maxpval'.
#' @seealso \code{\link{gscreen}}, \code{\link{gsegment}}
#' @keywords ~wilcoxon ~Mann-Whitney
#' @examples
#' 
#' gdb.init_examples()
#' gwilcox("dense_track", 100000, 1000, maxpval = 0.01,
#'         what2find = 1)
#' 
#' @export gwilcox
gwilcox <- function(expr = NULL, winsize1 = NULL, winsize2 = NULL, maxpval = 0.05, onetailed = TRUE, what2find = 1, intervals = NULL, iterator = NULL, intervals.set.out = NULL) {
	if (is.null(substitute(expr)) || is.null(winsize1) || is.null(winsize2))
		stop("Usage: gwilcox(expr, winsize1, winsize2, maxpval = 0.05, onetailed = TRUE, what2find = 1 (-1=lows, 0=lows/highs, 1=highs), intervals = ALLGENOME, iterator = NULL, intervals.set.out = NULL)", call. = F)
	.gcheckroot()

	if (is.null(intervals))
		intervals <- get("ALLGENOME")

	exprstr <- do.call(.gexpr2str, list(substitute(expr)), envir = parent.frame())
	.iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())
	intervals.set.out <- do.call(.gexpr2str, list(substitute(intervals.set.out)), envir = parent.frame())
	if (!is.null(intervals.set.out))
		fullpath <- .gintervals.check_new_set(intervals.set.out)
	
	if (!onetailed)
		maxpval <- maxpval / 2
	
	# intervals can be NULL if piped with gscreen and the latter returns NULL
	success <- FALSE
	res <- NULL
	tryCatch({
		if (!is.null(intervals)) {
			res <- .gcall("gwilcox", exprstr, intervals, winsize1, winsize2, qnorm(maxpval), onetailed, as.integer(what2find), .iterator, intervals.set.out, new.env(parent = parent.frame()))
			if (!is.null(intervals.set.out) && .gintervals.is_bigset(intervals.set.out, F) && !.gintervals.needs_bigset(intervals.set.out))
				.gintervals.big2small(intervals.set.out)
		}
		success <- TRUE
	},
	finally = {
		if (!success && !is.null(intervals.set.out))
			unlink(fullpath, recursive = TRUE)
	})

	# refresh the list of GINTERVS, etc.
	if (!is.null(intervals.set.out)) {
		.gdb.add_intervals.set(intervals.set.out)
		retv <- 0 # suppress return value
	} else
		res
}



#' Calculates quantiles of a track expression for bins
#' 
#' Calculates quantiles of a track expression for bins.
#' 
#' This function is a binned version of 'gquantiles'. For each iterator
#' interval the value of 'bin_expr' is calculated and assigned to the
#' corresponding bin determined by 'breaks'. The quantiles of 'expr' are
#' calculated then separatedly for each bin.
#' 
#' The bins can be multi-dimensional depending on the number of
#' 'bin_expr'-'breaks' pairs.
#' 
#' The range of bins is determined by 'breaks' argument. For example:
#' 'breaks=c(x1, x2, x3, x4)' represents three different intervals (bins): (x1,
#' x2], (x2, x3], (x3, x4].
#' 
#' If 'include.lowest' is 'TRUE' the the lowest value will be included in the
#' first interval, i.e. in [x1, x2].
#' 
#' @param bin_expr a track expression that determines the bin
#' @param breaks breaks that define the bins
#' @param expr track expression for which quantiles are calculated
#' @param percentiles an array of percentiles of quantiles in [0, 1] range
#' @param intervals genomic scope for which the function is applied.
#' @param include.lowest if 'TRUE', the lowest value of the range determined by
#' breaks is included
#' @param iterator track expression iterator. If 'NULL' iterator is determined
#' implicitly based on track expressions.
#' @param band track expression band. If 'NULL' no band is used.
#' @return Multi-dimensional array representing quantiles for each percentile
#' and bin.
#' @seealso \code{\link{gquantiles}}, \code{\link{gintervals.quantiles}},
#' \code{\link{gdist}}
#' @keywords ~quantiles ~percentiles
#' @examples
#' 
#' gdb.init_examples()
#' gbins.quantiles("dense_track", c(0, 0.2, 0.4, 2), "sparse_track",
#'                 percentiles = c(0.2, 0.5),
#'                 intervals = gintervals(1),
#'                 iterator = "dense_track")
#' 
#' @export gbins.quantiles
gbins.quantiles <- function(..., expr = NULL, percentiles = 0.5, intervals = get("ALLGENOME"), include.lowest = FALSE, iterator = NULL, band = NULL) {
	args <- as.list(substitute(list(...)))[-1L]

	if (length(args) >= 0 && length(args) %% 2 != 0)
		expr <- args[[length(args)]]
		
	if (length(args) < 2 || is.null(substitute(expr)))
		stop("Usage: gbins.quantiles([bin_expr, breaks]+, expr, percentiles = 0.5, intervals = ALLGENOME, include.lowest = FALSE, iterator = NULL, band = NULL)", call. = F)
	.gcheckroot()

	exprs <- c()
	breaks <- list()

	exprs <- append(exprs, do.call(.gexpr2str, list(substitute(expr)), envir = parent.frame()))
	for (i in (0 : ((length(args) - 1) / 2 - 1))) {
		exprs <- append(exprs, do.call(.gexpr2str, list(args[[i * 2 + 1]]), envir = parent.frame()))
		breaks[[length(breaks) + 1]] <- eval.parent(args[[i * 2 + 2]])
	}

	.iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())

	res <- .gcall("gbins_quantiles", exprs, breaks, include.lowest, percentiles, intervals, .iterator, band, new.env(parent = parent.frame()))
	attr(res, "breaks") <- breaks
	res
}



#' Calculates summary statistics of a track expression for bins
#' 
#' Calculates summary statistics of a track expression for bins.
#' 
#' This function is a binned version of 'gsummary'. For each iterator interval
#' the value of 'bin_expr' is calculated and assigned to the corresponding bin
#' determined by 'breaks'. The summary statistics of 'expr' are calculated then
#' separatedly for each bin.
#' 
#' The bins can be multi-dimensional depending on the number of
#' 'bin_expr'-'breaks' pairs.
#' 
#' The range of bins is determined by 'breaks' argument. For example:
#' 'breaks=c(x1, x2, x3, x4)' represents three different intervals (bins): (x1,
#' x2], (x2, x3], (x3, x4].
#' 
#' If 'include.lowest' is 'TRUE' the the lowest value will be included in the
#' first interval, i.e. in [x1, x2].
#' 
#' @param bin_expr a track expression that determines the bin
#' @param breaks breaks that define the bins
#' @param expr track expression for which summary statistics is calculated
#' @param intervals genomic scope for which the function is applied
#' @param include.lowest if 'TRUE', the lowest value of the range determined by
#' breaks is included
#' @param iterator track expression iterator. If 'NULL' iterator is determined
#' implicitly based on track expressions.
#' @param band track expression band. If 'NULL' no band is used.
#' @return Multi-dimensional array representing summary statistics for each
#' bin.
#' @seealso \code{\link{gsummary}}, \code{\link{gintervals.summary}},
#' \code{\link{gdist}}
#' @keywords ~summary
#' @examples
#' 
#' gdb.init_examples()
#' gbins.summary("dense_track", c(0, 0.2, 0.4, 2), "sparse_track",
#'               intervals = gintervals(1), iterator = "dense_track")
#' 
#' @export gbins.summary
gbins.summary <- function(..., expr = NULL, intervals = get("ALLGENOME"), include.lowest = FALSE, iterator = NULL, band = NULL) {
	args <- as.list(substitute(list(...)))[-1L]

	if (length(args) >= 0 && length(args) %% 2 != 0)
		expr <- args[[length(args)]]
		
	if (length(args) < 2 || is.null(substitute(expr)))
		stop("Usage: gbins.summary([expr, breaks]+, expr, intervals = ALLGENOME, include.lowest = FALSE, iterator = NULL, band = NULL)", call. = F)
	.gcheckroot()

	exprs <- c()
	breaks <- list()

	exprs <- append(exprs, do.call(.gexpr2str, list(substitute(expr)), envir = parent.frame()))
	for (i in (0 : ((length(args) - 1) / 2 - 1))) {
		exprs <- append(exprs, do.call(.gexpr2str, list(args[[i * 2 + 1]]), envir = parent.frame()))
		breaks[[length(breaks) + 1]] <- eval.parent(args[[i * 2 + 2]])
	}

	.iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())

	res <- .gcall("gbins_summary", exprs, breaks, include.lowest, intervals, .iterator, band, new.env(parent = parent.frame()))
	attr(res, "breaks") <- breaks
	res
}



#' Runs R commands on a cluster
#' 
#' Runs R commands on a cluster that supports SGE.
#' 
#' This function runs R commands on a cluster by distributing them among
#' cluster nodes. It must run on a machine that supports Sun Grid Engine (SGE).
#' The order in which the commands are executed can not be guaranteed,
#' therefore the commands must be inter-independent.
#' 
#' Optional flags to 'qsub' command can be passed through 'opt.flags'
#' parameter. Users are strongly recommended to use only '-l' flag as other
#' flags might interfere with those that are already used (-terse, -S, -o, -e,
#' -V). For additional information please refer to the manual of 'qsub'.
#' 
#' The maximal number of simultaneously submitted jobs is controlled by
#' 'max.jobs'.
#' 
#' Set 'debug' argument to 'TRUE to allow additional report prints.
#' 
#' 'gcluster.run' launches R on the cluster nodes to execute the commands. 'R'
#' argument specifies how R executable should be invoked.
#' 
#' @param ... R commands
#' @param opt.flags optional flags for qsub command
#' @param max.jobs maximal number of simultaneously submitted jobs
#' @param debug if 'TRUE', additional reports are printed
#' @param R command that launches R
#' @return Return value ('retv') is a list, such that 'retv[[i]]' represents
#' the result of the run of command number 'i'. Each result consists of 4
#' fields that can be accessed by 'retv[[i]]$FIELDNAME':
#' 
#' \tabular{ll}{ \emph{FIELDNAME} \tab \emph{DESCRIPTION}\cr exit.status \tab
#' Exit status of the command. Possible values: 'success', 'failure' or
#' 'interrupted'.\cr retv \tab Return value of the command.\cr stdout \tab
#' Standard output of the command.\cr stderr \tab Standard error of the
#' command. }
#' @keywords ~cluster
#' @examples
#' 
#' gdb.init_examples()
#' v <- 17
#' gcluster.run(
#'   gsummary("dense_track + v"),
#'   {
#'     intervs <- gscreen("dense_track > 0.1", gintervals(1, 2))
#'     gsummary("sparse_track", intervs)
#'   },
#'   gsummary("rects_track")
#' )
#' 
#' @export gcluster.run
gcluster.run <- function(..., opt.flags = "", max.jobs = 400, debug = FALSE, R = "R") {
	commands <- as.list(substitute(list(...))[-1L])

	if (length(commands) < 1)
		stop("Usage: gculster.run(..., opt.flags = \"\" max.jobs = 400, debug = FALSE)", call. = F)

	if (!length(system("which qsub", ignore.stderr = T, intern = T)))
		stop("gcluster.run must run on a host that supports Sun Grid Engine (qsub)", call. = F)
	
	.gcheckroot()

	tmp.dirname <- ""

	submitted.jobs <- c()
	
	tryCatch({
		tmp.dirname <- tempfile(pattern = "", tmpdir = paste(get("GROOT"), "/tmp", sep = ""))
		if (!dir.create(tmp.dirname, recursive = T, mode = "0777"))
			stop(sprintf("Failed to create a directory %s", tmp.dirname), call. = F)

	    # save the environment + options
		# parent.frame() is the environment of the caller
		cat("Preparing for distribution...\n")
		save(.GLIBDIR, file  = paste(tmp.dirname, "libdir", sep = "/"))
		vars <- ls(all.names = TRUE, envir = parent.frame())
		envir <- parent.frame()
		while (!identical(envir, .GlobalEnv)) {
			envir <- parent.env(envir)
			vars <- union(vars, ls(all.names = TRUE, envir = envir))
		}
		save(list = vars, file = paste(tmp.dirname, "envir", sep = "/"), envir = parent.frame())
		.GSGECMD <- commands
		save(.GSGECMD, file  = paste(tmp.dirname, "commands", sep = "/"))
		opts <- options()
		save(opts, file = paste(tmp.dirname, "opts", sep = "/"))

		cat("Running the commands...\n")
		completed.jobs <- c()
		progress <- -1
		repeat {
		    # submit the commands
			num.running.jobs <- length(submitted.jobs) - length(completed.jobs)
			if (length(submitted.jobs) < length(commands) && num.running.jobs < max.jobs) {
				istart <- length(submitted.jobs) + 1
				iend <- min(length(commands), istart + (max.jobs - num.running.jobs) - 1)
				
				for (i in istart : iend) {
					out.file <- sprintf("%s/%d.out", tmp.dirname, i)
					err.file <- sprintf("%s/%d.err", tmp.dirname, i)
					script <- paste(get(".GLIBDIR"), "exec", "sgjob.sh", sep = "/")
					command <- sprintf("unset module; qsub -terse -S /bin/bash -o %s -e %s -V %s %s %d '%s' '%s'",
						out.file, err.file, opt.flags, script, i, tmp.dirname, R)
					jobid <- system(command, intern = TRUE)
					if (length(jobid) != 1)
						stop("Failed to run qsub", call. = FALSE)
					if (debug)
						cat(sprintf("\tSubmitted job %d (id: %s)\n", i, jobid))
					submitted.jobs <- c(submitted.jobs, jobid)
				}
			}

		    # monitor progress
			Sys.sleep(3)
			running.jobs <- .gcluster.running.jobs(submitted.jobs)
			
			old.completed.jobs <- completed.jobs
			completed.jobs <- setdiff(submitted.jobs, running.jobs)
			if (debug) {
				delta.jobs <- setdiff(completed.jobs, old.completed.jobs)
				if (length(delta.jobs) > 0) {
					for (jobid in delta.jobs)
						cat(sprintf("\tJob %d (id: %s) completed\n", match(jobid, submitted.jobs), jobid))
				}

				if (!length(running.jobs) && length(submitted.jobs) == length(commands))
					break

				new.progress <- length(completed.jobs)
				if (new.progress != progress) {
					progress <- new.progress
					cat(sprintf("\t%d job(s) still in progress\n", length(commands) - progress))
				}
			} else {
				if (!length(running.jobs) && length(submitted.jobs) == length(commands))
					break
				
				new.progress <- as.integer(100 * length(completed.jobs) / length(commands))
				if (new.progress != progress) {
					progress <- new.progress
					cat(sprintf("%d%%...", progress))
				} else
    				cat(".")
			}
		}
		if (!debug && progress != -1 && progress != 100)
			cat("100%\n")
	},
			 
	interrupt = function(interrupt) {
		cat("\n")
		stop("Command interrupted!", call. = FALSE)
	},
			 
    finally = {
		# We should perform now cleanup. If the user presses again Ctr+C "finaly" statement will be interrupted and the cleanup will
		# be incomplete. Unfortunately there is no way to block interrupts up until "finally" is done.
		# The only way is to solve the problem is to run "finally" in a loop.
		cleanup.finished <- FALSE
		while (!cleanup.finished) {
			tryCatch({
				# kill still running jobs
				if (length(submitted.jobs) > 0) {
					# pack the answer
					running.jobs <- .gcluster.running.jobs(submitted.jobs)
					answer <- c()
					for (i in 1 : length(commands)) {
						res <- list()
						res$exit.status <- NA
						res$retv <- NA
						res$stdout <- NA
						res$stderr <- NA
						
						if (submitted.jobs[i] %in% running.jobs)
							res$exit.status <- "interrupted"
						else {
							fname <- sprintf("%s/%d.retv", tmp.dirname, i)
							if (file.exists(fname)) {
								load(fname)
								res$exit.status <- "success"
								res$retv <- retv
							} else
								res$exit.status <- "failure"
						}

						out.file <- sprintf("%s/%d.out", tmp.dirname, i)
						if (file.exists(out.file)) {
							f <- file(out.file, "rc")
							res$stdout <- readChar(f, 1e6)
							close(f)
						}
					
						err.file <- sprintf("%s/%d.err", tmp.dirname, i)
						if (file.exists(err.file)) {
							f <- file(err.file, "rc")
							res$stderr <- readChar(f, 1e6)
							close(f)
						}

						answer[[i]] <- res
					}
					for (job in running.jobs)
						system(sprintf("qdel %s", job), ignore.stderr = T, intern = T)
					
					unlink(tmp.dirname, recursive = TRUE)
					return (answer)
				}
				unlink(tmp.dirname, recursive = TRUE)
				cleanup.finished <- TRUE
			},
			interrupt = function(interrupt) {}
			)
		}
	})
}



#' Changes current working directory in Genomic Database
#' 
#' Changes current working directory in Genomic Database.
#' 
#' This function changes the current working directory in Genomic Database (not
#' to be confused with shell's current working directory). The list of database
#' objects - tracks, intervals, track variables - is rescanned recursively
#' under 'dir'. Object names are updated with the respect to the new current
#' working directory. Example: a track named 'subdir.dense' will be referred as
#' 'dense' once current working directory is set to 'subdir'. All virtual
#' tracks are removed.
#' 
#' @param dir directory path
#' @return None.
#' @seealso \code{\link{gdb.init}}, \code{\link{gdir.cwd}},
#' \code{\link{gdir.create}}, \code{\link{gdir.rm}},
#' \code{\link{gset_input_mode}}
#' @keywords ~db ~data ~database ~cd ~dir ~directory ~folder
#' @examples
#' 
#' gdb.init_examples()
#' gdir.cd("subdir")
#' gtrack.ls()
#' gdir.cd("..")
#' gtrack.ls()
#' 
#' @export gdir.cd
gdir.cd <- function(dir = NULL) {
	if (is.null(dir))
		stop("Usage: gdir.cd(dir)", call. = F)

	success <- FALSE
	oldgwd <- get("GWD")

	tryCatch({
		.gdir.cd(dir, TRUE)
		success <- TRUE
	},
	finally = {
		if (!success)
			.gdir.cd(oldgwd, TRUE)
	})
}



#' Creates a new directory in Genomic Database
#' 
#' Creates a new directory in Genomic Database.
#' 
#' This function creates a new directory in Genomic Database. Creates only the
#' last element in the specified path.
#' 
#' @param dir directory path
#' @param showWarnings see 'dir.create'
#' @param mode see 'dir.create'
#' @return None.
#' @note A new directory cannot be created within an existing track directory.
#' @seealso \code{\link{dir.create}}, \code{\link{gdb.init}},
#' \code{\link{gdir.cwd}}, \code{\link{gdir.rm}}
#' @keywords ~db ~data ~database ~dir ~directory ~folder ~create
#' @export gdir.create
gdir.create <- function(dir = NULL, showWarnings = TRUE, mode = "0777") {
	if (is.null(dir))
		stop("Usage: gdir.create(dir, showWarnings = TRUE, mode = \"0777\")", call. = F)

	oldwd <- getwd()
	setwd(get("GWD"))
	tryCatch({
		d <- dirname(dir)
	
		if (!file.exists(d))
			stop(sprintf("Path %s does not exist.\nNote: recursive directory creation is forbidden.", d), call. = F)

		t <- .gfindtrackinpath(d)
		if (!is.null(t))
			stop(sprintf("Cannot create a directory within a track %s", t), call. = F)

		if (length(grep("\\.track$", basename(dir))) > 0)
			stop("gdir.create cannot create track directories", call. = F)

		dir.create(dir, showWarnings = showWarnings, recursive = FALSE, mode = mode)
	},
			 interrupt = function(interrupt){ setwd(oldwd) },
			 finally = { setwd(oldwd) })
}



#' Returns the current working directory in Genomic Database
#' 
#' Returns the absolute path of the current working directory in Genomic
#' Database.
#' 
#' This function returns the absolute path of the current working directory in
#' Genomic Database (not to be confused with shell's current working
#' directory).
#' 
#' @return A character string of the path.
#' @seealso \code{\link{gdb.init}}, \code{\link{gdir.cd}},
#' \code{\link{gdir.create}}, \code{\link{gdir.rm}}
#' @keywords ~db ~data ~database ~cwd ~pwd ~dir ~directory ~folder
#' @export gdir.cwd
gdir.cwd <- function() {
	.gcheckroot()
	get("GWD")
}



#' Deletes a directory from Genomic Database
#' 
#' Deletes a directory from Genomic Database.
#' 
#' This function deletes a directory from Genomic Database. If 'recursive' is
#' 'TRUE', the directory is deleted with all the files/directories it contains.
#' If the directory contains tracks or intervals, the user is prompted to
#' confirm the deletion. Set 'force' to 'TRUE' to suppress the prompt.
#' 
#' @param dir directory path
#' @param recursive if 'TRUE', the directory is deleted recursively
#' @param force if 'TRUE', supresses user confirmation of tracks/intervals
#' removal
#' @return None.
#' @seealso \code{\link{gdb.init}}, \code{\link{gdir.create}},
#' \code{\link{gdir.cd}}, \code{\link{gdir.cwd}}
#' @keywords ~db ~data ~database ~dir ~directory ~folder ~rm
#' @export gdir.rm
gdir.rm <- function(dir = NULL, recursive = FALSE, force = FALSE) {
	if (is.null(dir))
		stop("Usage: gdir.rm(dir, recursive = FALSE, force = FALSE)", call. = F)
	
	oldwd <- getwd()
	setwd(get("GWD"))
	tryCatch({
		if (!file.exists(dir)) {
			if (force)
				return(invisible())
			stop(sprintf("Directory %s does not exist", dir), call. = F)
		}
	
		r <- file.info(dir)
		if (r[names(r) == "isdir"] != 1)
			stop(sprintf("%s is not a directory", dir), call. = F)

		t <- .gfindtrackinpath(dir)
		if (!is.null(t))
			stop(sprintf("Directory %s belongs to track %s", dir, t), call. = F)

		answer <- "Y"

		if (recursive && !force) {
			res <- .gcall("gfind_tracks_n_intervals", dir, new.env(parent = parent.frame()), silent = TRUE)
			tracks <- res[[1]]
			intervals <- res[[2]]
	
			if (!force && length(tracks) + length(intervals) > 0) {
				cat(sprintf("Directory %s contains tracks or intervals. Are you still sure you want to delete it (Y/N)? ", dir))
				answer <- toupper(readLines(n = 1))
			}
		}
	
		if (answer == "Y" || answer == "YES") {
			if (recursive)
				unlink(dir, recursive)
			else
				file.remove(dir)
			
			if (file.exists(dir))
				stop("Failed to remove the directory", call. = F)
		}
		gdb.reload()
	},
			 interrupt = function(interrupt){ setwd(oldwd) },
			 finally = { setwd(oldwd) })
}



#' Creates a set of 1D intervals
#' 
#' Creates a set of 1D intervals.
#' 
#' This function returns a set of one-dimensional intervals. The returned value
#' can be used in all functions that accept 'intervals' argument.
#' 
#' One-dimensional intervals is a data frame whose first three columns are
#' 'chrom', 'start' and 'end'. Each row of the data frame represents a genomic
#' interval of the specified chromosome in the range of [start, end).
#' Additional columns can be presented in 1D intervals object yet these columns
#' must be added after the three obligatory ones.
#' 
#' If 'strands' argument is not 'NULL' an additional column "strand" is added
#' to the intervals. The possible values of a strand can be '1' (plus strand),
#' '-1' (minus strand) or '0' (unknown).
#' 
#' @param chroms chromosomes - an array of strings with or without "chr"
#' prefixes or an array of integers (like: '1' for "chr1")
#' @param starts an array of start coordinates
#' @param ends an array of end coordinates. If '-1' chromosome size is assumed.
#' @param strands 'NULL' or an array consisting of '-1', '0' or '1' values
#' @return A data frame representing the intervals.
#' @seealso \code{\link{gintervals.2d}}, \code{\link{gintervals.force_range}}
#' @keywords ~intervals
#' @examples
#' 
#' gdb.init_examples()
#' 
#' ## the following 3 calls produce identical results
#' gintervals(1)
#' gintervals("1")
#' gintervals("chrX")
#' 
#' gintervals(1, 1000)
#' gintervals(c("chr2", "chrX"), 10, c(3000, 5000))
#' 
#' @export gintervals
gintervals <- function(chroms = NULL, starts = 0, ends = -1, strands = NULL) {
	if (is.null(chroms))
		stop("Usage: gintervals(chroms, starts = 0, ends = -1, strands = NULL)", call. = F)
	.gcheckroot()

	intervals <- .gintervals(chroms, starts, ends, strands)
	.gcall("gintervsort", intervals, new.env(parent = parent.frame()))
}



#' Creates a set of 2D intervals
#' 
#' Creates a set of 2D intervals.
#' 
#' This function returns a set of two-dimensional intervals. The returned value
#' can be used in all functions that accept 'intervals' argument.
#' 
#' Two-dimensional intervals is a data frame whose first six columns are
#' 'chrom1', 'start1', 'end1', 'chrom2', 'start2' and 'end2'. Each row of the
#' data frame represents two genomic intervals from two chromosomes in the
#' range of [start, end). Additional columns can be presented in 2D intervals
#' object yet these columns must be added after the six obligatory ones.
#' 
#' @param chroms1 chromosomes1 - an array of strings with or without "chr"
#' prefixes or an array of integers (like: '1' for "chr1")
#' @param starts1 an array of start1 coordinates
#' @param ends1 an array of end1 coordinates. If '-1' chromosome size is
#' assumed.
#' @param chroms2 chromosomes2 - an array of strings with or without "chr"
#' prefixes or an array of integers (like: '1' for "chr1"). If 'NULL',
#' 'chroms2' is assumed to be equal to 'chroms1'.
#' @param starts2 an array of start2 coordinates
#' @param ends2 an array of end2 coordinates. If '-1' chromosome size is
#' assumed.
#' @return A data frame representing the intervals.
#' @seealso \code{\link{gintervals}}, \code{\link{gintervals.force_range}}
#' @keywords ~intervals
#' @examples
#' 
#' gdb.init_examples()
#' 
#' ## the following 3 calls produce identical results
#' gintervals.2d(1)
#' gintervals.2d("1")
#' gintervals.2d("chrX")
#' 
#' gintervals.2d(1, 1000, 2000, "chrX", 400, 800)
#' gintervals.2d(c("chr2", "chrX"), 10, c(3000, 5000), 1)
#' 
#' @export gintervals.2d
gintervals.2d <- function(chroms1 = NULL, starts1 = 0, ends1 = -1, chroms2 = NULL, starts2 = 0, ends2 = -1) {
	if (is.null(chroms1))
		stop("Usage: gintervals.2d(chroms1, starts1 = 0, ends1 = -1, chroms2 = NULL, starts2 = 0, ends2 = -1)", call. = F)
	.gcheckroot()

	if (is.null(chroms2))
		chroms2 <- chroms1

	intervals1 <- .gintervals(chroms1, starts1, ends1, NULL)
	intervals2 <- .gintervals(chroms2, starts2, ends2, NULL)

	intervals <- data.frame(chrom1 = intervals1$chrom, start1 = intervals1$start, end1 = intervals1$end,
							chrom2 = intervals2$chrom, start2 = intervals2$start, end2 = intervals2$end)

	.gcall("gintervsort", intervals, new.env(parent = parent.frame()))
}



#' Returns 2D intervals that cover the whole genome
#' 
#' Returns 2D intervals that cover the whole genome.
#' 
#' This function returns a set of two-dimensional intervals that cover the
#' whole genome as it is defined by 'chrom_sizes.txt' file.
#' 
#' @return A data frame representing the intervals.
#' @seealso \code{\link{gintervals.2d}}
#' @keywords ~genome ~chromosome ~chromosomes ~ALLGENOME
#' @export gintervals.2d.all
gintervals.2d.all <- function() {
	.gcheckroot()
	get("ALLGENOME")[[2]]
}



#' Returns 1D intervals that cover the whole genome
#' 
#' Returns 1D intervals that cover the whole genome.
#' 
#' This function returns a set of one-dimensional intervals that cover the
#' whole genome as it is defined by 'chrom_sizes.txt' file.
#' 
#' @return A data frame representing the intervals.
#' @seealso \code{\link{gintervals}}
#' @keywords ~genome ~chromosome ~chromosomes ~ALLGENOME
#' @export gintervals.all
gintervals.all <- function() {
	.gcheckroot()
	get("ALLGENOME")[[1]]
}



#' Intersects two-dimenstional intervals with a band
#' 
#' Intersects two-dimenstional intervals with a band.
#' 
#' This function intersects each two-dimensional interval from 'intervals' with
#' 'band'. If the intersection is not empty, the interval is shrunk to the
#' minimal rectangle that contains the band and added to the return value.
#' 
#' If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
#' set. Use this parameter if the result size exceeds the limits of the
#' physical memory.
#' 
#' @param intervals two-dimensional intervals
#' @param band track expression band. If 'NULL' no band is used.
#' @param intervals.set.out intervals set name where the function result is
#' optionally outputed
#' @return If 'intervals.set.out' is 'NULL' a data frame representing the
#' intervals.
#' @seealso \code{\link{gintervals.2d}}, \code{\link{gintervals.intersect}}
#' @keywords ~band ~intersect
#' @examples
#' 
#' gdb.init_examples()
#' gintervals.2d.band_intersect(gintervals.2d(1), c(10000, 20000))
#' 
#' @export gintervals.2d.band_intersect
gintervals.2d.band_intersect <- function(intervals = NULL, band = NULL, intervals.set.out = NULL) {
	if (is.null(intervals))
		stop("Usage: gintervals.2d.band_intersect(intervals, band = NULL, intervals.set.out = NULL)", call. = F)

	intervals.set.out <- do.call(.gexpr2str, list(substitute(intervals.set.out)), envir = parent.frame())

	if (!is.null(intervals.set.out))
		fullpath <- .gintervals.check_new_set(intervals.set.out)

	# intervals can be NULL if gextract is piped with gscreen and the latter returns NULL
	success <- FALSE
	res <- NULL
	tryCatch({
		if (!is.null(intervals)) {
			res <- .gcall("ginterv_intersectband", intervals, band, intervals.set.out, new.env(parent = parent.frame()))
			if (!is.null(intervals.set.out) && .gintervals.is_bigset(intervals.set.out, F) && !.gintervals.needs_bigset(intervals.set.out))
				.gintervals.big2small(intervals.set.out)
		}
		success <- TRUE
	},
	finally = {
		if (!success && !is.null(intervals.set.out))
			unlink(fullpath, recursive = TRUE)
	})

	# refresh the list of GINTERVS, etc.
	if (!is.null(intervals.set.out)) {
		.gdb.add_intervals.set(intervals.set.out)
		retv <- 0 # suppress return value
	} else
		res
}



#' Converts intervals to canonic form
#' 
#' Converts intervals to canonic form.
#' 
#' This function converts 'intervals' into a "canonic" form: properly sorted
#' with no overlaps. The result can be used later in the functions that require
#' the intervals to be in canonic form. Use 'unify_touching_intervals' to
#' control whether the intervals that touch each other (i.e. the end coordinate
#' of one equals to the start coordinate of the other) are unified.
#' 'unify_touching_intervals' is ignored if two-dimensional intervals are used.
#' 
#' Since 'gintervals.canonic' unifies overlapping or touching intervals, the
#' number of the returned intervals might be less than the number of the
#' original intervals. To allow the user to find the origin of the new interval
#' 'mapping' attribute is attached to the result. It maps between the original
#' intervals and the resulted intervals. Use 'attr(retv_of_gintervals.canonic,
#' "mapping")' to retrieve the map.
#' 
#' @param intervals intervals to be converted
#' @param unify_touching_intervals if 'TRUE', touching one-dimensional
#' intervals are unified
#' @return A data frame representing the canonic intervals and an attribute
#' 'mapping' that maps the original intervals to the resulted ones.
#' @seealso \code{\link{gintervals}}, \code{\link{gintervals.2d}}
#' @keywords ~intervals ~canonic
#' @examples
#' 
#' gdb.init_examples()
#' 
#' ## Create intervals manually by using 'data.frame'.
#' ## Note that we add an additional column 'data'.
#' ## Return value:
#' ##   chrom start   end data
#' ## 1  chr1 11000 12000   10
#' ## 2  chr1   100   200   20
#' ## 3  chr1 10000 13000   30
#' ## 4  chr1 10500 10600   40
#' intervs <- data.frame(chrom = "chr1",
#'                       start = c(11000, 100, 10000, 10500),
#'                       end = c(12000, 200, 13000, 10600),
#'                       data = c(10, 20, 30, 40))
#' 
#' ## Convert the intervals into the canonic form.
#' ## The function discards any columns besides chrom, start and end.
#' ## Return value:
#' ##  chrom start   end
#' ## 1  chr1   100   200
#' ## 2  chr1 10000 13000
#' res <- gintervals.canonic(intervs)
#' 
#' ## By inspecting mapping attribute we can see how the new
#' ## intervals were created: "2 1 2 2" means that the first
#' ## interval in the result was created from the second interval in
#' ## the original set (we look for the indices in mapping where "1"
#' ## appears). Likewise the second interval in the result was
#' ## created from 3 intervals in the original set. Their indices are
#' ## 1, 3 and 4 (once again we look for the indices in mapping where
#' ## "2" appears).
#' ## Return value:
#' ## 2 1 2 2
#' attr(res, "mapping")
#' 
#' ## Finally (and that is the most useful part of 'mapping'
#' ## attribute): we add a new column 'data' to our result which is
#' ## the mean value of the original data column. The trick is done
#' ## using 'tapply' on par with 'mapping' attribute. For example,
#' ## 20.00000 equals is a result of 'mean(intervs[2,]$data' while
#' ## 26.66667 is a result of 'mean(intervs[c(1,3,4),]$data)'.
#' ## 'res' after the following call:
#' ##  chrom start   end     data
#' ## 1  chr1   100   200 20.00000
#' ## 2  chr1 10000 13000 26.66667
#' res$data <- tapply(intervs$data, attr(res, "mapping"), mean)
#' 
#' @export gintervals.canonic
gintervals.canonic <- function(intervals = NULL, unify_touching_intervals = TRUE) {
	if (is.null(intervals))
		stop("Usage: gintervals.canonic(intervals, unify_touching_intervals = TRUE)", call. = F)

	res <- .gcall("gintervcanonic", intervals, unify_touching_intervals, new.env(parent = parent.frame()))
	res
}



#' Calculates difference of two intervals sets
#' 
#' Returns difference of two sets of intervals.
#' 
#' This function returns a genomic space that is covered by 'intervals1' but
#' not covered by 'intervals2'.
#' 
#' If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
#' set. Use this parameter if the result size exceeds the limits of the
#' physical memory.
#' 
#' @param intervals1,intervals2 set of one-dimensional intervals
#' @param intervals.set.out intervals set name where the function result is
#' optionally outputed
#' @return If 'intervals.set.out' is 'NULL' a data frame representing the
#' intervals.
#' @seealso \code{\link{gintervals}}, \code{\link{gintervals.intersect}},
#' \code{\link{gintervals.union}}
#' @keywords ~diff
#' @examples
#' 
#' gdb.init_examples()
#' 
#' intervs1 <- gscreen("dense_track > 0.15")
#' intervs2 <- gscreen("dense_track < 0.2")
#' 
#' ## 'res3' equals to 'res4'
#' res3 <- gintervals.diff(intervs1, intervs2)
#' res4 <- gscreen("dense_track >= 0.2")
#' 
#' @export gintervals.diff
gintervals.diff <- function(intervals1 = NULL, intervals2 = NULL, intervals.set.out = NULL) {
	if (is.null(intervals1) || is.null(intervals2))
		stop("Usage: gintervals.diff(intervals1, intervals2, intervals.set.out = NULL)", call. = F)

	intervals.set.out <- do.call(.gexpr2str, list(substitute(intervals.set.out)), envir = parent.frame())

	if (.gintervals.is_bigset(intervals1) || .gintervals.is_bigset(intervals2) || !is.null(intervals.set.out)) {
		res <- NULL

		FUN <- function(intervals, intervals.set.out, envir) {
			intervals1 <- intervals[[1]]
			intervals2 <- intervals[[2]]
			chrom_res <- .gcall("gintervdiff", intervals1, intervals2, new.env(parent = parent.frame()))
			if (!is.null(chrom_res) && nrow(chrom_res) > 0) {
				if (is.null(intervals.set.out)) {
					assign("res", c(get("res", envir = envir), list(chrom_res)), envir = envir)
					.gverify_max_data_size(sum(unlist(lapply(get("res", envir), nrow))), arguments = "intervals.set.out")
				}
			}
			chrom_res
		}

		.gintervals.apply(gintervals.chrom_sizes(intervals1), list(intervals1, intervals2), intervals.set.out, FUN, intervals.set.out, environment())

		if (!is.null(res))
			res <- do.call(.grbind, res)   # much faster than calling rbind incrementally in FUN

		if (is.null(intervals.set.out)) {
			if (!is.null(res) && nrow(res))
				res
			else
				NULL
		} else
			retv <- 0 # suppress return value
	} else {
		res <- .gcall("gintervdiff", intervals1, intervals2, new.env(parent = parent.frame()))
		res
	}
}



#' Tests for a named intervals set existence
#' 
#' Tests for a named intervals set existence.
#' 
#' This function returns 'TRUE' if a named intervals set exists in Genomic
#' Database.
#' 
#' @param intervals.set name of an intervals set
#' @return 'TRUE' if a named intervals set exists. Otherwise 'FALSE'.
#' @seealso \code{\link{gintervals.ls}}, \code{\link{gintervals.load}},
#' \code{\link{gintervals.rm}}, \code{\link{gintervals.save}},
#' \code{\link{gintervals}}, \code{\link{gintervals.2d}}
#' @keywords ~intervals
#' @examples
#' 
#' gdb.init_examples()
#' gintervals.exists("annotations")
#' 
#' @export gintervals.exists
gintervals.exists <- function(intervals.set = NULL) {
	if (is.null(substitute(intervals.set)))
		stop("Usage: gintervals.exists(intervals.set)", call. = F)
	.gcheckroot()

	intervals.set <- do.call(.gexpr2str, list(substitute(intervals.set)), envir = parent.frame())
	!is.na(match(intervals.set, get("GINTERVS")))
}



#' Limits intervals to chromosomal range
#' 
#' Limits intervals to chromosomal range.
#' 
#' This function enforces the intervals to be within the chromosomal range [0,
#' chrom length) by altering the intervals' boundaries. Intervals that lay
#' entirely outside of the chromosomal range are eliminated. The new intervals
#' are returned.
#' 
#' If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
#' set. Use this parameter if the result size exceeds the limits of the
#' physical memory.
#' 
#' @param intervals intervals
#' @param intervals.set.out intervals set name where the function result is
#' optionally outputed
#' @return If 'intervals.set.out' is 'NULL' a data frame representing the
#' intervals.
#' @seealso \code{\link{gintervals}}, \code{\link{gintervals.2d}},
#' \code{\link{gintervals.canonic}}
#' @keywords ~intervals
#' @examples
#' 
#' gdb.init_examples()
#' intervs <- data.frame(chrom = "chr1",
#'                       start = c(11000, -100, 10000, 10500),
#'                       end = c(12000, 200, 13000000, 10600))
#' gintervals.force_range(intervs)
#' 
#' @export gintervals.force_range
gintervals.force_range <- function(intervals = NULL, intervals.set.out = NULL) {
	if (is.null(substitute(intervals)))
		stop("Usage: gintervals.force_range(intervals, intervals.set.out = NULL)", call. = F)
	.gcheckroot()

	intervals.set.out <- do.call(.gexpr2str, list(substitute(intervals.set.out)), envir = parent.frame())

	res <- NULL
	FUN <- function(intervals, intervals.set.out, envir) {
		intervals <- intervals[[1]]
		if (.gintervals.is1d(intervals)) {
			intervals$start <- pmax(intervals$start, 0)
			intervals$end <- pmin(intervals$end, gintervals.all()[match(intervals$chrom, gintervals.all()$chrom), ]$end)
			intervals <- intervals[!is.na(intervals$end) & intervals$start < intervals$end, ]
		} else {
			intervals$start1 <- pmax(intervals$start1, 0)
			intervals$end1 <- pmin(intervals$end1, gintervals.all()[match(intervals$chrom1, gintervals.all()$chrom), ]$end)
			intervals$start2 <- pmax(intervals$start2, 0)
			intervals$end2 <- pmin(intervals$end2, gintervals.all()[match(intervals$chrom2, gintervals.all()$chrom), ]$end)
			intervals <- intervals[!is.na(intervals$end1) & !is.na(intervals$end2) & intervals$start1 < intervals$end1 & intervals$start2 < intervals$end2, ]
		}
		if (is.null(intervals.set.out)) {
			assign("res", c(get("res", envir = envir), list(intervals)), envir = envir)
			.gverify_max_data_size(sum(unlist(lapply(get("res", envir), nrow))), arguments = "intervals.set.out")
		}
		intervals
	}

	.gintervals.apply(gintervals.chrom_sizes(intervals), intervals, intervals.set.out, FUN, intervals.set.out, environment())

	if (!is.null(res))
		res <- do.call(.grbind, res)   # much faster than calling rbind incrementally in FUN

	if (is.null(intervals.set.out)) {
		if (!is.null(res) && nrow(res))
			res
		else
			NULL
	} else
		retv <- 0 # suppress return value
}



#' Imports genes and annotations from files
#' 
#' Imports genes and annotations from files.
#' 
#' This function reads a definition of genes from 'genes.file' and returns four
#' sets of intervals: TSS, exons, 3utr and 5utr. In addition to the regular
#' intervals columns 'strand' column is added. It contains '1' values for '+'
#' strands and '-1' values for '-' strands.
#' 
#' If annotation file 'annots.file' is given then annotations are attached too
#' to the intervals. The names of the annotations as they would appear in the
#' return value must be specified in 'annots.names' argument.
#' 
#' Both 'genes.file' and 'annots.file' can be either a file path or URL in a
#' form of 'ftp://[address]/[file]'. Files that these arguments point to can be
#' zipped or unzipped.
#' 
#' Examples of 'genes.file' and 'annots.file' can be found here:
#' 
#' \code{ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz}
#' \code{ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/kgXref.txt.gz}
#' 
#' If a few intervals overlap (for example: two TSS regions) they are all
#' unified to an interval that covers the whole overlapping region. 'strand'
#' value is set to '0' if two or more of the overlapping intervals have
#' different strands. The annotations of the overlapping intervals are
#' concatenated to a single character string separated by semicolons. Identical
#' values of overlapping intervals' annotation are eliminated.
#' 
#' @param genes.file name or URL of file that contains genes
#' @param annots.file name of URL file that contains annotations. If 'NULL' no
#' annotations are imported
#' @param annots.names annotations names
#' @return A list of four intervals sets named 'tss', 'exons', 'utr3' and
#' 'utr5'. 'strand' column and annotations are attached to the intevals.
#' @seealso \code{\link{gintervals}}, \code{\link{gdb.create}}
#' @keywords ~intervals ~import ~genes
#' @export gintervals.import_genes
gintervals.import_genes <- function(genes.file = NULL, annots.file = NULL, annots.names = NULL) {
	if (is.null(genes.file))
		stop("Usage: gintervals.import_genes(genes.file, annots.file = NULL, annots.names = NULL)", call. = F)

	tmp.dirname <- tempfile(pattern = "", tmpdir = paste(get("GROOT"), "/downloads", sep = ""))
	if (!dir.create(tmp.dirname, recursive = T, mode = "0777"))
		stop(sprintf("Failed to create a directory %s", tmp.dirname), call. = F)

	files <- list(genes.file, annots.file)
	file.types <- c("genes.file", "annots.file")

	tryCatch({
		for (i in 1 : length(files)) {
			if (is.null(files[[i]]))
				next

        	protocol <- "ftp://"
        	if (substr(files[[i]], 1, nchar(protocol)) == protocol) {
        	    # ftp
        		f <- gwget(files[[i]], tmp.dirname)
				if (length(f) != 1)
					stop(sprintf("More than one file matches %s argument", file.types[i]), call. = F)
        		files[[i]] <- paste(tmp.dirname, "/", f, sep = "")
        	}

    		if (length(grep("^.+\\.gz$", files[[i]], perl = T))) {
    			f.unzipped <- basename(gsub("^(.+)\\.gz$", "\\1", files[[i]], perl = T))
    			f.unzipped <- paste(tmp.dirname, "/", f.unzipped, sep = "")
    			cmd <- paste("/bin/sh -c \"gunzip -q -c", files[[i]], ">", f.unzipped, "\"")
    			if (system(cmd))
    				stop(sprintf("Command failed: %s", cmd), call. = F)
    			files[[i]] <- f.unzipped
    		}
		}

		res <- .gcall("gintervals_import_genes", files[[1]], files[[2]], annots.names, new.env(parent = parent.frame()))
		res
	}, finally = {
		unlink(tmp.dirname, recursive = TRUE) }
	)
}



#' Calculates an intersection of two sets of intervals
#' 
#' Calculates an intersection of two sets of intervals.
#' 
#' This function returns intervals that represent a genomic space which is
#' achieved by intersection of 'intervals1' and 'intervals2'.
#' 
#' If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
#' set. Use this parameter if the result size exceeds the limits of the
#' physical memory.
#' 
#' @param intervals1,intervals2 set of intervals
#' @param intervals.set.out intervals set name where the function result is
#' optionally outputed
#' @return If 'intervals.set.out' is 'NULL' a data frame representing the
#' intersection of intervals.
#' @seealso \code{\link{gintervals.2d.band_intersect}},
#' \code{\link{gintervals.diff}}, \code{\link{gintervals.union}},
#' \code{\link{gintervals}}, \code{\link{gintervals.2d}}
#' @keywords ~intersect
#' @examples
#' 
#' gdb.init_examples()
#' 
#' intervs1 <- gscreen("dense_track > 0.15")
#' intervs2 <- gscreen("dense_track < 0.2")
#' 
#' ## 'intervs3' and 'intervs4' are identical
#' intervs3 <- gintervals.intersect(intervs1, intervs2)
#' intervs4 <- gscreen("dense_track > 0.15 & dense_track < 0.2")
#' 
#' @export gintervals.intersect
gintervals.intersect <- function(intervals1 = NULL, intervals2 = NULL, intervals.set.out = NULL) {
	if (is.null(intervals1) || is.null(intervals2))
		stop("Usage: gintervals.intersect(intervals1, intervals2, intervals.set.out = NULL)", call. = F)

	intervals.set.out <- do.call(.gexpr2str, list(substitute(intervals.set.out)), envir = parent.frame())

	if (.gintervals.is_bigset(intervals1) || .gintervals.is_bigset(intervals2) || !is.null(intervals.set.out)) {
		res <- NULL

		FUN <- function(intervals, intervals.set.out, envir) {
			intervals1 <- intervals[[1]]
			intervals2 <- intervals[[2]]
			chrom_res <- .gcall("gintervintersect", intervals1, intervals2, new.env(parent = parent.frame()))
			if (!is.null(chrom_res) && nrow(chrom_res) > 0) {
				if (is.null(intervals.set.out)) {
					assign("res", c(get("res", envir = envir), list(chrom_res)), envir = envir)
					.gverify_max_data_size(sum(unlist(lapply(get("res", envir), nrow))), arguments = "intervals.set.out")
				}
			}
			chrom_res
		}

		chroms1 <- gintervals.chrom_sizes(intervals1)
		chroms1$size <- NULL
		chroms2 <- gintervals.chrom_sizes(intervals2)
		chroms2$size <- NULL
		.gintervals.apply(merge(chroms1, chroms2), list(intervals1, intervals2), intervals.set.out, FUN, intervals.set.out, environment())

		if (!is.null(res))
			res <- do.call(.grbind, res)   # much faster than calling rbind incrementally in FUN

		if (is.null(intervals.set.out)) {
			if (!is.null(res) && nrow(res))
				res
			else
				NULL
		} else
			retv <- 0 # suppress return value
	} else {
		res <- .gcall("gintervintersect", intervals1, intervals2, new.env(parent = parent.frame()))
		res
	}
}



#' Returns number of intervals per chromosome
#' 
#' Returns number of intervals per chromosome (or chromosome pair).
#' 
#' This function returns number of intervals per chromosome (for 1D intervals)
#' or chromosome pair (for 2D intervals).
#' 
#' @param intervals intervals set
#' @return Data frame representing number of intervals per chromosome (for 1D
#' intervals) or chromosome pair (for 2D intervals).
#' @seealso \code{\link{gintervals.load}}, \code{\link{gintervals.save}},
#' \code{\link{gintervals.exists}}, \code{\link{gintervals.ls}},
#' \code{\link{gintervals}}, \code{\link{gintervals.2d}}
#' @keywords ~intervals
#' @examples
#' 
#' gdb.init_examples()
#' gintervals.chrom_sizes("annotations")
#' 
#' @export gintervals.chrom_sizes
gintervals.chrom_sizes <- function(intervals = NULL) {
	if (is.null(intervals))
		stop("Usage: gintervals.chrom_sizes(intervals)", call. = F)
	.gcheckroot()

	if (.gintervals.is_bigset(intervals)) {
		if (.gintervals.big.is1d(intervals)) {
			stats <- .gintervals.big.meta(intervals)$stats
			res <- stats[, match(c("chrom", "size"), colnames(stats))]
		} else {
			stats <- .gintervals.big.meta(intervals)$stats
			res <- stats[, match(c("chrom1", "chrom2", "size"), colnames(stats))]
		}
	} else
		res <- .gcall("gintervals_chrom_sizes", .gintervals.load_ext(intervals), new.env(parent = parent.frame()))

	if (nrow(res) > 1)
		rownames(res) <- 1 : nrow(res)
	res
}



#' Tests for big intervals set
#' 
#' Tests for big intervals set.
#' 
#' This function tests whether 'intervals.set' is a big intervals set.
#' Intervals set is big if it is stored in big intervals set format and given
#' the current limits it cannot be fully loaded into memory.
#' 
#' Memory limit is controlled by 'gmax.data.size' option (see:
#' 'getOption("gmax.data.size")').
#' 
#' @param intervals.set name of an intervals set
#' @return 'TRUE' if intervals set is big, otherwise 'FALSE'.
#' @seealso \code{\link{gintervals.load}}, \code{\link{gintervals.save}},
#' \code{\link{gintervals.exists}}, \code{\link{gintervals.ls}}
#' @keywords ~intervals
#' @examples
#' 
#' gdb.init_examples()
#' gintervals.is.bigset("annotations")
#' 
#' @export gintervals.is.bigset
gintervals.is.bigset <- function(intervals.set = NULL) {
	if (is.null(intervals.set))
		stop("Usage: gintervals.is.bigset(intervals.set)", call. = F)
	.gcheckroot()

	.gintervals.is_bigset(intervals.set) && !.gintervals.loadable(intervals.set)
}



#' Converts intervals from another assembly
#' 
#' Converts intervals from another assembly to the current one.
#' 
#' This function converts 'intervals' from another assembly to the current one.
#' Chain file instructs how the conversion of coordinates should be done. It
#' can be either a name of a chain file or a data frame in the same format as
#' returned by 'gintervals.load_chain' function.
#' 
#' The converted intervals are returned. An additional column named
#' 'intervalID' is added to the resulted data frame. For each interval in the
#' resulted intervals it indicates the index of the original interval.
#' 
#' @param intervals intervals from another assembly
#' @param chain name of chain file or data frame as returned by
#' 'gintervals.load_chain'
#' @return A data frame representing the converted intervals.
#' @seealso \code{\link{gintervals.load_chain}}, \code{\link{gtrack.liftover}},
#' \code{\link{gintervals}}
#' @keywords ~intervals ~liftover ~chain
#' @examples
#' 
#' gdb.init_examples()
#' chainfile <- paste(GROOT, "data/test.chain", sep = "/")
#' intervs <- data.frame(chrom = "chr25", start = c(0, 7000),
#'                       end = c(6000, 20000))
#' gintervals.liftover(intervs, chainfile)
#' 
#' @export gintervals.liftover
gintervals.liftover <- function(intervals = NULL, chain = NULL) {
	if (is.null(intervals) || is.null(chain))
		stop("Usage: gintervals.liftover(intervals, chain)", call. = F)
	.gcheckroot()

	if (is.character(chain))
		chain.intervs <- gintervals.load_chain(chain)
	else
		chain.intervs <- f
	
	.gcall("gintervs_liftover", intervals, chain.intervs, new.env(parent = parent.frame()))
}



#' Loads a named intervals set
#' 
#' Loads a named intervals set.
#' 
#' This function loads and returns intervals stored in a named intervals set.
#' 
#' If intervals set contains 1D intervals and 'chrom' is not 'NULL' only the
#' intervals of the given chromosome are returned.
#' 
#' Likewise if intervals set contains 2D intervals and 'chrom1', 'chrom2' are
#' not 'NULL' only the intervals of the given pair of chromosomes are returned.
#' 
#' For big intervals sets 'chrom' parameter (1D case) / 'chrom1', 'chrom2'
#' parameters (2D case) must be specified. In other words: big intervals sets
#' can be loaded only by chromosome or chromosome pair.
#' 
#' @param intervals.set name of an intervals set
#' @param chrom chromosome for 1D intervals set
#' @param chrom1 first chromosome for 2D intervals set
#' @param chrom2 second chromosome for 2D intervals set
#' @return A data frame representing the intervals.
#' @seealso \code{\link{gintervals.save}}, \code{\link{gintervals.is.bigset}},
#' \code{\link{gintervals.exists}}, \code{\link{gintervals.ls}},
#' \code{\link{gintervals}}, \code{\link{gintervals.2d}}
#' @keywords ~intervals
#' @examples
#' 
#' gdb.init_examples()
#' gintervals.load("annotations")
#' 
#' @export gintervals.load
gintervals.load <- function(intervals.set = NULL, chrom = NULL, chrom1 = NULL, chrom2 = NULL) {
	.gintervals.load_ext(intervals.set, chrom, chrom1, chrom2, TRUE)
}



#' Loads assembly conversion table from a chain file
#' 
#' Loads assembly conversion table from a chain file.
#' 
#' This function reads a file in 'chain' format and returns assembly conversion
#' table that can be used in 'gtrack.liftover' and 'gintervals.liftover'.
#' 
#' Note: chain file might map a few different source intervals into a single
#' target one. These ambigous mappings are not presented in the data frame
#' returned by 'gintervals.load_chain'.
#' 
#' @param file name of chain file
#' @return A data frame representing assembly conversion table.
#' @seealso \code{\link{gintervals.liftover}}, \code{\link{gtrack.liftover}}
#' @keywords ~intervals ~liftover ~chain
#' @examples
#' 
#' gdb.init_examples()
#' chainfile <- paste(GROOT, "data/test.chain", sep = "/")
#' gintervals.load_chain(chainfile)
#' 
#' @export gintervals.load_chain
gintervals.load_chain <- function(file = NULL) {
	if (is.null(file))
		stop("Usage: gintervals.load_chain(file)", call. = F)
	.gcall("gchain2interv", file, new.env(parent = parent.frame()))
}



#' Returns a list of named intervals sets
#' 
#' Returns a list of named intervals sets in Genomic Database.
#' 
#' This function returns a list of named intervals sets that match the pattern
#' (see 'grep'). If called without any arguments all named intervals sets are
#' returned.
#' 
#' @param pattern,ignore.case,perl,fixed,useBytes see 'grep'
#' @return An array that contains the names of intervals sets.
#' @seealso \code{\link{grep}}, \code{\link{gintervals.exists}},
#' \code{\link{gintervals.load}}, \code{\link{gintervals.save}},
#' \code{\link{gintervals.rm}}, \code{\link{gintervals}},
#' \code{\link{gintervals.2d}}
#' @keywords ~intervals ~ls
#' @examples
#' 
#' gdb.init_examples()
#' gintervals.ls()
#' gintervals.ls(pattern = "annot*")
#' 
#' @export gintervals.ls
gintervals.ls <- function(pattern = "", ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE) {
	.gcheckroot()
	grep(pattern, get("GINTERVS"), value = TRUE, ignore.case = ignore.case, perl = perl, fixed = fixed, useBytes = useBytes)
}



#' Applies a function to values of track expressions
#' 
#' Applies a function to values of track expressions for each interval.
#' 
#' This function evaluates track expressions for each interval from
#' 'intervals'. The resulted vectors are passed then as arguments to 'FUN'.
#' 
#' If the intervals are one-dimensional and have an additional column named
#' 'strand' whose value is '-1', the values of the track expression are placed
#' to the vector in reverse order.
#' 
#' The current interval index (1-based) is stored in 'GAPPLY.INTERVID' variable
#' that is available during the execution of 'gintervals.mapply'. There is no
#' guarantee about the order in which the intervals are processed. Do not rely
#' on any specific order and use 'GITERATOR.INTERVID' variable to detect the
#' current interval id.
#' 
#' If 'enable.gapply.intervals' is 'TRUE', an additional variable
#' 'GAPPLY.INTERVALS' is defined during the execution of 'gintervals.mapply'.
#' This variable stores the current iterator intervals prior to track
#' expression evaluation. Please note that setting 'enable.gapply.intervals' to
#' 'TRUE' might severely affect the run-time of the function.
#' 
#' Note: all the changes made in R environment by 'FUN' will be void if
#' multitasking mode is switched on. One should also refrain from performing
#' any other operations in 'FUN' that might be not "thread-safe" such as
#' updating files, etc. Please switch off multitasking ('options(gmultitasking
#' = FALSE)') if you wish to perform such operations.
#' 
#' If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
#' set. Use this parameter if the result size exceeds the limits of the
#' physical memory.
#' 
#' @param FUN function to apply, found via ‘match.fun’
#' @param expr track expressions whose values are used as arguments for 'FUN'
#' @param intervals intervals for which track expressions are calculated
#' @param enable.gapply.intervals if 'TRUE', then a variable 'GAPPLY.INTERVALS'
#' is available
#' @param iterator track expression iterator. If 'NULL' iterator is determined
#' implicitly based on track expressions.
#' @param intervals.set.out intervals set name where the function result is
#' optionally outputed
#' @return If 'intervals.set.out' is 'NULL' a data frame representing intervals
#' with an additional column that contains the return values of 'FUN'.
#' @seealso \code{\link{mapply}}
#' @keywords ~apply ~mapply
#' @examples
#' 
#' gdb.init_examples()
#' gintervals.mapply(max, "dense_track",
#'                   gintervals(c(1, 2), 0, 10000))
#' gintervals.mapply(function(x, y){ max(x + y) }, "dense_track",
#'                   "sparse_track", gintervals(c(1, 2), 0, 10000),
#'                   iterator = "sparse_track")
#' 
#' @export gintervals.mapply
gintervals.mapply <- function(FUN = NULL, ..., intervals = NULL, enable.gapply.intervals = F, iterator = NULL, band = NULL, intervals.set.out = NULL) {
	assign("GINTERVID", -1, envir = .GlobalEnv)
	args <- as.list(substitute(list(...)))[-1L]
	if (is.null(intervals) && length(args) < 2 || !is.null(intervals) && length(args) < 1)
		stop("Usage: gintervals.mapply(FUN, [expr]+, intervals, enable.gapply.intervals = FALSE, iterator = NULL, intervals.set.out = NULL)", call. = F)
	.gcheckroot()

	if (is.null(intervals)) {
		intervals <- eval.parent(args[[length(args)]])
		args <- args[1 : (length(args) - 1)]
	}

	tracks <- c()
	for (track in args)
		tracks <- c(tracks, do.call(.gexpr2str, list(track), envir = parent.frame()))

	.iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())

	if (exists("GAPPLY.INTERVALS"))
		remove(list = "GAPPLY.INTERVALS", envir = .GlobalEnv)

	intervals.set.out <- do.call(.gexpr2str, list(substitute(intervals.set.out)), envir = parent.frame())

	if (.gintervals.is_bigset(intervals) || !is.null(intervals.set.out)) {
		res <- NULL

		INTERVALS_FUN <- function(intervals, intervals.set.out, envir) {
			intervals <- intervals[[1]]
			chrom_res <- .gcall("gmapply", intervals, FUN, tracks, enable.gapply.intervals, .iterator, band, FALSE, new.env(parent = parent.frame()))
			if (!is.null(chrom_res) && nrow(chrom_res) > 0) {
				if (is.null(intervals.set.out)) {
					assign("res", c(get("res", envir = envir), list(chrom_res)), envir = envir)
					.gverify_max_data_size(sum(unlist(lapply(get("res", envir), nrow))), arguments = "intervals.set.out")
				}
			}
			chrom_res
		}

		.gintervals.apply(gintervals.chrom_sizes(intervals), intervals, intervals.set.out, INTERVALS_FUN, intervals.set.out, environment())

		if (!is.null(res))
			res <- do.call(.grbind, res)   # much faster than calling rbind incrementally in FUN

		if (is.null(intervals.set.out)) {
			if (!is.null(res) && nrow(res))
				res
			else
				NULL
		} else
			retv <- 0 # suppress return value
	} else {
		if (.ggetOption("gmultitasking"))
			.gcall("gmapply_multitask", intervals, FUN, tracks, enable.gapply.intervals, .iterator, band, TRUE, new.env(parent = parent.frame()))
		else
			.gcall("gmapply", intervals, FUN, tracks, enable.gapply.intervals, .iterator, band, TRUE, new.env(parent = parent.frame()))
	}
}



#' Finds neighbors between two sets of intervals
#' 
#' Finds neighbors between two sets of intervals.
#' 
#' This function finds for each interval in 'intervals1' the closest
#' 'maxneighbors' intervals from 'intervals2'.
#' 
#' For 1D intervals the distance must fall in the range of ['mindist',
#' 'maxdist']. If 'intervals2' contains a 'strand' column the distance can be
#' positive or negative depending on the 'strand' value and the position of
#' interval2 relatively to interval1. If 'strand' column is missing the
#' distance is always positive.
#' 
#' For 2D intervals two distances are calculated and returned for each axis.
#' The distances must fall in the range of ['mindist1', 'maxdist1'] for axis 1
#' and ['mindist2', 'maxdist2'] for axis 2. For selecting the closest
#' 'maxneighbors' intervals Manhattan distance is used (i.e. dist1+dist2).
#' 
#' If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
#' set. Use this parameter if the result size exceeds the limits of the
#' physical memory.
#' 
#' @param intervals1,intervals2 intervals
#' @param maxneighbors maximal number of neighbors
#' @param mindist,maxdist distance range for 1D intervals
#' @param mindist1,maxdist1,mindist2,maxdist2 distance range for 2D intervals
#' @param na.if.notfound if 'TRUE' return 'NA' interval if no matching
#' neighbors were found, otherwise omit the interval in the answer
#' @param intervals.set.out intervals set name where the function result is
#' optionally outputed
#' @return If 'intervals.set.out' is 'NULL' a data frame containing the pairs
#' of intervals from 'intervals1', intervals from 'intervals2' and an
#' additional column named 'dist' ('dist1' and 'dist2' for 2D intervals)
#' representing the distance between the corresponding intervals. If
#' 'na.if.notfound' is 'TRUE', the data frame contains all the intervals from
#' 'intervals1' including those for which no matching neighbor was found. For
#' the latter intervals an 'NA' neighboring interval is stated. If
#' 'na.if.notfound' is 'FALSE', the data frame contains only intervals from
#' 'intervals1' for which matching neighbor(s) was found.
#' @seealso \code{\link{gintervals}},
#' @keywords ~intervals ~annotate ~nearest ~neighbor ~neighbors
#' @examples
#' 
#' gdb.init_examples()
#' intervs1 <- giterator.intervals("dense_track",
#'                                 gintervals(1, 0, 4000),
#'                                 iterator = 233)
#' intervs2 <- giterator.intervals("sparse_track",
#'                                 gintervals(1, 0, 2000))
#' gintervals.neighbors(intervs1, intervs2, 10, mindist = -300,
#'                      maxdist = 500)
#' intervs2$strand = c(1, 1, -1, 1)
#' gintervals.neighbors(intervs1, intervs2, 10, mindist = -300,
#'                      maxdist = 500)
#' 
#' @export gintervals.neighbors
gintervals.neighbors <- function(intervals1 = NULL, intervals2 = NULL, maxneighbors = 1, mindist = -1e+09, maxdist = 1e+09,
	                             mindist1 = -1e+09, maxdist1 = 1e+09, mindist2 = -1e+09, maxdist2 = 1e+09,
								 na.if.notfound = FALSE, intervals.set.out = NULL)
{
	if (is.null(intervals1) || is.null(intervals2))
		stop(paste("Usage: gintervals.neighbors(intervals1, intervals2, maxneighbors = 1, mindist = -1e+09, maxdist = 1e+09, ",
			"mindist1 = -1e+09, maxdist1 = 1e+09, mindist2 = -1e+09, maxdist2 = 1e+09, na.if.notfound = FALSE, intervals.set.out = NULL)",
			sep = ""), call. = F)

	if (is.null(colnames)) {
		intervals1name <- deparse(substitute(intervals1), width.cutoff = 500)[1]
		intervals2name <- deparse(substitute(intervals2), width.cutoff = 500)[1]
	}

	intervals.set.out <- do.call(.gexpr2str, list(substitute(intervals.set.out)), envir = parent.frame())

	if (.gintervals.is_bigset(intervals1) || .gintervals.is_bigset(intervals2) || !is.null(intervals.set.out)) {
		res <- NULL

		FUN <- function(intervals, intervals.set.out, envir) {
			intervals1 <- intervals[[1]]
			intervals2 <- intervals[[2]]
			chrom_res <- .gcall("gfind_neighbors", intervals1, intervals2, maxneighbors, mindist, maxdist, mindist1, maxdist1, mindist2, maxdist2, na.if.notfound, FALSE, new.env(parent = parent.frame()))
			if (!is.null(chrom_res) && is.null(intervals.set.out)) {
				assign("res", c(get("res", envir = envir), list(chrom_res)), envir = envir)
				.gverify_max_data_size(sum(unlist(lapply(get("res", envir), nrow))), arguments = "intervals.set.out")
			}
			chrom_res
		}

		if (na.if.notfound)
			.gintervals.apply(gintervals.chrom_sizes(intervals1), list(intervals1, intervals2), intervals.set.out, FUN, intervals.set.out, environment())
		else {
			chroms1 <- gintervals.chrom_sizes(intervals1)
			chroms1$size <- NULL
			chroms2 <- gintervals.chrom_sizes(intervals2)
			chroms2$size <- NULL
			.gintervals.apply(merge(chroms1, chroms2), list(intervals1, intervals2), intervals.set.out, FUN, intervals.set.out, environment())
		}

		if (!is.null(res))
			res <- do.call(.grbind, res)   # much faster than calling rbind incrementally in FUN

		if (is.null(intervals.set.out)) {
			if (!is.null(res) && nrow(res)) 
				res
			else
				NULL
		} else
			retv <- 0 # suppress return value
	} else {
		intervals1 <- .gintervals.load_ext(intervals1)
		intervals2 <- .gintervals.load_ext(intervals2)
		res <- .gcall("gfind_neighbors", intervals1, intervals2, maxneighbors, mindist, maxdist, mindist1, maxdist1, mindist2, maxdist2, na.if.notfound, TRUE, new.env(parent = parent.frame()))
		res
	}
}



#' Calculates quantiles of a track expression for intervals
#' 
#' Calculates quantiles of a track expression for intervals.
#' 
#' This function calculates quantiles of 'expr' for each interval in
#' 'intervals'.
#' 
#' If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
#' set. Use this parameter if the result size exceeds the limits of the
#' physical memory.
#' 
#' @param expr track expression for which quantiles are calculated
#' @param percentiles an array of percentiles of quantiles in [0, 1] range
#' @param intervals set of intervals
#' @param iterator track expression iterator. If 'NULL' iterator is determined
#' implicitly based on track expressions.
#' @param band track expression band. If 'NULL' no band is used.
#' @param intervals.set.out intervals set name where the function result is
#' optionally outputed
#' @return If 'intervals.set.out' is 'NULL' a set of intervals with additional
#' columns representing quantiles for each percentile.
#' @seealso \code{\link{gquantiles}}, \code{\link{gbins.quantiles}}
#' @keywords ~quantiles ~percentiles
#' @examples
#' 
#' gdb.init_examples()
#' intervs <- gintervals(c(1, 2), 0, 5000)
#' gintervals.quantiles("dense_track",
#'                      percentiles = c(0.5, 0.3, 0.9), intervs)
#' 
#' @export gintervals.quantiles
gintervals.quantiles <- function(expr = NULL, percentiles = 0.5, intervals = NULL, iterator = NULL, band = NULL, intervals.set.out = NULL) {
	if (is.null(substitute(expr)) || is.null(intervals))
		stop("Usage: gintervals.quantiles(expr, percentiles = 0.5, intervals, iterator = NULL, band = NULL, intervals.set.out = NULL)", call. = F)
	.gcheckroot()

	exprstr <- do.call(.gexpr2str, list(substitute(expr)), envir = parent.frame())
	.iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())
	intervals.set.out <- do.call(.gexpr2str, list(substitute(intervals.set.out)), envir = parent.frame())

	if (!is.null(intervals.set.out))
		fullpath <- .gintervals.check_new_set(intervals.set.out)

	success <- FALSE
	res <- NULL
	tryCatch({
		if (.ggetOption("gmultitasking"))
			res <- .gcall("gintervals_quantiles_multitask", intervals, exprstr, percentiles, .iterator, band, intervals.set.out, new.env(parent = parent.frame()))
		else
			res <- .gcall("gintervals_quantiles", intervals, exprstr, percentiles, .iterator, band, intervals.set.out, new.env(parent = parent.frame()))

		if (!is.null(intervals.set.out) && .gintervals.is_bigset(intervals.set.out, F) && !.gintervals.needs_bigset(intervals.set.out))
			.gintervals.big2small(intervals.set.out)

		success <- TRUE
	},
	finally = {
		if (!success && !is.null(intervals.set.out))
			unlink(fullpath, recursive = TRUE)
	})

	# refresh the list of GINTERVS, etc.
	if (is.null(intervals.set.out))
		res
	else {
		.gdb.add_intervals.set(intervals.set.out)
		retv <- 0 # suppress return value
	}
}



#' Combines several sets of intervals
#' 
#' Combines several sets of intervals into one set.
#' 
#' This function combines several intervals sets into one set. It works in a
#' similar manner as 'rbind' yet it is faster. Also it supports intervals sets
#' that are stored in files including the big intervals sets.
#' 
#' If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
#' set. If the format of the output intervals is set to be "big" (determined
#' implicitly based on the result size and options), the order of the resulted
#' intervals is altered as they are sorted by chromosome (or chromosomes pair -
#' for 2D).
#' 
#' @param intervals intervals set
#' @param intervals.set.out intervals set name where the function result is
#' optionally outputed
#' @return If 'intervals.set.out' is 'NULL' a data frame combining intervals
#' sets.
#' @seealso \code{\link{gintervals}}, \code{\link{gintervals.2d}},
#' \code{\link{gintervals.canonic}}
#' @keywords ~rbind
#' @examples
#' 
#' gdb.init_examples()
#' 
#' intervs1 <- gextract("sparse_track", gintervals(c(1, 2), 1000, 4000))
#' intervs2 <- gextract("sparse_track", gintervals(c(2, "X"), 2000, 5000))
#' gintervals.save("testintervs", intervs2)
#' gintervals.rbind(intervs1, "testintervs")
#' gintervals.rm("testintervs", force = TRUE)
#' 
#' @export gintervals.rbind
gintervals.rbind <- function(..., intervals.set.out = NULL) {
	intervals <- list(...)
	if (!length(intervals))
		stop("Usage: gintervals.rbind([intervals]+, intervals.set.out = NULL)", call. = F)
	.gcheckroot()

	intervals.set.out <- do.call(.gexpr2str, list(substitute(intervals.set.out)), envir = parent.frame())

	res <- NULL
	if (any(unlist(lapply(intervals, function(intervals) { .gintervals.is_bigset(intervals) }))) || !is.null(intervals.set.out)) {
		if (is.null(intervals.set.out)) {
			FUN <- function(intervals, intervals.set.out, envir) {
				assign("res", c(get("res", envir = envir), intervals), envir = envir)
				.gverify_max_data_size(sum(unlist(lapply(get("res", envir), nrow))), arguments = "intervals.set.out")
				intervals[[1]]
			}

			# preserve the order of intervals inside the answer
			lapply(intervals, f <- function(intervals) { .gintervals.apply(gintervals.chrom_sizes(intervals), intervals, NULL, FUN, NULL, parent.frame(2)) })
			if (!is.null(res))
				res <- do.call(.grbind, res)   # much faster than calling rbind incrementally in FUN
		} else {
			FUN <- function(intervals, intervals.set.out, envir) {
				intervals <- do.call(.grbind, intervals)
				intervals
			}

			# use for .gintervals.apply chromosomes from all intervals
			chroms <- NULL
			chroms <- lapply(intervals, gintervals.chrom_sizes)
			chroms <- do.call(rbind, chroms)
			if (.gintervals.is1d(intervals[[1]])) {
				chroms <- factor(chroms$chrom, levels(chroms$chrom))
				chroms <- unique(chroms)
				chroms <- sort(chroms)
				chroms <- data.frame(chrom = chroms)
			} else {
				chroms <- data.frame(chrom1 = chroms$chrom1, chrom2 = chroms$chrom2)
				chroms <- unique(chroms)
				chroms <- chroms[with(chroms, order(chrom1, chrom2)), ]
			}

			.gintervals.apply(chroms, intervals, intervals.set.out, FUN, intervals.set.out, environment())
		}
	} else {
		intervals <- lapply(intervals, .gintervals.load_ext)
		res <- do.call(.grbind, intervals)   # much faster than calling rbind incrementally in FUN
	}

	if (is.null(intervals.set.out)) {
		if (!is.null(res) && nrow(res))
			res
		else
			NULL
	} else
		retv <- 0 # suppress return value
}



#' Deletes a named intervals set
#' 
#' Deletes a named intervals set.
#' 
#' This function deletes a named intervals set from the Genomic Database. By
#' default 'gintervals.rm' requires the user to interactively confirm the
#' deletion. Set 'force' to 'TRUE' to suppress the user prompt.
#' 
#' @param intervals.set name of an intervals set
#' @param force if 'TRUE', supresses user confirmation of a named intervals set
#' removal
#' @return None.
#' @seealso \code{\link{gintervals.save}}, \code{\link{gintervals.exists}},
#' \code{\link{gintervals.ls}}, \code{\link{gintervals}},
#' \code{\link{gintervals.2d}}
#' @keywords ~intervals
#' @examples
#' 
#' gdb.init_examples()
#' intervs <- gintervals(c(1, 2))
#' gintervals.save("testintervs", intervs)
#' gintervals.ls()
#' gintervals.rm("testintervs", force = TRUE)
#' gintervals.ls()
#' 
#' @export gintervals.rm
gintervals.rm <- function(intervals.set = NULL, force = FALSE) {
	if (is.null(substitute(intervals.set)))
		stop("Usage: gintervals.rm(intervals.set, force = FALSE)", call. = F)
	.gcheckroot()

	intervals.set <- do.call(.gexpr2str, list(substitute(intervals.set)), envir = parent.frame())
	
	# check whether intervals.set appears among GINTERVS
	if (is.na(match(intervals.set, get("GINTERVS")))) {
		if (force)
			return(invisible())
		stop(sprintf("Intervals set %s does not exist", intervals.set), call. = F)
	}

	answer <- "N"
	if (force)
		answer <- "Y"
	else {
		str <- sprintf("Are you sure you want to delete intervals set %s (Y/N)? ", intervals.set)
		cat(str)
		answer <- toupper(readLines(n = 1))
	}
	
	if (answer == "Y" || answer == "YES") {
		fname <- sprintf("%s.interv", paste(get("GWD"), gsub("\\.", "/", intervals.set), sep = "/"))

		# remove the intervals set
		unlink(fname, recursive = TRUE)
		
		if (file.exists(fname))
			cat(sprintf("Failed to delete intervals set %s\n", intervals.set))
		else
			# refresh the list of GINTERVS, etc.
			.gdb.rm_intervals.set(intervals.set)
	}
}



#' Creates a named intervals set
#' 
#' Saves intervals to a named intervals set.
#' 
#' This function saves 'intervals' as a named intervals set.
#' 
#' @param intervals.set name of an intervals set
#' @param intervals intervals
#' @return None.
#' @seealso \code{\link{gintervals.rm}}, \code{\link{gintervals.load}},
#' \code{\link{gintervals.exists}}, \code{\link{gintervals.ls}},
#' \code{\link{gintervals}}, \code{\link{gintervals.2d}}
#' @keywords ~intervals
#' @examples
#' 
#' gdb.init_examples()
#' intervs <- gintervals(c(1, 2))
#' gintervals.save("testintervs", intervs)
#' gintervals.ls()
#' gintervals.rm("testintervs", force = TRUE)
#' 
#' @export gintervals.save
gintervals.save <- function(intervals.set.out = NULL, intervals = NULL) {
	if (is.null(substitute(intervals.set.out)) || is.null(intervals))
		stop("Usage: gintervals.save(intervals.set.out, intervals)", call. = F)
	.gcheckroot()
	
	intervals.set.out <- do.call(.gexpr2str, list(substitute(intervals.set.out)), envir = parent.frame())
	.gintervals.apply(gintervals.chrom_sizes(intervals), intervals, intervals.set.out, function(intervs, ...) { intervs[[1]] })
	retv <- NULL
}



#' Updates a named intervals set
#' 
#' Updates a named intervals set.
#' 
#' This function replaces all intervals of given chromosome (or chromosome
#' pair) within 'intervals.set' with 'intervals'. Chromosome is specified by
#' 'chrom' for 1D intervals set or 'chrom1', 'chrom2' for 2D intervals set.
#' 
#' If 'intervals' is 'NULL' all intervals of given chromosome are removed from
#' 'intervals.set'.
#' 
#' @param intervals.set name of an intervals set
#' @param intervals intervals or 'NULL'
#' @param chrom chromosome for 1D intervals set
#' @param chrom1 first chromosome for 2D intervals set
#' @param chrom2 second chromosome for 2D intervals set
#' @return None.
#' @seealso \code{\link{gintervals.save}}, \code{\link{gintervals.load}},
#' \code{\link{gintervals.exists}}, \code{\link{gintervals.ls}}
#' @keywords ~intervals
#' @examples
#' 
#' gdb.init_examples()
#' intervs <- gscreen("sparse_track > 0.2",
#'                    gintervals(c(1,2), 0, 10000))
#' gintervals.save("testintervs", intervs)
#' gintervals.load("testintervs")
#' gintervals.update("testintervs", intervs[intervs$chrom == "chr2",][1:5,], chrom = 2)
#' gintervals.load("testintervs")
#' gintervals.update("testintervs", NULL, chrom = 2)
#' gintervals.load("testintervs")
#' gintervals.rm("testintervs", force = TRUE)
#' 
#' @export gintervals.update
gintervals.update <- function(intervals.set = NULL, intervals = "", chrom = NULL, chrom1 = NULL, chrom2 = NULL) {
	if (is.null(substitute(intervals.set)) || identical(intervals, ""))
		stop("Usage: gintervals.update(intervals.set, intervals, chrom = NULL, chrom1 = NULL, chrom2 = NULL)", call. = F)
	.gcheckroot()

	if (identical(intervals.set, intervals))
		return(retv <- NULL)

	if (is.null(chrom) && is.null(chrom1) && is.null(chrom2))
		stop("Chromosome must be specified in chrom (for 2D intervals: chrom1, chrom2) parameter", call. = F)

	if (!is.null(chrom)) {
		chrom <- .gchroms(chrom)
		if (length(chrom) > 1)
			stop("chrom parameter should mark only one chromosome")
	}

	if (!is.null(chrom1)) {
		chrom1 <- .gchroms(chrom1)
		if (length(chrom1) > 1)
			stop("chrom1 parameter should mark only one chromosome")
	}

	if (!is.null(chrom2)) {
		chrom2 <- .gchroms(chrom2)
		if (length(chrom2) > 1)
			stop("chrom2 parameter should mark only one chromosome")
	}

	if (!is.null(chrom) && !is.null(chrom1))
		stop("Cannot use chrom and chrom1 parameters in the same call", call. = F)

	if (!is.null(chrom) && !is.null(chrom2))
		stop("Cannot use chrom and chrom2 parameters in the same call", call. = F)

	if (!is.character(intervals.set) || length(intervals.set) != 1)
		stop("Invalid format of intervals.set parameter", call. = F)

	 if (is.na(match(intervals.set, get("GINTERVS"))))
		stop(sprintf("Intervals set %s does not exist", intervals.set), call. = F)

	path <- gsub(".", "/", intervals.set, fixed = T)
	path <- paste(path, ".interv", sep = "")
	fullpath <- paste(get("GWD"), path, sep = "/")

	if (!is.null(intervals)) {
		if (!is.null(chrom))
			intervals <- .gintervals.load_ext(intervals, chrom = chrom)
		else
			intervals <- .gintervals.load_ext(intervals, chrom1 = chrom1, chrom2 = chrom2)
	}

	# big: update stats (including delete), save chrom (or delete), convert to small if needed
	if (.gintervals.is_bigset(intervals.set)) {
		is1d <- .gintervals.big.is1d(intervals.set)
		meta <- .gintervals.big.meta(intervals.set)
		stats <- meta$stats
		zeroline <- meta$zeroline

		if (!is.null(intervals) && !identical(sapply(intervals, "class"), sapply(zeroline, "class")))
			stop(sprintf("Cannot update intervals set %s: columns differ", intervals.set), call. = F)

		if (is1d) {
			if (is.null(chrom))
				stop("chrom parameter is not specified", call. = F)
			idx <- which(stats$chrom == chrom)
			if (length(idx) > 0)
				stats <- stats[-idx,]
			if (!is.null(intervals) && nrow(intervals)) {
				stat <- .gcall("gintervals_stats", intervals, new.env(parent = parent.frame()))
				stats <- rbind(stats, data.frame(chrom = chrom, stat))
				stats <- stats[order(stats$chrom), ]
			}
			.gintervals.big.save(fullpath, intervals, chrom = chrom)
		} else {
			if (is.null(chrom1) || is.null(chrom2))
				stop("chrom1 and chrom2 parameters must be specified", call. = F)
			idx <- which(stats$chrom1 == chrom1 & stats$chrom2 == chrom2)
			if (length(idx) > 0)
				stats <- stats[-idx,]
			if (!is.null(intervals) && nrow(intervals)) {
				stat <- .gcall("gintervals_stats", intervals, new.env(parent = parent.frame()))
				stats <- rbind(stats, data.frame(chrom1 = chrom1, chrom2 = chrom2, stat))
				stats <- stats[order(stats$chrom1, stats$chrom2), ]
			}
			.gintervals.big.save(fullpath, intervals, chrom1 = chrom1, chrom2 = chrom2)
		}

		if (nrow(stats) > 1)
			rownames(stats) <- 1 : nrow(stats)
		.gintervals.big.save_meta(fullpath, stats, zeroline)

		if (!.gintervals.needs_bigset(intervals.set))
			.gintervals.big2small(intervals.set)
	}

	# small: load all, update in place (including delete), save back, convert to big if needed
	else {
		tgt.intervals <- .gintervals.load_ext(intervals.set)
		is1d <- .gintervals.is1d(intervals.set)

		if (!is.null(intervals) && !identical(sapply(intervals, "class"), sapply(tgt.intervals, "class")))
			stop(sprintf("Cannot update intervals set %s: columns differ", intervals.set), call. = F)

		if (is1d) {
			if (is.null(chrom))
				stop("chrom parameter is not specified", call. = F)
			idx <- which(tgt.intervals$chrom == chrom)
			if (length(idx) > 0)
				tgt.intervals <- tgt.intervals[-idx,]
			if (!is.null(intervals) && nrow(intervals)) {
				tgt.intervals <- .grbind(tgt.intervals, intervals)
				tgt.intervals <- tgt.intervals[order(tgt.intervals$chrom), ]
			}
		} else {
			if (is.null(chrom1) || is.null(chrom2))
				stop("chrom1 and chrom2 parameters must be specified", call. = F)
			idx <- which(tgt.intervals$chrom1 == chrom1 & tgt.intervals$chrom2 == chrom2)
			if (length(idx) > 0)
				tgt.intervals <- tgt.intervals[-idx,]
			if (!is.null(intervals) && nrow(intervals)) {
				tgt.intervals <- .grbind(tgt.intervals, intervals)
				tgt.intervals <- tgt.intervals[order(tgt.intervals$chrom1, tgt.intervals$chrom2), ]
			}
		}
		if (.gintervals.needs_bigset(tgt.intervals))
			.gintervals.small2big(intervals.set, tgt.intervals)
		else
			.gintervals.save_file(fullpath, tgt.intervals)
	}

	retv <- 0 # suppress return value
}



#' Calculates summary statistics of track expression for intervals
#' 
#' Calculates summary statistics of track expression for intervals.
#' 
#' This function returns summary statistics of a track expression for each
#' interval 'intervals': total number of bins, total number of bins whose value
#' is NaN, min, max, sum, mean and standard deviation of the values.
#' 
#' If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
#' set. Use this parameter if the result size exceeds the limits of the
#' physical memory.
#' 
#' @param expr track expression
#' @param intervals set of intervals
#' @param iterator track expression iterator. If 'NULL' iterator is determined
#' implicitly based on track expression.
#' @param band track expression band. If 'NULL' no band is used.
#' @param intervals.set.out intervals set name where the function result is
#' optionally outputed
#' @return If 'intervals.set.out' is 'NULL' a set of intervals with additional
#' columns representing summary statistics for each percentile and interval.
#' @seealso \code{\link{gsummary}}, \code{\link{gbins.summary}}
#' @keywords ~summary ~statistics
#' @examples
#' 
#' gdb.init_examples()
#' intervs <- gintervals(c(1, 2), 0, 5000)
#' gintervals.summary("dense_track", intervs)
#' 
#' @export gintervals.summary
gintervals.summary <- function(expr = NULL, intervals = NULL, iterator = NULL, band = NULL, intervals.set.out = NULL) {
	if (is.null(substitute(expr)) || is.null(intervals))
		stop("Usage: gintervals.summary(expr, intervals, iterator = NULL, band = NULL, intervals.set.out = NULL)", call. = F)
	.gcheckroot()

	exprstr <- do.call(.gexpr2str, list(substitute(expr)), envir = parent.frame())
	.iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())
	intervals.set.out <- do.call(.gexpr2str, list(substitute(intervals.set.out)), envir = parent.frame())

	if (!is.null(intervals.set.out))
		fullpath <- .gintervals.check_new_set(intervals.set.out)
	
	# intervals can be NULL if gextract is piped with gscreen and the latter returns NULL
	success <- FALSE
	res <- NULL
	tryCatch({
		if (!is.null(intervals)) {
			res <- .gcall("gintervals_summary", exprstr, intervals, .iterator, band, intervals.set.out, new.env(parent = parent.frame()))
			if (!is.null(intervals.set.out) && .gintervals.is_bigset(intervals.set.out, F) && !.gintervals.needs_bigset(intervals.set.out))
				.gintervals.big2small(intervals.set.out)
		}
		success <- TRUE
	},
	finally = {
		if (!success && !is.null(intervals.set.out))
			unlink(fullpath, recursive = TRUE)
	})

	# refresh the list of GINTERVS, etc.
	if (!is.null(intervals.set.out)) {
		.gdb.add_intervals.set(intervals.set.out)
		retv <- 0 # suppress return value
	} else
		res
}



#' Calculates a union of two sets of intervals
#' 
#' Calculates a union of two sets of intervals.
#' 
#' This function returns intervals that represent a genomic space covered by
#' either 'intervals1' or 'intervals2'.
#' 
#' If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
#' set. Use this parameter if the result size exceeds the limits of the
#' physical memory.
#' 
#' @param intervals1,intervals2 set of one-dimensional intervals
#' @param intervals.set.out intervals set name where the function result is
#' optionally outputed
#' @return If 'intervals.set.out' is 'NULL' a data frame representing the union
#' of intervals.
#' @seealso \code{\link{gintervals.intersect}}, \code{\link{gintervals.diff}},
#' \code{\link{gintervals}}, \code{\link{gintervals.2d}}
#' @keywords ~union
#' @examples
#' 
#' gdb.init_examples()
#' 
#' intervs1 <- gscreen("dense_track > 0.15 & dense_track < 0.18")
#' intervs2 <- gscreen("dense_track >= 0.18 & dense_track < 0.2")
#' 
#' ## 'intervs3' and 'intervs4' are identical
#' intervs3 <- gintervals.union(intervs1, intervs2)
#' intervs4 <- gscreen("dense_track > 0.15 & dense_track < 0.2")
#' 
#' @export gintervals.union
gintervals.union <- function(intervals1 = NULL, intervals2 = NULL, intervals.set.out = NULL) {
	if (is.null(intervals1) || is.null(intervals2))
		stop("Usage: gintervals.union(intervals1, intervals2, intervals.set.out = NULL)", call. = F)

	if (.gintervals.is_bigset(intervals1) || .gintervals.is_bigset(intervals2) || !is.null(intervals.set.out)) {
		res <- NULL

		FUN <- function(intervals, intervals.set.out, envir) {
			intervals1 <- intervals[[1]]
			intervals2 <- intervals[[2]]
			chrom_res <- .gcall("gintervunion", intervals1, intervals2, new.env(parent = parent.frame()))
			if (!is.null(chrom_res) && nrow(chrom_res) > 0) {
				if (is.null(intervals.set.out)) {
					assign("res", c(get("res", envir = envir), list(chrom_res)), envir = envir)
					.gverify_max_data_size(sum(unlist(lapply(get("res", envir), nrow))), arguments = "intervals.set.out")
				}
			}
			chrom_res
		}

		# use for .gintervals.apply chromosomes from both intervals1 and intervals2
		chroms <- NULL
		if (.gintervals.is1d(intervals1)) {
			chroms <- rbind(gintervals.chrom_sizes(intervals1), gintervals.chrom_sizes(intervals2))
			chroms <- factor(chroms$chrom, levels(chroms$chrom))
			chroms <- unique(chroms)
			chroms <- sort(chroms)
			chroms <- data.frame(chrom = chroms)
		} else {
			chroms <- rbind(gintervals.chrom_sizes(intervals1), gintervals.chrom_sizes(intervals2))
			chroms <- data.frame(chrom1 = chroms$chrom1, chrom2 = chroms$chrom2)
			chroms <- unique(chroms)
			chroms <- chroms[with(chroms, order(chrom1, chrom2)), ]
		}
		.gintervals.apply(chroms, list(intervals1, intervals2), intervals.set.out, FUN, intervals.set.out, environment())

		if (!is.null(res))
			res <- do.call(.grbind, res)   # much faster than calling rbind incrementally in FUN

		if (is.null(intervals.set.out)) {
			if (!is.null(res) && nrow(res))
				res
			else
				NULL
		} else
			retv <- 0 # suppress return value
	} else {
		res <- .gcall("gintervunion", intervals1, intervals2, new.env(parent = parent.frame()))
		res
	}
}



#' Creates a cartesian-grid iterator
#' 
#' Creates a cartesian grid two-dimensional iterator that can be used by any
#' function that accepts an iterator argument.
#' 
#' This function creates and returns a cartesian grid two-dimensional iterator
#' that can be used by any function that accepts an iterator argument.
#' 
#' Assume 'centers1' and 'centers2' to be the central points of each interval
#' from 'intervals1' and 'intervals2', and 'C1', 'C2' to be two points from
#' 'centers1', 'centers2' accordingly. Assume also that the values in
#' 'expansion1' and 'expansion2' are unique and sorted.
#' 
#' 'giterator.cartesian_grid' creates a set of all possible unique and
#' non-overlapping two-dimensional intervals of form: '(chrom1, start1, end1,
#' chrom2, start2, end2)'. Each '(chrom1, start1, end1)' is created by taking a
#' point 'C1' - '(chrom1, coord1)' and converting it to 'start1' and 'end1'
#' such that 'start1 == coord1+E1[i]', 'end1 == coord1+E1[i+1]', where 'E1[i]'
#' is one of the sorted 'expansion1' values. Overlaps between rectangles or
#' expansion beyond the limits of chromosome are avoided.
#' 
#' 'min.band.idx' and 'max.band.idx' parameters control whether a pair of 'C1'
#' and 'C2' is skipped or not. If both of these parameters are not 'NULL' AND
#' if both 'C1' and 'C2' share the same chromosome AND the delta of indices of
#' 'C1' and 'C2' ('C1 index - C2 index') lays within '[min.band.idx,
#' max.band.idx]' range - only then the pair will be used to create the
#' intervals. Otherwise 'C1-C2' pair is filtered out. Note: if 'min.band.idx'
#' and 'max.band.idx' are not 'NULL', i.e. band indices filtering is applied,
#' then 'intervals2' parameter must be set to 'NULL'.
#' 
#' @param intervals1 one-dimensional intervals
#' @param expansion1 an array of integers that define expansion around
#' intervals1 centers
#' @param intervals2 one-dimensional intervals. If 'NULL' then 'intervals2' is
#' considered to be equal to 'intervals1'
#' @param expansion2 an array of integers that define expansion around
#' intervals2 centers. If 'NULL' then 'expansion2' is considered to be equal to
#' 'expansion1'
#' @param min.band.idx,max.band.idx integers that limit iterator intervals to
#' band
#' @return A list containing the definition of cartesian iterator.
#' @seealso \code{\link{giterator.intervals}}
#' @keywords ~iterator ~cartesian
#' @examples
#' 
#' gdb.init_examples()
#' 
#' intervs1 <- gintervals(c(1, 1, 2), c(100, 300, 200),
#'                        c(300, 500, 300))
#' intervs2 <- gintervals(c(1, 2, 2), c(400, 1000, 3000),
#'                        c(800, 2000, 4000))
#' itr <- giterator.cartesian_grid(intervs1, c(-20, 100), intervs2,
#'                                 c(-40, -10, 50))
#' giterator.intervals(iterator = itr)
#' 
#' itr <- giterator.cartesian_grid(intervs1, c(-20, 50, 100))
#' giterator.intervals(iterator = itr)
#' 
#' itr <- giterator.cartesian_grid(intervs1, c(-20, 50, 100),
#'                                 min.band.idx = -1,
#'                                 max.band.idx = 0)
#' giterator.intervals(iterator = itr)
#' 
#' @export giterator.cartesian_grid
giterator.cartesian_grid <- function(intervals1 = NULL, expansion1 = NULL, intervals2 = NULL, expansion2 = NULL, min.band.idx = NULL, max.band.idx = NULL) {
	if (is.null(intervals1) || is.null(expansion1))
		stop("Usage: giterator.cartesian_grid(intervals1, expansion1, intervals2 = NULL, expansion2 = NULL, min.band.idx = NULL, max.band.idx = NULL)", call. = F)

	use.band.idx.limit <- !is.null(min.band.idx) && !is.null(max.band.idx)
	if (use.band.idx.limit) {
		if (min.band.idx > max.band.idx)
			stop("min.band.idx exceeds max.band.idx", call. = F)

		if (!is.null(intervals2))
			stop("band.idx limit can only be used when intervals2 is set to NULL", call. = F)
	} else {
		min.band.idx <- 0
		max.band.idx <- 0
	}
	
	r <- list(intervals1 = intervals1, intervals2 = intervals2, expansion1 = expansion1, expansion2 = expansion2,
			  band.idx = c(min.band.idx, max.band.idx, use.band.idx.limit))
	class(r) <- "cartesian.grid"
	.gcall("gcheck_iterator", r, new.env(parent = parent.frame()))
	r
}



#' Returns iterator intervals
#' 
#' Returns iterator intervals given track expression, scope, iterator and band.
#' 
#' This function returns a set of intervals used by the iterator intervals for
#' the given track expression, genomic scope, iterator and band. Some functions
#' accept an iterator without accepting a track expression (like
#' 'gtrack.create_pwm_energy'). These functions generate the values for each
#' iterator interval by themselves. Use set 'expr' to 'NULL' to simulate the
#' work of these functions.
#' 
#' If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
#' set. Use this parameter if the result size exceeds the limits of the
#' physical memory.
#' 
#' @param expr track expression
#' @param intervals genomic scope
#' @param iterator track expression iterator. If 'NULL' iterator is determined
#' implicitly based on track expression.
#' @param band track expression band. If 'NULL' no band is used.
#' @param intervals.set.out intervals set name where the function result is
#' optionally outputed
#' @return If 'intervals.set.out' is 'NULL' a data frame representing iterator
#' intervals.
#' @seealso \code{\link{giterator.cartesian_grid}}
#' @keywords ~iterator ~intervals
#' @examples
#' 
#' gdb.init_examples()
#' 
#' ## iterator is set implicitly to bin size of 'dense' track
#' giterator.intervals("dense_track", gintervals(1, 0, 200))
#' 
#' ## iterator = 30
#' giterator.intervals("dense_track", gintervals(1, 0, 200), 30)
#' 
#' ## iterator is an intervals set named 'annotations'
#' giterator.intervals("dense_track", ALLGENOME, "annotations")
#' 
#' ## iterator is set implicitly to intervals of 'array_track' track
#' giterator.intervals("array_track", gintervals(1, 0, 200))
#' 
#' ## iterator is a rectangle 100000 by 50000
#' giterator.intervals("rects_track",
#'                     gintervals.2d(chroms1 = 1, chroms2 = "chrX"),
#'                     c(100000, 50000))
#' 
#' @export giterator.intervals
giterator.intervals <- function(expr = NULL, intervals = get("ALLGENOME"), iterator = NULL, band = NULL, intervals.set.out = NULL) {
	if (is.null(substitute(expr)) && is.null(substitute(iterator)))
		stop("Usage: giterator.intervals(expr = NULL, intervals = ALLGENOME, iterator = NULL, band = NULL, intervals.set.out = NULL)", call. = F)
	
	if (is.null(substitute(expr)))
		exprstr <- "0"
	else
		exprstr <- do.call(.gexpr2str, list(substitute(expr)), envir = parent.frame())

	.iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())
	intervals.set.out <- do.call(.gexpr2str, list(substitute(intervals.set.out)), envir = parent.frame())

	if (!is.null(intervals.set.out))
		fullpath <- .gintervals.check_new_set(intervals.set.out)

	# intervals can be NULL if gextract is piped with gscreen and the latter returns NULL
	success <- FALSE
	res <- NULL
	tryCatch({
		if (!is.null(intervals)) {
			res <- .gcall("giterator_intervals", exprstr, intervals, .iterator, band, intervals.set.out, new.env(parent = parent.frame()))

			if (!is.null(intervals.set.out) && .gintervals.is_bigset(intervals.set.out, F) && !.gintervals.needs_bigset(intervals.set.out))
				.gintervals.big2small(intervals.set.out)
		}

		success <- TRUE
	},
	finally = {
		if (!success && !is.null(intervals.set.out))
			unlink(fullpath, recursive = TRUE)
	})

	# refresh the list of GINTERVS, etc.
	if (!is.null(intervals.set.out)) {
		.gdb.add_intervals.set(intervals.set.out)
		retv <- 0 # suppress return value
	} else
		res
}



#' Returns DNA sequences
#' 
#' Returns DNA sequences for given intervals
#' 
#' This function returns an array of sequence strings for each interval from
#' 'intervals'. If intervals contain an additional 'strand' column and its
#' value is '-1', the reverse-complementary sequence is returned.
#' 
#' @param intervals intervals for which DNA sequence is returned
#' @return An array of character strings representing DNA sequence.
#' @seealso \code{\link{gextract}}
#' @keywords ~extract ~DNA ~sequence
#' @examples
#' 
#' gdb.init_examples()
#' intervs <- gintervals(c(1, 2), 10000, 10020)
#' gseq.extract(intervs)
#' 
#' @export gseq.extract
gseq.extract <- function(intervals = NULL) {
	if (is.null(intervals))
		stop("Usage: gseq.extract(intervals)", call. = F);
	.gcheckroot()

	res <- .gcall("gseqread", intervals, new.env(parent = parent.frame()))
	res
}



#' Creates a 'Rectangles' track from intervals and values
#' 
#' Creates a 'Rectangles' track from intervals and values.
#' 
#' This function creates a new 'Rectangles' (two-dimensional) track with values
#' at given intervals. 'description' is added as a track attribute.
#' 
#' @param track track name
#' @param description a character string description
#' @param intervals a set of two-dimensional intervals
#' @param values an array of numeric values - one for each interval
#' @return None.
#' @seealso \code{\link{gtrack.create}}, \code{\link{gtrack.create_sparse}},
#' \code{\link{gtrack.smooth}}, \code{\link{gtrack.modify}},
#' \code{\link{gtrack.rm}}, \code{\link{gtrack.info}},
#' \code{\link{gdir.create}}, \code{\link{gtrack.attr.get}}
#' @keywords ~create ~track
#' @examples
#' 
#' gdb.init_examples()
#' intervs1 <- gintervals.2d(1, (1 : 5) * 200, (1 : 5) * 200 + 100,
#'                           1, (1 : 5) * 300, (1 : 5) * 300 + 200)
#' intervs2 <- gintervals.2d("X", (7 : 10) * 100, (7 : 10) * 100 + 50,
#'                           2, (1 : 4) * 200, (1 : 4) * 200 + 130)
#' intervs <- rbind(intervs1, intervs2)
#' gtrack.2d.create("test_rects", "Test 2d track", intervs,
#'                  runif(dim(intervs)[1], 1, 100))
#' gextract("test_rects", ALLGENOME)
#' gtrack.rm("test_rects", force = TRUE)
#' 
#' @export gtrack.2d.create
gtrack.2d.create <- function(track = NULL, description = NULL, intervals = NULL, values = NULL) {
	if (is.null(substitute(track)) || is.null(description) || is.null(intervals) || is.null(values))
		stop("Usage: gtrack.2d.create(track, description, intervals, values)", call. = F)
	.gcheckroot()

	trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
	intervalsstr <- deparse(substitute(intervals), width.cutoff = 500)[1]
	valuesstr <- deparse(substitute(values), width.cutoff = 500)[1]
	trackdir <- sprintf("%s.track", paste(get("GWD"), gsub("\\.", "/", trackstr), sep = "/"))

	direxisted <- file.exists(trackdir)

	if (!is.na(match(trackstr, get("GTRACKS"))))
		stop(sprintf("Track %s already exists", trackstr), call. = F)
	
	.gconfirmtrackcreate(trackstr)
	success <- FALSE
	tryCatch({
		.gcall("gtrack_create_track2d", trackstr, intervals, values, new.env(parent = parent.frame()), silent = TRUE)
		.gdb.add_track(trackstr)
		.gtrack.attr.set(trackstr, "created.by", sprintf("gtrack.2d.create(%s, description, %s, %s)", trackstr, intervalsstr, valuesstr), T)
		.gtrack.attr.set(trackstr, "created.date", date(), T)
		.gtrack.attr.set(trackstr, "description", description, T)
		success <- TRUE
	},
			 finally = {
				 if (!success && !direxisted) {
					 unlink(trackdir, recursive = TRUE)
					 .gdb.rm_track(trackstr)
				 }
			 }
	)
	retv <- 0 # suppress return value
}



#' Creates a 2D track from tab-delimited file
#' 
#' Creates a 2D track from tab-delimited file(s).
#' 
#' This function creates a 2D track track from one or more tab-delimited files.
#' Each file must start with a header desribing the columns. The first 6
#' columns must have the following names: 'chrom1', 'start1', 'end1', 'chrom2',
#' 'start2', 'end2'. The last column is designated for the value and it may
#' have an arbitrary name. The header is followed by a list of intervals and a
#' value for each interval. Overlapping intervals are forbidden.
#' 
#' One can learn about the format of the tab-delimited file by running
#' 'gextract' function on a 2D track with a 'file' parameter set to the name of
#' the file.
#' 
#' If all the imported intervals represent a point (i.e. end == start + 1) a
#' 'Points' track is created otherwise it is a 'Rectangles' track.
#' 
#' 'description' is added as a track attribute.
#' 
#' Note: temporary files are created in the directory of the track during the
#' run of the function. A few of them need to be kept simultaneously open. If
#' the number of chromosomes and / or intervals is particularly high, a few
#' thousands files might be needed to be opened simultaneously. Some operating
#' systems limit the number of open files per user, in which case the function
#' might fail with "Too many open files" or similar error. The workaround could
#' be:
#' 
#' 1. Increase the limit of simultaneously opened files (the way varies
#' depending on your operating system). 2. Increase the value of
#' 'gmax.data.size' option. Higher values of 'gmax.data.size' option will
#' increased memory usage of the function but create fewer temporary files.
#' 
#' @param track track name
#' @param description a character string description
#' @param file vector of file paths
#' @return None.
#' @seealso \code{\link{gtrack.rm}}, \code{\link{gtrack.info}},
#' \code{\link{gdir.create}}
#' @keywords ~track
#' @export gtrack.2d.import
gtrack.2d.import <- function(track = NULL, description = NULL, file = NULL) {
	if (is.null(substitute(track)) || is.null(description) || is.null(file))
		stop("Usage: gtrack.2d.import(track, description, file)", call. = F)

	.gcheckroot()

	trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
	trackdir <- sprintf("%s.track", paste(get("GWD"), gsub("\\.", "/", trackstr), sep = "/"))

	direxisted <- file.exists(trackdir)

	if (!is.na(match(trackstr, get("GTRACKS"))))
		stop(sprintf("Track %s already exists", trackstr), call. = F)
	
	.gconfirmtrackcreate(trackstr)
	retv <- 0
	success <- FALSE

	tryCatch({
		.gcall("gtrack_2d_import", trackstr, file, new.env(parent = parent.frame()), silent = TRUE)
		.gdb.add_track(trackstr)
		.gtrack.attr.set(trackstr, "created.by", 
				 sprintf("gtrack.2d.import(%s, description, c(\"%s\"))", trackstr, paste(file, collapse="\", \"")), T)
		.gtrack.attr.set(trackstr, "created.date", date(), T)
		.gtrack.attr.set(trackstr, "description", description, T)
		success <- TRUE
	},
			 finally = {
				 if (!success && !direxisted) {
					 unlink(trackdir, recursive = TRUE)
					 .gdb.rm_track(trackstr)
				 }
			 }
	 )
	retv <- 0 # suppress return value
}



#' Creates a track from a file of inter-genomic contacts
#' 
#' Creates a track from a file of inter-genomic contacts.
#' 
#' This function creates a 'Points' (two-dimensional) track from contacts
#' files. If 'allow.duplicates' is 'TRUE' duplicated contacts are allowed and
#' summed up, otherwise an error is reported.
#' 
#' Contacts (coord1, coord2) within the same chromosome are automatically
#' doubled to include also '(coord2, coord1)' unless 'coord1' equals to
#' 'coord2'.
#' 
#' Contacts may come in one or more files.
#' 
#' If 'fends' is 'NULL' contacts file is expected to be in "intervals-value"
#' tab-separated format. The file starts with a header defining the column
#' names. The first 6 columns must have the following names: 'chrom1',
#' 'start1', 'end1', 'chrom2', 'start2', 'end2'. The last column is designated
#' for the value and it may have an arbitrary name. The header is followed by a
#' list of intervals and a value for each interval. An interval of form
#' (chrom1, start1, end1, chrom2, start2, end2) is added as a point (X, Y) to
#' the resulted track where X = (start1 + end1) / 2 and Y = (start2 + end2) /
#' 2.
#' 
#' One can see an example of "intervals-value" format by running 'gextract'
#' function on a 2D track with a 'file' parameter set to the name of the file.
#' 
#' If 'fends' is not 'NULL' contacts file is expected to be in "fends-value"
#' tab-separated format. It should start with a header containing at least 3
#' column names 'fend1', 'fend2' and 'count' in arbitrary order followed by
#' lines each defining a contact between two fragment ends.
#' 
#' \tabular{lll}{ COLUMN \tab VALUE \tab DESCRIPTION\cr fend1 \tab Integer \tab
#' ID of the first fragment end \cr fend2 \tab Integer \tab ID of the second
#' fragment end \cr count \tab Numeric \tab Value associated with the contact
#' \cr }
#' 
#' A fragment ends file is also in tab-separated format. It should start with a
#' header containing at least 3 column names 'fend', 'chr' and 'coord' in
#' arbitrary order followed by lines each defining a single fragment end.
#' 
#' \tabular{lll}{ COLUMN \tab VALUE \tab DESCRIPTION\cr fend \tab Unique
#' integer \tab ID of the fragment end \cr chr \tab Chromosome name \tab Can be
#' specified with or without "chr" prefix, like: "X" or "chrX" \cr coord \tab
#' Integer \tab Coordinate\cr }
#' 
#' 'description' is added as a track attribute.
#' 
#' Note: temporary files are created in the directory of the track during the
#' run of the function. A few of them need to be kept simultaneously open. If
#' the number of chromosomes and / or contacts is particularly high, a few
#' thousands files might be needed to be opened simultaneously. Some operating
#' systems limit the number of open files per user, in which case the function
#' might fail with "Too many open files" or similar error. The workaround could
#' be:
#' 
#' 1. Increase the limit of simultaneously opened files (the way varies
#' depending on your operating system). 2. Increase the value of
#' 'gmax.data.size' option. Higher values of 'gmax.data.size' option will
#' increased memory usage of the function but create fewer temporary files.
#' 
#' @param track track name
#' @param description a character string description
#' @param contacts vector of contacts files
#' @param fends name of fragment ends file
#' @param allow.duplicates if 'TRUE' duplicated contacts are allowed
#' @return None.
#' @seealso \code{\link{gtrack.2d.import}}, \code{\link{gtrack.rm}},
#' \code{\link{gtrack.info}}, \code{\link{gdir.create}}
#' @keywords ~contacts ~fragment ~track
#' @export gtrack.2d.import_contacts
gtrack.2d.import_contacts <- function(track = NULL, description = NULL, contacts = NULL, fends = NULL, allow.duplicates = T) {
	if (is.null(substitute(track)) || is.null(description) || is.null(contacts))
		stop("Usage: gtrack.2d.import_contacts(track, description, contacts, fends = NULL, allow.duplicates = TRUE)", call. = F)
	.gcheckroot()

	trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
	trackdir <- sprintf("%s.track", paste(get("GWD"), gsub("\\.", "/", trackstr), sep = "/"))

	direxisted <- file.exists(trackdir)

	if (!is.na(match(trackstr, get("GTRACKS"))))
		stop(sprintf("Track %s already exists", trackstr), call. = F)
	
	.gconfirmtrackcreate(trackstr)
	success <- FALSE
	tryCatch({
		.gcall("gtrack_import_contacts", trackstr, contacts, fends, allow.duplicates, new.env(parent = parent.frame()), silent = TRUE)
		.gdb.add_track(trackstr)
		.gtrack.attr.set(trackstr, "created.by",
		    sprintf("gtrack.2d.import_contacts(\"%s\", description, c(\"%s\"), \"%s\", %s)",
			    trackstr, paste(contacts, collapse = "\", \""), ifelse(is.null(fends), "NULL", fends), allow.duplicates),
			T)
		.gtrack.attr.set(trackstr, "created.date", date(), T)
		.gtrack.attr.set(trackstr, "description", description, T)
		success <- TRUE
	},
			 finally = {
				 if (!success && !direxisted) {
					 unlink(trackdir, recursive = TRUE)
					 .gdb.rm_track(trackstr)
				 }
			 }
	)
	retv <- 0 # suppress return value
}



#' Returns values from 'Array' track
#' 
#' Returns values from 'Array' track.
#' 
#' This function returns the column values of an 'Array' track in the genomic
#' scope specified by 'intervals'. 'slice' parameter determines which columns
#' should appear in the result. The columns can be indicated by their names or
#' their indices. If 'slice' is 'NULL' the values of all track columns are
#' returned.
#' 
#' The order inside the result might not be the same as the order of intervals.
#' An additional column 'intervalID' is added to the return value. Use this
#' column to refer to the index of the original interval from the supplied
#' 'intervals'.
#' 
#' If 'file' parameter is not 'NULL' the result is saved to a tab-delimited
#' text file (without 'intervalID' column) rather than returned to the user.
#' This can be especially useful when the result is too big to fit into the
#' physical memory.  The resulted file can be used as an input for
#' 'gtrack.array.import' function.
#' 
#' If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
#' set. Similarly to 'file' parameter 'intervals.set.out' can be useful to
#' overcome the limits of the physical memory.
#' 
#' @param track track name
#' @param slice a vector of column names or column indices or 'NULL'
#' @param intervals genomic scope for which the function is applied
#' @param file file name where the function result is to be saved. If 'NULL'
#' result is returned to the user.
#' @param intervals.set.out intervals set name where the function result is
#' optionally outputed
#' @return If 'file' and 'intervals.set.out' are 'NULL' a set of intervals with
#' additional columns for 'Array' track column values and 'columnID'.
#' @seealso \code{\link{gextract}}, \code{\link{gtrack.array.get_colnames}},
#' \code{\link{gtrack.array.import}}
#' @keywords ~extract ~array
#' @examples
#' 
#' gdb.init_examples()
#' gtrack.array.extract("array_track", c("col3", "col5"),
#'                      gintervals(1, 0, 2000))
#' 
#' @export gtrack.array.extract
gtrack.array.extract <- function(track = NULL, slice = NULL, intervals = NULL, file = NULL, intervals.set.out = NULL) {
	if (is.null(substitute(track)) || is.null(intervals))
		stop("Usage: gtrack.array.extract(track, slice, intervals, file = NULL, intervals.set.out = NULL)", call. = F)
	.gcheckroot()

	trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
	slice <- .gslice(trackstr, slice)

	intervals.set.out <- do.call(.gexpr2str, list(substitute(intervals.set.out)), envir = parent.frame())
	if (!is.null(intervals.set.out))
		fullpath <- .gintervals.check_new_set(intervals.set.out)

	# intervals can be NULL if the function is piped with gscreen and the latter returns NULL
	success <- FALSE
	res <- NULL
	tryCatch({
		if (!is.null(intervals)) {
			res <- .gcall("garrayextract", trackstr, slice$slice, slice$colnames, file, intervals, intervals.set.out, new.env(parent = parent.frame()))

			if (!is.null(intervals.set.out) && .gintervals.is_bigset(intervals.set.out, F) && !.gintervals.needs_bigset(intervals.set.out))
				.gintervals.big2small(intervals.set.out)
		}

		success <- TRUE
	},
	finally = {
		if (!success && !is.null(intervals.set.out))
			unlink(fullpath, recursive = TRUE)
	})

	# refresh the list of GINTERVS, etc.
	if (!is.null(intervals.set.out)) {
		.gdb.add_intervals.set(intervals.set.out)
		retv <- 0 # suppress return value
	} else if (!is.null(file))
		retv <- 0 # suppress return value
	else
		res
}



#' Returns column names of array track
#' 
#' Returns column names of array track.
#' 
#' This function returns the column names of an array track.
#' 
#' @param track track name
#' @return A character vector with column names.
#' @seealso \code{\link{gtrack.array.set_colnames}},
#' \code{\link{gtrack.array.extract}}, \code{\link{gvtrack.array.slice}},
#' \code{\link{gtrack.info}}
#' @keywords ~array ~columns
#' @examples
#' 
#' gtrack.array.get_colnames("array_track")
#' 
#' @export gtrack.array.get_colnames
gtrack.array.get_colnames <- function(track = NULL) {
	if (is.null(substitute(track)))
		stop("Usage: gtrack.array.get_colnames(track)", call. = F)

	trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
	names(.gtrack.array.get_colnames(trackstr))
}



#' Creates an array track from array tracks or files
#' 
#' Creates an array track from array tracks or files.
#' 
#' This function creates a new 'Array' track from one or more "sources". Each
#' source can be either another 'Array' track or a tab-delimited file that
#' contains one-dimentional intervals and column values that should be added to
#' the newly created track. One can learn about the exact format of the file by
#' running 'gtrack.array.extract' or 'gextract' functions with a 'file'
#' parameter and inspecting the output file.
#' 
#' There might be more than one source used to create the new track. In that
#' case the new track will contain the columns from all the sources. The
#' equally named columns are merged. Intervals that appear in one source but
#' not in the other are added and the values for the missing columns are set to
#' NaN. Intervals with all NaN values are not added. Partial overlaps between
#' two intervals from different sources are forbidden.
#' 
#' 'description' is added as a track attribute.
#' 
#' @param track name of the newly created track
#' @param description a character string description
#' @param src array track or name of a tab-delimited file
#' @return None.
#' @seealso \code{\link{gextract}}, \code{\link{gtrack.array.extract}},
#' \code{\link{gtrack.array.set_colnames}}, \code{\link{gtrack.rm}},
#' \code{\link{gtrack.info}}, \code{\link{gdir.create}}
#' @keywords ~array ~import ~create ~track
#' @examples
#' 
#' f1 <- tempfile()
#' gextract("sparse_track", gintervals(1, 5000, 20000), file = f1)
#' f2 <- tempfile()
#' gtrack.array.extract("array_track", c("col2", "col3", "col4"),
#'                      gintervals(1, 0, 20000), file = f2)
#' f3 <- tempfile()
#' gtrack.array.extract("array_track", c("col1", "col3"),
#'                      gintervals(1, 0, 20000), file = f3)
#' 
#' gtrack.array.import("test_track1", "Test array track 1", f1, f2)
#' gtrack.array.extract("test_track1", NULL, ALLGENOME)
#' 
#' gtrack.array.import("test_track2", "Test array track 2",
#'                     "test_track1", f3)
#' gtrack.array.extract("test_track2", NULL, ALLGENOME)
#' 
#' gtrack.rm("test_track1", TRUE)
#' gtrack.rm("test_track2", TRUE)
#' unlink(c(f1, f2, f3))
#' 
#' @export gtrack.array.import
gtrack.array.import <- function(track = NULL, description = NULL, ...) {
	args <- as.list(substitute(list(...)))[-1L]
	if (is.null(substitute(track)) || is.null(description) || !length(args))
		stop("Usage: gtrack.array.import(track, description, [src]+)", call. = F)
	.gcheckroot()

	trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())

	srcs <- c()
	colnames <- list()
	for (src in args) {
	    src <- do.call(.gexpr2str, list(src), envir = parent.frame())
		srcs <- c(srcs, src)
		if (is.na(match(src, get("GTRACKS"))))
			colnames[[length(colnames) + 1]] <- as.character(NULL)
		else {
			if (.gcall_noninteractive(gtrack.info, src)$type != "array")
				stop(sprintf("Track %s: only array tracks can be used as a source", src), call. = F)
			colnames[[length(colnames) + 1]] <- names(.gtrack.array.get_colnames(src))
		}
	}

	trackdir <- sprintf("%s.track", paste(get("GWD"), gsub("\\.", "/", trackstr), sep = "/"))

	direxisted <- file.exists(trackdir)

	if (!is.na(match(trackstr, get("GTRACKS"))))
		stop(sprintf("Track %s already exists", trackstr), call. = F)
	
	.gconfirmtrackcreate(trackstr)
	success <- FALSE
	tryCatch({
		colnames <- .gcall("garrays_import", trackstr, srcs, colnames, new.env(parent = parent.frame()), silent = TRUE)
		.gdb.add_track(trackstr)
		.gtrack.array.set_colnames(trackstr, colnames, FALSE)
		created.by <- sprintf("gtrack.array.import(\"%s\", description, src = c(\"%s\"))", trackstr, paste(srcs, collapse = "\", \""))
		.gtrack.attr.set(trackstr, "created.by", created.by, T)
		.gtrack.attr.set(trackstr, "created.date", date(), T)
		.gtrack.attr.set(trackstr, "description", description, T)
		success <- TRUE
	},
			 finally = {
				 if (!success && !direxisted) {
					 unlink(trackdir, recursive = TRUE)
					 .gdb.rm_track(trackstr)
				 }
			 }
	)
	retv <- 0 # suppress return value
}



#' Sets column names of array track
#' 
#' Sets column names of array track.
#' 
#' This sets the column names of an array track.
#' 
#' @param track track name
#' @param track vector of column names
#' @return None.
#' @seealso \code{\link{gtrack.array.get_colnames}},
#' \code{\link{gtrack.array.extract}}, \code{\link{gvtrack.array.slice}},
#' \code{\link{gtrack.info}}
#' @keywords ~array ~columns
#' @examples
#' 
#' old.names <- gtrack.array.get_colnames("array_track")
#' new.names <- paste("modified", old.colnames, sep = "_")
#' gtrack.array.set_colnames("array_track", new.names)
#' gtrack.array.get_colnames("array_track")
#' gtrack.array.set_colnames("array_track", old.names)
#' gtrack.array.get_colnames("array_track")
#' 
#' @export gtrack.array.set_colnames
gtrack.array.set_colnames <- function(track = NULL, names = NULL) {
	if (is.null(substitute(track)) || is.null(names))
		stop("Usage: gtrack.array.set_colnames(track, names)", call. = F)

	trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
	.gtrack.array.set_colnames(trackstr, names, TRUE)
}



#' Converts a track to the most current format
#' 
#' Converts a track (if needed) to the most current format.
#' 
#' This function converts a track to the most current format. It should be used
#' if a track created by an old version of the library cannot be read anymore
#' by the newer version. The old track is given by 'src.track'. After
#' conversion a new track 'tgt.track' is created. If 'tgt.track' is 'NULL' the
#' source track is overwritten.
#' 
#' @param src.track source track name
#' @param tgt.track target track name. If 'NULL' the source track is
#' overwritten.
#' @return None
#' @seealso \code{\link{gtrack.create}}, \code{\link{gtrack.2d.create}},
#' \code{\link{gtrack.create_sparse}}
#' @keywords ~track ~convert
#' @export gtrack.convert
gtrack.convert <- function(src.track = NULL, tgt.track = NULL) {
	if (is.null(substitute(src.track)))
		stop("Usage: gtrack.convert(src.track, tgt.track = NULL)", call. = F)
	.gcheckroot()

	src.trackstr <- do.call(.gexpr2str, list(substitute(src.track)), envir = parent.frame())
	if (is.na(match(src.trackstr, get("GTRACKS"))))
		stop(sprintf("Track %s does not exist", src.trackstr), call. = F)
	
	tgt.trackstr <- ""
	if (is.null(substitute(tgt.track))) {
		tgt.trackstr <- paste(src.trackstr, "_converted", sep = "")
		counter <- 2
		while (!is.na(match(tgt.trackstr, get("GTRACKS")))) {
			tgt.trackstr <- paste(src.trackstr, "_converted", counter, sep = "")
			counter <- counter + 1
		}
	} else {
		tgt.trackstr <- do.call(.gexpr2str, list(substitute(tgt.track)), envir = parent.frame())
		.gconfirmtrackcreate(tgt.trackstr)
	}

	src.dirname <- sprintf("%s.track", paste(get("GWD"), gsub("\\.", "/", src.trackstr), sep = "/"))
	tgt.dirname <- sprintf("%s.track", paste(get("GWD"), gsub("\\.", "/", tgt.trackstr), sep = "/"))
	success <- FALSE

	tryCatch({
		.gcall("gtrackconvert", src.trackstr, tgt.trackstr, new.env(parent = parent.frame()), silent = TRUE)

		# copy all supplimentary data of a track (vars, etc.)
		if (!system(sprintf("cp -r -u %s/. %s", src.dirname, tgt.dirname))) {
			# if tgt track is null move it to the source track
			if (is.null(substitute(tgt.track))) {
				unlink(src.dirname, recursive = TRUE)
				success <- TRUE
				file.rename(tgt.dirname, src.dirname)
			}
		} else {
			msg = sprintf("Failed to copy some or all track supplementary data from %s to %s", src.dirname, tgt.dirname)
			if (is.null(substitute(tgt.track))) {
				msg <- paste(msg,
							 sprintf("Track %s will remain unchanged.\nA new converted track named %s was created without supplementary data.",
									 src.trackstr, tgt.trackstr),
							 sep = "\n")
			}
			warning(msg, call. = F)
		}
		success <- TRUE
	},
			 finally = {
				 if (!success)
					 unlink(tgt.dirname, recursive = TRUE)
				 .gdb.rm_track(tgt.trackstr)
			 }
	)
	retv <- 0 # suppress return value
}



#' Creates a track from a track expression
#' 
#' Creates a track from a track expression.
#' 
#' This function creates a new track named track. The values of the track are
#' determined by evaluation of 'expr' - a numeric track expression. The type of
#' the new track is determined by the type of the iterator. 'Fixed bin',
#' 'Sparse' or 'Rectangles' track can be created accordingly. 'description' is
#' added as a track attribute.
#' 
#' @param track track name
#' @param description a character string description
#' @param expr track expression
#' @param iterator track expression iterator. If 'NULL' iterator is determined
#' implicitly based on track expression.
#' @param band track expression band. If 'NULL' no band is used.
#' @return None.
#' @seealso \code{\link{gtrack.2d.create}}, \code{\link{gtrack.create_sparse}},
#' \code{\link{gtrack.smooth}}, \code{\link{gtrack.modify}},
#' \code{\link{gtrack.rm}}, \code{\link{gtrack.info}},
#' \code{\link{gdir.create}}
#' @keywords ~create ~track
#' @examples
#' 
#' gdb.init_examples()
#' 
#' ## Creates a new track that is a sum of values from 'dense' and
#' ## 2 * non-nan values of 'sparse' track. The new track type is
#' ## Dense with a bin size that equals to '70'.
#' gtrack.create("mixed_track", "Test track",
#'               "dense_track +
#'               replace(sparse_track, is.nan(sparse_track), 0) * 2",
#'               iterator = 70)
#' gtrack.info("mixed_track")
#' gtrack.rm("mixed_track", force = T)
#' 
#' @export gtrack.create
gtrack.create <- function(track = NULL, description = NULL, expr = NULL, iterator = NULL, band = NULL) {
	if (is.null(substitute(track)) || is.null(description) || is.null(substitute(expr)))
		stop("Usage: gtrack.create(track, description, expr, iterator = NULL, band = NULL)", call. = F)
	.gcheckroot()

	trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
	exprstr <- do.call(.gexpr2str, list(substitute(expr)), envir = parent.frame())
	.iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())
	trackdir <- sprintf("%s.track", paste(get("GWD"), gsub("\\.", "/", trackstr), sep = "/"))

	direxisted <- file.exists(trackdir)

	if (!is.na(match(trackstr, get("GTRACKS"))))
		stop(sprintf("Track %s already exists", trackstr), call. = F)
	
	.gconfirmtrackcreate(trackstr)
	success <- FALSE
	tryCatch({
		if (.ggetOption("gmultitasking"))
			.gcall("gtrackcreate_multitask", trackstr, exprstr, .iterator, band, new.env(parent = parent.frame()), silent = TRUE)
		else
			.gcall("gtrackcreate", trackstr, exprstr, .iterator, band, new.env(parent = parent.frame()), silent = TRUE)
        .gdb.add_track(trackstr)
		.gtrack.attr.set(trackstr, "created.by",
				 sprintf("gtrack.create(%s, description, %s, iterator=%s)", trackstr, exprstr, deparse(substitute(iterator), width.cutoff = 500)[1]), T)
		.gtrack.attr.set(trackstr, "created.date", date(), T)
		.gtrack.attr.set(trackstr, "description", description, T)
		success <- TRUE
	},
			 finally = {
				 if (!success && !direxisted) {
					 unlink(trackdir, recursive = TRUE)
                     .gdb.rm_track(trackstr)
				 }
			 }
	)
	retv <- 0 # suppress return value
}



#' Creates a new track from PSSM energy function
#' 
#' Creates a new track from PSSM energy function.
#' 
#' This function creates a new track with values of a PSSM energy function.
#' PSSM parameters (nucleotide probability per position and ploarization) are
#' determined by 'pssmset' key and data files ('pssmset.key' and
#' 'pssmset.data'). These two files must be located in 'GROOT/pssms' directory.
#' The type of the created track is determined by the type of the iterator.
#' 'description' is added as a track attribute.
#' 
#' @param track track name
#' @param description a character string description
#' @param pssmset name of PSSM set: 'pssmset.key' and 'pssmset.data' must be
#' presented in 'GROOT/pssms' directory
#' @param pssmid PSSM id
#' @param prior prior
#' @param iterator track expression iterator for the newly created track
#' @return None.
#' @seealso \code{\link{gtrack.create}}, \code{\link{gtrack.2d.create}},
#' \code{\link{gtrack.create_sparse}}, \code{\link{gtrack.smooth}},
#' \code{\link{gtrack.modify}}, \code{\link{gtrack.rm}},
#' \code{\link{gtrack.info}}, \code{\link{gdir.create}}
#' @keywords ~energy ~pssm ~pwm ~track
#' @examples
#' 
#' gdb.init_examples()
#' gtrack.create_pwm_energy("pwm_energy_track", "Test track","pssm",
#'                          3, 0.01, iterator = 100)
#' gextract("pwm_energy_track", gintervals(1, 0, 1000))
#' gtrack.rm("pwm_energy_track", force = TRUE)
#' 
#' @export gtrack.create_pwm_energy
gtrack.create_pwm_energy <- function(track = NULL, description = NULL, pssmset = NULL, pssmid = NULL, prior = NULL, iterator = NULL) {
	if (is.null(substitute(track)) || is.null(description) || is.null(pssmset) || is.null(pssmid) || is.null(prior) || is.null(iterator))
		stop("Usage: gtrack.create_pwm_energy(track, description, pssmset, pssmid, prior, iterator)", call. = F)
	.gcheckroot()

	trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
	trackdir <- sprintf("%s.track", paste(get("GWD"), gsub("\\.", "/", trackstr), sep = "/"))
	direxisted <- file.exists(trackdir)
	.iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())

	if (!is.na(match(trackstr, get("GTRACKS"))))
		stop(sprintf("Track %s already exists", trackstr), call. = F)
	
	.gconfirmtrackcreate(trackstr)
	success <- FALSE
	tryCatch({
		if (.ggetOption("gmultitasking"))
			.gcall("gcreate_pwm_energy_multitask", trackstr, pssmset, pssmid, prior, .iterator, new.env(parent = parent.frame()), silent = TRUE)
		else
			.gcall("gcreate_pwm_energy", trackstr, pssmset, pssmid, prior, .iterator, new.env(parent = parent.frame()), silent = TRUE)
		.gdb.add_track(trackstr)
		.gtrack.attr.set(trackstr, "created.by",
			sprintf("gtrack.create_pwm_energy(%s, description, \"%s\", %g, %g, iterator=%s)",
			trackstr, pssmset, as.numeric(pssmid), as.numeric(prior), deparse(substitute(iterator), width.cutoff = 500)[1]), T)
		.gtrack.attr.set(trackstr, "created.date", date(), T)
		.gtrack.attr.set(trackstr, "description", description, T)
		success <- TRUE
	},
			 finally = {
				 if (!success && !direxisted) {
					 unlink(trackdir, recursive = TRUE)
					 .gdb.rm_track(trackstr)
				 }
			 }
	)
	retv <- 0 # suppress return value
}



#' Creates a 'Sparse' track from intervals and values
#' 
#' Creates a 'Sparse' track from intervals and values.
#' 
#' This function creates a new 'Sparse' track with values at given intervals.
#' 'description' is added as a track attribute.
#' 
#' @param track track name
#' @param description a character string description
#' @param intervals a set of one-dimensional intervals
#' @param values an array of numeric values - one for each interval
#' @return None.
#' @seealso \code{\link{gtrack.create}}, \code{\link{gtrack.2d.create}},
#' \code{\link{gtrack.smooth}}, \code{\link{gtrack.modify}},
#' \code{\link{gtrack.rm}}, \code{\link{gtrack.info}},
#' \code{\link{gdir.create}}
#' @keywords ~create ~sparse ~track
#' @examples
#' 
#' gdb.init_examples()
#' intervs <- gintervals.load("annotations")
#' gtrack.create_sparse("test_sparse", "Test track", intervs,
#'                      1 : dim(intervs)[1])
#' gextract("test_sparse", ALLGENOME)
#' gtrack.rm("test_sparse", force = TRUE)
#' 
#' @export gtrack.create_sparse
gtrack.create_sparse <- function(track = NULL, description = NULL, intervals = NULL, values = NULL) {
	if (is.null(substitute(track)) || is.null(description) || is.null(intervals) || is.null(values))
		stop("Usage: gtrack.create_sparse(track, description, intervals, values)", call. = F)
	.gcheckroot()

	trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
	intervalsstr <- deparse(substitute(intervals), width.cutoff = 500)[1]
	valuesstr <- deparse(substitute(values), width.cutoff = 500)[1]
	trackdir <- sprintf("%s.track", paste(get("GWD"), gsub("\\.", "/", trackstr), sep = "/"))

	direxisted <- file.exists(trackdir)

	if (!is.na(match(trackstr, get("GTRACKS"))))
		stop(sprintf("Track %s already exists", trackstr), call. = F)
	
	.gconfirmtrackcreate(trackstr)
	success <- FALSE
	tryCatch({
		.gcall("gtrack_create_sparse", trackstr, intervals, values, new.env(parent = parent.frame()), silent = TRUE)
		.gdb.add_track(trackstr)
		.gtrack.attr.set(trackstr, "created.by", sprintf("gtrack.create_sparse(%s, description, %s, %s)", trackstr, intervalsstr, valuesstr), T)
		.gtrack.attr.set(trackstr, "created.date", date(), T)
		.gtrack.attr.set(trackstr, "description", description, T)
		success <- TRUE
	},
			 finally = {
				 if (!success && !direxisted) {
					 unlink(trackdir, recursive = TRUE)
					 .gdb.rm_track(trackstr)
				 }
			 }
	)
	retv <- 0 # suppress return value
}



#' Tests for a track existence
#' 
#' Tests for a track existence.
#' 
#' This function returns 'TRUE' if a track exists in Genomic Database.
#' 
#' @param track track name
#' @return 'TRUE' if a track exists. Otherwise 'FALSE'.
#' @seealso \code{\link{gtrack.ls}}, \code{\link{gtrack.info}},
#' \code{\link{gtrack.create}}, \code{\link{gtrack.rm}}
#' @keywords ~track
#' @examples
#' 
#' gdb.init_examples()
#' gtrack.exists("dense_track")
#' 
#' @export gtrack.exists
gtrack.exists <- function(track = NULL) {
	if (is.null(substitute(track)))
		stop("Usage: gtrack.exists(track)", call. = F)
	.gcheckroot()

	trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
	!is.na(match(trackstr, get("GTRACKS")))
}



#' Returns track attributes values
#' 
#' Returns track attributes values.
#' 
#' This function returns a data frame that contains track attributes values.
#' Column names of the data frame consist of the attribute names, row names
#' contain the track names.
#' 
#' The list of required tracks is specified by 'tracks' argument. If 'tracks'
#' is 'NULL' the attribute values of all existing tracks are returned.
#' 
#' Likewise the list of required attributes is controled by 'attrs' argument.
#' If 'attrs' is 'NULL' all attribute values of the specified tracks are
#' returned. The columns are also sorted then by "popularity" of an attribute,
#' i.e. the number of tracks containing this attribute. This sorting is not
#' applied if 'attrs' is not 'NULL'.
#' 
#' Empty character string in a table cell marks a non-existing attribute.
#' 
#' @param tracks a vector of track names or 'NULL'
#' @param attrs a vector of attribute names or 'NULL'
#' @return A data frame containing track attributes values.
#' @seealso \code{\link{gtrack.attr.import}}, \code{\link{gtrack.attr.get}},
#' \code{\link{gtrack.attr.set}}
#' @keywords ~attr ~attribute
#' @examples
#' 
#' gdb.init_examples()
#' gtrack.attr.export()
#' gtrack.attr.export(tracks = c("sparse_track", "dense_track"))
#' gtrack.attr.export(attrs = "created.by")
#' 
#' @export gtrack.attr.export
gtrack.attr.export <- function(tracks = NULL, attrs = NULL) {
	.gcheckroot()

	if (.ggetOption(".ginteractive", FALSE)) {
		trackstr <- do.call(.gexpr2str, list(substitute(tracks)), envir = parent.frame())
		if (!is.na(match(trackstr, get("GTRACKS"))))
			tracks <- trackstr
	}

	if (!is.null(attrs))
	    tracks <- unique(tracks)
	if (!is.null(attrs))
	    attrs <- unique(attrs)
		
	if (is.null(tracks))
	    tracks <- get("GTRACKS")
	else {
    	idx <- which(!(tracks %in% get("GTRACKS")))[1]
	    if (!is.na(idx))
	        stop(sprintf("Track %s does not exist", tracks[idx]), call. = F)
	}

	.gcall("gget_tracks_attrs", tracks, attrs, new.env(parent = parent.frame()))
}



#' Imports track attributes values
#' 
#' Imports track attributes values.
#' 
#' This function makes imports attribute values contained in a data frame
#' 'table'. The format of a table is similar to the one returned by
#' 'gtrack.attr.export'. The values of the table must be character strings.
#' Column names of the table should specify the attribute names, while row
#' names should contain the track names.
#' 
#' The specified attributes of the specified tracks are modified. If an
#' attribute value is an empty string this attribute is removed from the track.
#' 
#' If 'remove.others' is 'TRUE' all non-readonly attributes that do not appear
#' in the table are removed, otherwise they are preserved unchanged.
#' 
#' Error is reported on an attempt to modify a value of a read-only attribute.
#' 
#' @param table a data frame containing attribute values
#' @param remove.others specifies what to do with the attributes that are not
#' in the table
#' @return None.
#' @seealso \code{\link{gtrack.attr.import}}, \code{\link{gtrack.attr.set}},
#' \code{\link{gtrack.attr.get}}, \code{\link{gdb.get_readonly_attrs}}
#' @keywords ~attr ~attribute
#' @examples
#' 
#' gdb.init_examples()
#' t <- gtrack.attr.export()
#' t$newattr <- as.character(1 : dim(t)[1])
#' gtrack.attr.import(t)
#' gtrack.attr.export(attrs = "newattr")
#' 
#' # roll-back the changes
#' t$newattr <- ""
#' gtrack.attr.import(t)
#' 
#' @export gtrack.attr.import
gtrack.attr.import <- function(table = NULL, remove.others = FALSE) {
	if (is.null(table))
	    stop("Usage: gtrack.attr.import(table, remove.others = FALSE)", call. = F)
	.gcheckroot()

	tracks <- rownames(table)
	attrs <- colnames(table)
	
	if (!is.data.frame(table) || any(dim(table) < 1) || !length(tracks) || !length(attrs) || any(is.na(tracks)) || any(is.na(attrs)) || any(attrs == ""))
        stop("Invalid format of attributes table", call. = F)

    idx <- which(!(tracks %in% get("GTRACKS")))[1]
    if (!is.na(idx))
        stop(sprintf("Track %s does not exist", tracks[idx]), call. = F)

	idx <- which(duplicated(tracks))[1]
	if (!is.na(idx))
	    stop(sprintf("Track %s appears more than once", tracks[idx]), call. = F)

	idx <- which(duplicated(attrs))[1]
	if (!is.na(idx))
	    stop(sprintf("Attribute %s appears more than once", attrs[idx]), call. = F)

	.gtrack.attr.import(table, remove.others, FALSE)
	retv <- 0 # suppress return value
}



#' Returns value of a track attribute
#' 
#' Returns value of a track attribute.
#' 
#' This function returns the value of a track attribute. If the attribute does
#' not exist an empty sting is returned.
#' 
#' @param track track name
#' @param attr attribute name
#' @return Track attribute value.
#' @seealso \code{\link{gtrack.attr.import}}, \code{\link{gtrack.attr.set}}
#' @keywords ~attr ~attribute
#' @examples
#' 
#' gdb.init_examples()
#' gtrack.attr.set("sparse_track", "test_attr", "value")
#' gtrack.attr.get("sparse_track", "test_var")
#' gtrack.attr.set("sparse_track", "test_attr", "")
#' 
#' @export gtrack.attr.get
gtrack.attr.get <- function(track = NULL, attr = NULL) {
	if (is.null(substitute(track)) || is.null(attr))
	    stop("Usage: gtrack.attr.get(track, attr)", call. = F)
	.gcheckroot()

	trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
	res <- gtrack.attr.export(trackstr, attr)
	res[1, 1]
}



#' Assigns value to a track attribute
#' 
#' Assigns value to a track attribute.
#' 
#' This function creates a track attribute and assigns 'value' to it. If the
#' attribute already exists its value is overwritten.
#' 
#' If 'value' is an empty string the attribute is removed.
#' 
#' Error is reported on an attempt to modify a value of a read-only attribute.
#' 
#' @param track track name
#' @param attr attribute name
#' @param value value
#' @return None.
#' @seealso \code{\link{gtrack.attr.get}}, \code{\link{gtrack.attr.import}},
#' \code{\link{gtrack.var.set}}, \code{\link{gdb.get_readonly_attrs}}
#' @keywords ~attr ~attribute
#' @examples
#' 
#' gdb.init_examples()
#' gtrack.attr.set("sparse_track", "test_attr", "value")
#' gtrack.attr.get("sparse_track", "test_var")
#' gtrack.attr.set("sparse_track", "test_attr", "")
#' 
#' @export gtrack.attr.set
gtrack.attr.set <- function(track = NULL, attr = NULL, value = NULL) {
	if (is.null(substitute(track)) || is.null(attr) || is.null(value))
	    stop("Usage: gtrack.attr.set(track, attr, value)", call. = F)
	.gcheckroot()

	trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
	.gtrack.attr.set(trackstr, attr, value, FALSE)
	retv <- 0 # suppress return value
}



#' Creates a track from WIG / BigWig / BedGraph / tab-delimited file
#' 
#' Creates a track from WIG / BigWig / BedGraph / tab-delimited file
#' 
#' This function creates a track from WIG / BigWig / BedGraph / tab-delimited
#' file.  One can learn about the format of the tab-delimited file by running
#' 'gextract' function on a 1D track with a 'file' parameter set to the name of
#' the file. Zipped files are supported (file name must have '.gz' or '.zip'
#' suffix).
#' 
#' If 'binsize' is 0 the resulted track is created in 'Sparse' format.
#' Otherwise the 'Dense' format is chosen with a bin size equal to 'binsize'.
#' The values that were not defined in input file file are substituted by
#' 'defval' value.
#' 
#' 'description' is added as a track attribute.
#' 
#' @param track track name
#' @param description a character string description
#' @param file file path
#' @param binsize bin size of the newly created 'Dense' track or '0' for a
#' 'Sparse' track
#' @param defval default track value
#' @return None.
#' @seealso \code{\link{gtrack.import_set}}, \code{\link{gtrack.rm}},
#' \code{\link{gtrack.info}}, \code{\link{gdir.create}}, \code{\link{gextract}}
#' @keywords ~wig ~bigwig ~bedgraph ~track
#' @export gtrack.import
gtrack.import <- function(track = NULL, description = NULL, file = NULL, binsize = NULL, defval = NaN) {
	if (is.null(substitute(track)) || is.null(description) || is.null(file))
		stop("Usage: gtrack.import(track, description, file, binsize, defval = NaN)", call. = F)

	.gcheckroot()

	trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
	trackdir <- sprintf("%s.track", paste(get("GWD"), gsub("\\.", "/", trackstr), sep = "/"))

	direxisted <- file.exists(trackdir)

	if (!is.na(match(trackstr, get("GTRACKS"))))
		stop(sprintf("Track %s already exists", trackstr), call. = F)
	
	.gconfirmtrackcreate(trackstr)
	retv <- 0
	success <- FALSE

	tmp.dirname <- ""
	file.original <- file

	tryCatch({
		report.progress <- FALSE

		if (length(grep("^.+\\.gz$", file, perl = T)) || length(grep("^.+\\.zip$", file, perl = T))) {
			cat("Unzipping...\n")
			report.progress <- TRUE
			tmp.dirname <- tempfile()
			if (!dir.create(tmp.dirname, recursive = T, mode = "0777"))
				stop(sprintf("Failed to create a directory %s", tmp.dirname), call. = F)

			file.noext <- basename(gsub("^(.+)\\.(.+)$", "\\1", file, perl = T))
			file.unzipped <- paste(tmp.dirname, "/", file.noext, sep = "")
			if (system(paste("/bin/sh -c \"gunzip -q -c", file, ">", file.unzipped, "\"")))
				stop(sprintf("Failed to unzip file %s", file), call. = F)
			file <- file.unzipped
		}

		# looks like all bigWig files start with "fc26" in their first two bytes
		if (length(grep("^.+\\.bw$", file, perl = T)) || length(grep("^.+\\.bigWig$", file, perl = T)) ||
			system(sprintf("od -x -N 2 \"%s\"", file), intern = TRUE)[1] == "0000000 fc26")
		{
			cat("Converting from BigWig to WIG...\n")
			report.progress <- TRUE
			if (tmp.dirname == "") {
				tmp.dirname <- tempfile()
				if (!dir.create(tmp.dirname, recursive = T, mode = "0777"))
					stop(sprintf("Failed to create a directory %s", tmp.dirname), call. = F)
			}

			file.noext <- basename(gsub("^(.+)\\.(.+)$", "\\1", file, perl = T))
			file.converted <- paste(tmp.dirname, "/", file.noext, ".wig", sep = "")
			if (system(paste(get(".GLIBDIR"), "/exec/bigWigToWig ", file, " ", file.converted, sep = "")))
				stop("Command failed", call. = F)
			file <- file.converted
		}

		if (report.progress)
			cat("Converting to track...\n")

		.gcall("gtrackimportwig", trackstr, file, binsize, defval, new.env(parent = parent.frame()), silent = TRUE)
		.gdb.add_track(trackstr)
		.gtrack.attr.set(trackstr, "created.by", 
				 sprintf("gtrack.import(%s, description, \"%s\", %d, %g)", trackstr, file.original, binsize, defval), T)
		.gtrack.attr.set(trackstr, "created.date", date(), T)
		.gtrack.attr.set(trackstr, "description", description, T)
		success <- TRUE
	},
			 finally = {
				 if (tmp.dirname != "")
					unlink(tmp.dirname, recursive = TRUE)

				 if (!success && !direxisted) {
					 unlink(trackdir, recursive = TRUE)
					 .gdb.rm_track(trackstr)
				 }
			 }
	 )
	retv <- 0 # suppress return value
}



#' Creates a track from a file of mapped sequences
#' 
#' Creates a track from a file of mapped sequences.
#' 
#' This function creates a track from a file of mapped sequences. The file can
#' be in SAM format or in a general TAB delimited text format where each line
#' describes a single read.
#' 
#' For a SAM file 'cols.order' must be set to 'NULL'.
#' 
#' For a general TAB delimited text format the following columns must be
#' presented in the file: sequence, chromosome, coordinate and strand. The
#' position of these columns should be specified in 'cols.order' argument. The
#' default value of 'cols.order' is an array of (9, 11, 13, 14) meaning that
#' sequence is expected to be found at column number 9, chromosome - at column
#' 11, coordinate - at column 13 and strand - at column 14. The column indices
#' are 1-based, i.e. the first column is referenced by 1. Chromosome needs a
#' prefix 'chr' e.g. 'chr1'. Valid strand values are '+' or 'F' for forward
#' strand and '-' or 'R' for the reverse strand.
#' 
#' Each read at given coordinate can be "expanded" to cover an interval rather
#' than a single point. The length of the interval is controlled by 'pileup'
#' argument. The direction of expansion depends on the strand value. If
#' 'pileup' is '0', no expansion is performed and the read is converted to a
#' single point. The track is created in sparse format. If 'pileup' is greater
#' than zero, the output track is in dense format. 'binsize' controls the bin
#' size of the dense track.
#' 
#' If 'remove.dups' is 'TRUE' the duplicated coordinates are counted only once.
#' 
#' 'description' is added as a track attribute.
#' 
#' 'gtrack.import_mappedseq' returns the statistics of the conversion process.
#' 
#' @param track track name
#' @param description a character string description
#' @param file name of mapped sequences file
#' @param pileup interval expansion
#' @param binsize bin size of a dense track
#' @param cols.order order of sequece, chromosome, coordinate and strand
#' columns in mapped sequences file or NULL if SAM file is used
#' @param remove.dups if 'TRUE' the duplicated coordinates are counted only
#' once.
#' @return A list of conversion process statistics.
#' @seealso \code{\link{gtrack.rm}}, \code{\link{gtrack.info}},
#' \code{\link{gdir.create}}
#' @keywords ~mapped ~sequence ~track
#' @export gtrack.import_mappedseq
gtrack.import_mappedseq <- function(track = NULL, description = NULL, file = NULL, pileup = 0, binsize = -1, cols.order = c(9, 11, 13, 14), remove.dups = TRUE) {
	if (is.null(substitute(track)) || is.null(description) || is.null(file))
		stop("Usage: gtrack.import_mappedseq(track, description, file, pileup = 0, binsize = -1, cols.order = c(9, 11, 13, 14), remove.dups = TRUE)", call. = F)
	.gcheckroot()

	trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
	trackdir <- sprintf("%s.track", paste(get("GWD"), gsub("\\.", "/", trackstr), sep = "/"))

	direxisted <- file.exists(trackdir)

	if (!is.na(match(trackstr, get("GTRACKS"))))
		stop(sprintf("Track %s already exists", trackstr), call. = F)
	
	.gconfirmtrackcreate(trackstr)
	retv <- 0
	success <- FALSE
	tryCatch({
		retv <- .gcall("gtrackimport_mappedseq", trackstr, file, pileup, binsize, cols.order, remove.dups, new.env(parent = parent.frame()), silent = TRUE)
		.gdb.add_track(trackstr)
		.gtrack.attr.set(trackstr, "created.by",
				 sprintf("gtrack.import_mappedseq(%s, description, \"%s\", pileup=%d, binsize=%d, remove.dups=%s)", trackstr, file, pileup, binsize, remove.dups), T)
		.gtrack.attr.set(trackstr, "created.date", date(), T)
		.gtrack.attr.set(trackstr, "description", description, T)
		success <- TRUE
	},
			 finally = {
				 if (!success && !direxisted) {
					 unlink(trackdir, recursive = TRUE)
					 .gdb.rm_track(trackstr)
				 }
			 }
			 )
	retv
}



#' Creates one or more tracks from multiple WIG / BigWig / BedGraph /
#' tab-delimited files on disk or FTP
#' 
#' Creates one or more tracks from WIG / BigWig / BedGraph / tab-delimited
#' files on disk or FTP.
#' 
#' This function is similar to 'gtrack.import' however unlike the latter it can
#' create multiple tracks. Additionaly the files can be fetched from an FTP
#' server.
#' 
#' The files are expected to be in WIG / BigWig / BedGraph / tab-delimited
#' formats. One can learn about the format of the tab-delimited file by running
#' 'gextract' function with a 'file' parameter set to the name of the file.
#' Zipped files are supported (file name must have '.gz' or '.zip' suffix).
#' 
#' Files are specified by 'path' argument. 'path' can be also a URL of an FTP
#' server in the form of 'ftp://[address]/[files]'. If 'path' is a URL, the
#' files are first downloaded from FTP server to a temporary directory and then
#' imported to tracks. The temporary directory is created at 'GROOT/downloads'.
#' 
#' Regardless whether 'path' is file path or to a URL, it can contain
#' wildcards. Hence multiple files can be imported (and downloaded) at once.
#' 
#' If 'binsize' is 0 the resulted tracks are created in 'Sparse' format.
#' Otherwise the 'Dense' format is chosen with a bin size equal to 'binsize'.
#' The values that were not defined in input file file are substituted by
#' 'defval' value.
#' 
#' The name of a each created track is of '[track.prefix][filename]' form,
#' where 'filename' is the name of the WIG file. For example, if 'track.prefix'
#' equals to "wigs."" and an input file name is 'mydata', a track named
#' 'wigs.mydata' is created. If 'track.prefix' is 'NULL' no prefix is appended
#' to the name of the created track.
#' 
#' Existing tracks are not overwritten and no new directories are automatically
#' created.
#' 
#' 'description' is added to the created tracks as an attribute.
#' 
#' 'gtrack.import_set' does not stop if an error occurs while importing a file.
#' It rather continues importing the rest of the files.
#' 
#' 'gtrack.import_set' returns the names of the files that were successfully
#' imported and those that failed.
#' 
#' @param path file path or URL (may contain wildcards)
#' @param description a character string description
#' @param binsize bin size of the newly created 'Dense' track or '0' for a
#' 'Sparse' track
#' @param track.prefix prefix for a track name
#' @param defval default track value
#' @return Names of files that were successfully imported and those that
#' failed.
#' @seealso \code{\link{gtrack.import}}, \code{\link{gwget}},
#' \code{\link{gtrack.rm}}, \code{\link{gtrack.info}},
#' \code{\link{gdir.create}}, \code{\link{gextract}}
#' @keywords ~wig ~bigwig ~bedgraph ~track
#' @export gtrack.import_set
gtrack.import_set <- function(description = NULL, path = NULL, binsize = NULL, track.prefix = NULL, defval = NaN) {
	.gcheckroot()
	
	if (is.null(description) || is.null(path) || is.null(binsize))
		stop("Usage: gtrack.import_set(description, path, binsize, track.prefix = NULL, defval = NaN)", call. = F)

	if (is.null(substitute(track.prefix)))
		track.prefix <- ""
	else
		track.prefix <- do.call(.gexpr2str, list(substitute(track.prefix)), envir = parent.frame())
	
	files <- c()
	tmp.dirname <- ""

	tryCatch({
		tmp.dirname <- tempfile(pattern = "", tmpdir = paste(get("GROOT"), "/downloads", sep = ""))
		if (!dir.create(tmp.dirname, recursive = T, mode = "0777"))
			stop(sprintf("Failed to create a directory %s", tmp.dirname), call. = F)
		protocol <- "ftp://"
		if (substr(path, 1, nchar(protocol)) == protocol) {
		    # ftp
			files <- gwget(path, tmp.dirname)

			if (!length(files))
				stop("No files downloaded. Exiting.", call. = F)

			files <- paste(tmp.dirname, "/", files, sep = "")
		} else
    		# local path
			files <- system(paste("/bin/sh -c \"ls -d -A", path, "\""), intern=T)

		files <- files[!file.info(files)$isdir]
		if (!length(files))
			stop("No files to import. Exiting.", call. = F)

		files.imported <- c()
	
		for (file in files) {
			tryCatch({
				cat(sprintf("Importing file %s\n", file))
				file.noext <- basename(gsub("^([^.]+)(\\..*)*$", "\\1", file, perl = T))
				trackstr <- paste(track.prefix, file.noext, sep = "")

				.gcall_noninteractive(gtrack.import, trackstr, description, file, binsize, defval)
				files.imported <- c(files.imported, file)
				success <- TRUE
			}, error = function(e) {
				msg <- as.character(e)
				if (msg == "Error: Command interrupted!\n")
					stop("Command interrupted!", call. = F)
				else
					cat(sprintf("%s\n", msg))
			})
		}

		files <- basename(files)
		if (length(files.imported))
			files.imported <- basename(files.imported)
		files.failed <- setdiff(files, files.imported)
		res <- new.env()
		if (length(files.failed))
			res$files.failed <- files.failed
		if (length(files.imported))
			res$files.imported <- files.imported
		as.list(res)
	}, finally = {
		unlink(tmp.dirname, recursive = TRUE) }
	)
}



#' Returns information about a track
#' 
#' Returns information about a track.
#' 
#' Returns information about the track (type, dimensions, size in bytes, etc.).
#' The fields in the returned value vary depending on the type of the track.
#' 
#' @param track track name
#' @return A list that contains track properties
#' @seealso \code{\link{gtrack.exists}}, \code{\link{gtrack.ls}}
#' @keywords ~track ~info ~property
#' @examples
#' 
#' gdb.init_examples()
#' gtrack.info("dense_track")
#' gtrack.info("rects_track")
#' 
#' @export gtrack.info
gtrack.info <- function(track = NULL) {
	if (is.null(substitute(track)))
		stop("Usage: gtrack.info(track)", call. = F)
	.gcheckroot()

	trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
	.gcall("gtrackinfo", trackstr, new.env(parent = parent.frame()))
}



#' Imports a track from another assembly
#' 
#' Imports a track from another assembly.
#' 
#' This function imports a track located in 'src.track.dir' of another assembly
#' to the current database. Chain file instructs how the conversion of
#' coordinates should be done. It can be either a name of a chain file or a
#' data frame in the same format as returned by 'gintervals.load_chain'
#' function. The name of the newly created track is specified by 'track'
#' argument and 'description' is added as a track attribute.
#' 
#' @param track name of a created track
#' @param description a character string description
#' @param src.track.dir path to the directory of the source track
#' @param chain name of chain file or data frame as returned by
#' 'gintervals.load_chain'
#' @return None.
#' @seealso \code{\link{gintervals.load_chain}},
#' \code{\link{gintervals.liftover}}
#' @keywords ~track ~liftover ~chain
#' @export gtrack.liftover
gtrack.liftover <- function(track = NULL, description = NULL, src.track.dir = NULL, chain = NULL) {
	if (is.null(substitute(track)) || is.null(description) || is.null(src.track.dir) || is.null(chain))
		stop("Usage: gtrack.liftover(track, description, src.track.dir, chain)", call. = F)
	.gcheckroot()

	trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
	trackdir <- sprintf("%s.track", paste(get("GWD"), gsub("\\.", "/", trackstr), sep = "/"))
	
	direxisted <- file.exists(trackdir)

	if (!is.na(match(trackstr, get("GTRACKS"))))
		stop(sprintf("Track %s already exists", trackstr), call. = F)

	if (is.character(chain))
		chain.intervs <- gintervals.load_chain(chain)
	else
		chain.intervs <- chain

	.gconfirmtrackcreate(trackstr)
	success <- FALSE
	tryCatch({
		.gcall("gtrack_liftover", trackstr, src.track.dir, chain.intervs, new.env(parent = parent.frame()), silent = TRUE)
		.gdb.add_track(trackstr)
		if (is.character(chain))
			.gtrack.attr.set(trackstr, "created.by", sprintf("gtrack.liftover(%s, description, \"%s\", \"%s\")", trackstr, src.track.dir, chain), T)
		else
			.gtrack.attr.set(trackstr, "created.by", sprintf("gtrack.liftover(%s, description, \"%s\", chain)", trackstr, src.track.dir), T)
		.gtrack.attr.set(trackstr, "created.date", date(), T)
		.gtrack.attr.set(trackstr, "description", description, T)
		success <- TRUE
	},
			 finally = {
				 if (!success && !direxisted) {
					 unlink(trackdir, recursive = TRUE)
					 .gdb.rm_track(trackstr)
				 }
			 }
	)
	retv <- 0 # suppress return value
}



#' Creates a new track from a lookup table based on track expression
#' 
#' Evaluates track expression and translates the values into bin indices that
#' are used in turn to retrieve values from a lookup table and create a track.
#' 
#' This function evaluates the track expression for all iterator intervals and
#' translates this value into an index based on the breaks. This index is then
#' used to address the lookup table and create with its values a new track.
#' More than one 'expr'-'breaks' pair can be used. In that case 'lookup_table'
#' is addressed in a multidimensional manner, i.e. 'lookup_table[i1, i2, ...]'.
#' 
#' The range of bins is determined by 'breaks' argument. For example: 'breaks =
#' c(x1, x2, x3, x4)' represents three different intervals (bins): (x1, x2],
#' (x2, x3], (x3, x4].
#' 
#' If 'include.lowest' is 'TRUE' the the lowest value is included in the first
#' interval, i.e. in [x1, x2].
#' 
#' 'force.binning' parameter controls what should be done when the value of
#' 'expr' exceeds the range determined by 'breaks'. If 'force.binning' is
#' 'TRUE' then values smaller than the minimal break will be translated to
#' index 1, and the values exceeding the maximal break will be translated to
#' index 'M-1' where 'M' is the number of breaks. If 'force.binning' is 'FALSE'
#' the out-of-range values will produce 'NaN' values.
#' 
#' Regardless of 'force.binning' value if the value of 'expr' is 'NaN' then the
#' value in the track would be 'NaN' too.
#' 
#' 'description' is added as a track attribute.
#' 
#' @param track track name
#' @param description a character string description
#' @param lookup_table a multi-dimensional array containing the values that are
#' returned by the function
#' @param expr track expression
#' @param breaks breaks that determine the bin
#' @param include.lowest if 'TRUE', the lowest value of the range determined by
#' breaks is included
#' @param force.binning if 'TRUE', the values smaller than the minimal break
#' will be translated to index 1, and the values that exceed the maximal break
#' will be translated to index N-1 where N is the number of breaks. If 'FALSE'
#' the out-of-range values will produce NaN values.
#' @param iterator track expression iterator. If 'NULL' iterator is determined
#' implicitly based on track expressions.
#' @param band track expression band. If 'NULL' no band is used.
#' @return None.
#' @seealso \code{\link{glookup}}, \code{\link{gtrack.2d.create}},
#' \code{\link{gtrack.create_sparse}}, \code{\link{gtrack.smooth}},
#' \code{\link{gtrack.modify}}, \code{\link{gtrack.rm}},
#' \code{\link{gtrack.info}}, \code{\link{gdir.create}}
#' @keywords ~lookup ~track
#' @examples
#' 
#' gdb.init_examples()
#' 
#' ## one-dimensional example
#' breaks1 = seq(0.1, 0.2, length.out = 6)
#' gtrack.lookup("lookup_track", "Test track", 1 : 5, "dense_track",
#'               breaks1)
#' gtrack.rm("lookup_track", force = TRUE)
#' 
#' ## two-dimensional example
#' t <- array(1 : 15, dim = c(5, 3))
#' breaks2 = seq(0.31, 0.37, length.out = 4)
#' gtrack.lookup("lookup_track", "Test track", t, "dense_track",
#'               breaks1, "2 * dense_track", breaks2)
#' gtrack.rm("lookup_track", force = TRUE)
#' 
#' @export gtrack.lookup
gtrack.lookup <- function(track = NULL, description = NULL, lookup_table = NULL, ..., include.lowest = FALSE, force.binning = TRUE, iterator = NULL, band = NULL) {
	args <- as.list(substitute(list(...)))[-1L]
	if (is.null(substitute(track)) || is.null(description) || is.null(lookup_table) || length(args) < 2 || length(args) %% 2 != 0)
		stop("Usage: gtrack.lookup(track, description, lookup_table, [expr, breaks]+, include.lowest = FALSE, force.binning = TRUE, iterator = NULL, band = NULL)", call. = F)
	.gcheckroot()

	trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
	trackdir <- sprintf("%s.track", paste(get("GWD"), gsub("\\.", "/", trackstr), sep = "/"))
	direxisted <- file.exists(trackdir)

	exprs <- c()
	breaks <- list()
	
	for (i in (0 : (length(args) / 2 - 1))) {
		exprs <- append(exprs, do.call(.gexpr2str, list(args[[i * 2 + 1]]), envir = parent.frame()))
		breaks[[length(breaks) + 1]] <- eval.parent(args[[i * 2 + 2]])
	}

	.iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())

	if (!is.na(match(trackstr, get("GTRACKS"))))
		stop(sprintf("Track %s already exists", trackstr), call. = F)
	
	.gconfirmtrackcreate(trackstr)
	success <- FALSE
	tryCatch({
		.gcall("gtrack_bintransform", trackstr, exprs, breaks, include.lowest, force.binning, lookup_table, .iterator, band, new.env(parent = parent.frame()), silent = TRUE)
		.gdb.add_track(trackstr)
		created.by = sprintf("gtrack.lookup(%s, description, lookup_table", trackstr)
		for (i in (1 : length(exprs)))
			created.by = sprintf("%s, %s, breaks%d", created.by, exprs[i], i)
		created.by = sprintf("%s, include.lowest = %s, force.binning = %s)", created.by, include.lowest, force.binning)
		.gtrack.attr.set(trackstr, "created.by", created.by, T)
		.gtrack.attr.set(trackstr, "created.date", date(), T)
		.gtrack.attr.set(trackstr, "description", description, T)
		success <- TRUE
	},
			 finally = {
				 if (!success && !direxisted) {
					 unlink(trackdir, recursive = TRUE)
					 .gdb.rm_track(trackstr)
				 }
			 }
	)
	retv <- 0 # suppress return value
}



#' Returns a list of track names
#' 
#' Returns a list of track names in Genomic Database.
#' 
#' This function returns a list of tracks whose name or track attribute value
#' match a pattern (see 'grep'). If called without any arguments all tracks are
#' returned.
#' 
#' If pattern is specified without a track attribute (i.e. in the form of
#' 'pattern') then filtering is applied to the track names. If pattern is
#' supplied with a track attribute (i.e. in the form of 'name = pattern') then
#' track attribute is matched against the pattern.
#' 
#' Multiple patterns are applied one after another. The resulted list of tracks
#' should match all the patterns.
#' 
#' @param ... these arguments are of either form 'pattern' or 'attribute =
#' pattern'
#' @param ignore.case,perl,fixed,useBytes see 'grep'
#' @return An array that contains the names of tracks that match the supplied
#' patterns.
#' @seealso \code{\link{grep}}, \code{\link{gtrack.exists}},
#' \code{\link{gtrack.create}}, \code{\link{gtrack.rm}}
#' @keywords ~intervals ~ls
#' @examples
#' 
#' gdb.init_examples()
#' 
#' # get all track names
#' gtrack.ls()
#' 
#' # get track names that match the pattern "den*"
#' gtrack.ls("den*")
#' 
#' # get track names whose "created.by" attribute match the pattern
#' # "create_sparse"
#' gtrack.ls(created.by = "create_sparse")
#' 
#' # get track names whose names match the pattern "den*" and whose
#' # "created.by" attribute match the pattern "track"
#' gtrack.ls("den*", created.by = "track")
#' 
#' @export gtrack.ls
gtrack.ls <- function(..., ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE) {
	.gcheckroot()

    args <- as.list(substitute(list(...)))[-1L]
    args <- list(...)

    tracks <- get("GTRACKS")

	if (is.null(tracks) || !length(tracks))
	    return (NULL)

    if (length(args) >= 1) {
	    attrs <- c()
		patterns <- c()
		
	    # first filter out file names (this filtering is faster than filtering by track variable)
        for (i in 1 : length(args)) {
            arg <- as.character(args[[i]])
            if (is.null(names(args)) || names(args)[i] == "")
                tracks <- grep(arg, tracks, value = TRUE, ignore.case = ignore.case, perl = perl, fixed = fixed, useBytes = useBytes)
			else {
			    attrs <- c(attrs, names(args)[i])
				patterns <- c(patterns, as.character(args[[i]]))
			}
		}
		
	    # filter out by attributes
		if (length(attrs)) {
		    attrs_table <- .gcall("gget_tracks_attrs", tracks, attrs, new.env(parent = parent.frame()))
			if (is.null(attrs_table))
			   return (NULL)

		    cols <- colnames(attrs_table)
	        for (i in 1 : length(attrs)) {
		        idx <- which(cols == attrs[i])[1]
				if (!is.na(idx)) {
				    attrs_table <- subset(attrs_table, grepl(patterns[i], attrs_table[, idx], ignore.case = ignore.case, perl = perl, fixed = fixed, useBytes = useBytes))
					if (!nrow(attrs_table))
					   return (NULL)
				}
			}
			tracks <- rownames(attrs_table)
		}
    }

    tracks
}



#' Modifies track contents
#' 
#' Modifies 'Dense' track contents.
#' 
#' This function modifies the contents of a 'Dense' track by the values of
#' 'expr'. 'intervals' argument controls which portion of the track is
#' modified. The iterator policy is set internally to the bin size of the
#' track.
#' 
#' @param track track name
#' @param expr track expression
#' @param intervals genomic scope for which track is modified
#' @return None.
#' @seealso \code{\link{gtrack.create}}, \code{\link{gtrack.rm}}
#' @keywords ~modify ~track
#' @examples
#' 
#' gdb.init_examples()
#' intervs <- gintervals(1, 300, 800)
#' gextract("dense_track", intervs)
#' gtrack.modify("dense_track", "dense_track * 2", intervs)
#' gextract("dense_track", intervs)
#' gtrack.modify("dense_track", "dense_track / 2", intervs)
#' 
#' @export gtrack.modify
gtrack.modify <- function(track = NULL, expr = NULL, intervals = NULL) {
	if (is.null(substitute(track)) || is.null(substitute(expr)))
		stop("Usage: gtrack.modify(track, expr, intervals = ALLGENOME)", call. = F)
	.gcheckroot()

	if (is.null(intervals))
		intervals <- get("ALLGENOME")
	
	trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
	exprstr <- do.call(.gexpr2str, list(substitute(expr)), envir = parent.frame())

	.gcall("gtrack_modify", trackstr, exprstr, intervals, iterator=trackstr, new.env(parent = parent.frame()))
	
	str <- sprintf("gtrack.modify(%s, %s, intervs)", trackstr, exprstr)
    created.by.str <- gtrack.attr.export(trackstr, "created.by")[1, 1]
	if (is.null(created.by.str))
		created.by.str <- str
	else
		created.by.str <- paste(created.by.str, str, sep = "\n");

	.gtrack.attr.set(trackstr, "created.by", created.by.str, T)
	
	retv <- 0 # suppress return value
}



#' Deletes a track
#' 
#' Deletes a track.
#' 
#' This function deletes a track from the Genomic Database. By default
#' 'gtrack.rm' requires the user to interactively confirm the deletion. Set
#' 'force' to 'TRUE' to suppress the user prompt.
#' 
#' @param track track name
#' @param force if 'TRUE', supresses user confirmation of a named track removal
#' @return None.
#' @seealso \code{\link{gtrack.exists}}, \code{\link{gtrack.ls}},
#' \code{\link{gtrack.create}}, \code{\link{gtrack.2d.create}},
#' \code{\link{gtrack.create_sparse}}, \code{\link{gtrack.smooth}}
#' @keywords ~track
#' @examples
#' 
#' gdb.init_examples()
#' gtrack.create("new_track", "Test track", "2 * dense_track")
#' gtrack.exists("new_track")
#' gtrack.rm("new_track", force = T)
#' gtrack.exists("new_track")
#' 
#' @export gtrack.rm
gtrack.rm <- function(track = NULL, force = FALSE) {
	if (is.null(substitute(track)))
		stop("Usage: gtrack.rm(track, force = FALSE)", call. = F)
	.gcheckroot()

	trackname <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
	
	# check whether track appears among GTRACKS
	if (is.na(match(trackname, get("GTRACKS")))) {
		if (force)
			return(invisible())
		stop(sprintf("Track %s does not exist", trackname), call. = F)
	}

	answer <- "N"
	if (force)
		answer <- "Y"
	else {
		str <- sprintf("Are you sure you want to delete track %s (Y/N)? ", trackname)
		cat(str)
		answer <- toupper(readLines(n = 1))
	}
	
	if (answer == "Y" || answer == "YES") {
		dirname <- sprintf("%s.track", paste(get("GWD"), gsub("\\.", "/", trackname), sep = "/"))

		# remove the track
		unlink(dirname, recursive = TRUE)

		if (file.exists(dirname))
			cat(sprintf("Failed to delete track %s\n", trackname))
		else
			# refresh the list of GTRACKS, etc.
			.gdb.rm_track(trackname)
	}
}



#' Creates a new track from smoothed values of track expression
#' 
#' Creates a new track from smoothed values of track expression.
#' 
#' This function creates a new 'Dense' track named 'track'. The values of the
#' track are results of smoothing the values of 'expr'.
#' 
#' Each track value at coordinate 'C' is determined by smoothing non 'NaN'
#' values of 'expr' over the window around 'C'. The window size is controlled
#' by 'winsize' and is given in coordinate units (not in number of bins),
#' defining the total regions to be considered when smoothing (on both sides of
#' the central point). Two different algorithms can be used for smoothing:
#' 
#' "MEAN" - an arithmetic average.
#' 
#' "LINEAR_RAMP" - a weighted arithmetic average, where the weights linearly
#' decrease as the distance from the center of the window increases.
#' 
#' 'weight_thr' determines the function behavior when some of the values in the
#' window are missing or 'NaN' (missing values may occur at the edges of each
#' chromosome when the window covers an area beyond chromosome boundaries).
#' 'weight_thr' sets the weight sum threshold below which smoothing algorithm
#' returns 'NaN' rather than a smoothing value based on non 'NaN' values in the
#' window.
#' 
#' 'smooth_nans' controls what would be the smoothed value if the central value
#' in the window is 'NaN'. If 'smooth_nans' is 'FALSE' then the smoothed value
#' is set to 'NaN' regardless of 'weight_thr' parameter. Otherwise it is
#' calculated normally.
#' 
#' 'description' is added as a track attribute.
#' 
#' Iterator policy must be of "fixed bin" type.
#' 
#' @param track track name
#' @param description a character string description
#' @param expr track expression
#' @param winsize size of smoothing window
#' @param weight_thr smoothing weight threshold
#' @param smooth_nans if 'FALSE' track value is always set to 'NaN' if central
#' window value is 'NaN', otherwise it is calculated from the rest of non 'NaN'
#' values
#' @param alg smoothing algorithm - "MEAN" or "LINEAR_RAMP"
#' @param iterator track expression iterator of 'Fixed bin' type
#' @return None.
#' @seealso \code{\link{gtrack.create}}, \code{\link{gtrack.2d.create}},
#' \code{\link{gtrack.create_sparse}}, \code{\link{gtrack.modify}},
#' \code{\link{gtrack.rm}}, \code{\link{gtrack.info}},
#' \code{\link{gdir.create}}
#' @keywords ~smooth ~track
#' @examples
#' 
#' gdb.init_examples()
#' gtrack.smooth("smoothed_track", "Test track", "dense_track", 500)
#' gextract("dense_track", "smoothed_track", gintervals(1, 0, 1000))
#' gtrack.rm("smoothed_track", force = TRUE)
#' 
#' @export gtrack.smooth
gtrack.smooth <- function(track = NULL, description = NULL, expr = NULL, winsize = NULL, weight_thr = 0, smooth_nans = F, alg = "LINEAR_RAMP", iterator = NULL) {
	if (is.null(substitute(track)) || is.null(description) || is.null(substitute(expr)) || is.null(winsize))
		stop("Usage: gtrack.smooth(track, description, expr, winsize, weight_thr = 0, smooth_nans = FALSE, alg = \"LINEAR_RAMP\" (\"LINEAR_RAMP\" | \"MEAN\"), iterator = NULL)",
			 call. = F)
	.gcheckroot()

	trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
	exprstr <- do.call(.gexpr2str, list(substitute(expr)), envir = parent.frame())
	trackdir <- sprintf("%s.track", paste(get("GWD"), gsub("\\.", "/", trackstr), sep = "/"))
	direxisted <- file.exists(trackdir)
	.iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())
	
	if (!is.na(match(trackstr, get("GTRACKS"))))
		stop(sprintf("Track %s already exists", trackstr), call. = F)
	
	.gconfirmtrackcreate(trackstr)
	success <- FALSE
	tryCatch({
		.gcall("gsmooth", trackstr, exprstr, winsize, weight_thr, smooth_nans, alg, .iterator, new.env(parent = parent.frame()), silent = TRUE)
		.gdb.add_track(trackstr)
		.gtrack.attr.set(trackstr, "created.by",
			sprintf("gtrack.smooth(%s, description, %s, %s, %s, %s, %s)", trackstr, exprstr, as.character(winsize), as.character(weight_thr), as.character(smooth_nans), alg), T)
		.gtrack.attr.set(trackstr, "created.date", date(), T)
		.gtrack.attr.set(trackstr, "description", description, T)
		success <- TRUE
	},
			 finally = {
				 if (!success && !direxisted) {
					 unlink(trackdir, recursive = TRUE)
					 .gdb.rm_track(trackstr)
				 }
			 }
	)
	retv <- 0 # suppress return value
}



#' Returns value of a track variable
#' 
#' Returns value of a track variable.
#' 
#' This function returns the value of a track variable. If the variable does
#' not exist an error is reported.
#' 
#' @param track track name
#' @param var track variable name
#' @return Track variable value.
#' @seealso \code{\link{gtrack.var.set}}, \code{\link{gtrack.var.ls}},
#' \code{\link{gtrack.var.rm}}
#' @keywords ~variable
#' @examples
#' 
#' gdb.init_examples()
#' gtrack.var.set("sparse_track", "test_var", 1 : 10)
#' gtrack.var.get("sparse_track", "test_var")
#' gtrack.var.rm("sparse_track", "test_var")
#' 
#' @export gtrack.var.get
gtrack.var.get <- function(track = NULL, var = NULL) {
	if (is.null(substitute(track)) || is.null(var))
		stop("Usage: gtrack.var.get(track, var)", call. = F)
	.gcheckroot()

	trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
	.gtrack.var.get(trackstr, var)
}



#' Returns a list of track variables for a track
#' 
#' Returns a list of track variables for a track.
#' 
#' This function returns a list of track variables of a track that match the
#' pattern (see 'grep'). If called without any arguments all track variables of
#' a track are returned.
#' 
#' @param track track name
#' @param pattern,ignore.case,perl,fixed,useBytes see 'grep'
#' @return An array that contains the names of track variables.
#' @seealso \code{\link{grep}}, \code{\link{gtrack.var.get}},
#' \code{\link{gtrack.var.set}}, \code{\link{gtrack.var.rm}}
#' @keywords ~variable ~ls
#' @examples
#' 
#' gdb.init_examples()
#' gtrack.var.ls("sparse_track")
#' gtrack.var.set("sparse_track", "test_var1", 1 : 10)
#' gtrack.var.set("sparse_track", "test_var2", "v")
#' gtrack.var.ls("sparse_track")
#' gtrack.var.ls("sparse_track", pattern = "2")
#' gtrack.var.rm("sparse_track", "test_var1")
#' gtrack.var.rm("sparse_track", "test_var2")
#' 
#' @export gtrack.var.ls
gtrack.var.ls <- function(track = NULL, pattern = "", ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE) {
	if (length(substitute(track)) != 1)
		stop("Usage: gtrack.var.ls(track, pattern = \"\", ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE)", call. = F)
	.gcheckroot()

	trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
	
	if (is.na(match(trackstr, get("GTRACKS"))))
		stop(sprintf("Track %s does not exist", trackstr), call. = F)
	
	trackdir <- sprintf("%s.track", paste(get("GWD"), gsub("\\.", "/", trackstr), sep = "/"))
	dirname <- paste(trackdir, "vars", sep = "/")
	options(warn = -1) # disable warnings since dir() on non dir or non existing dir produces warnings
	invisible(files <- dir(dirname))
	options(warn = 0) # restore the warning behavior
	if (length(files) > 0)
		grep(pattern, files, value = TRUE, ignore.case = ignore.case, perl = perl, fixed = fixed, useBytes = useBytes)
	else
		files
}



#' Deletes a track variable
#' 
#' Deletes a track variable.
#' 
#' This function deletes a track variable.
#' 
#' @param track track name
#' @param var track variable name
#' @return None.
#' @seealso \code{\link{gtrack.var.get}}, \code{\link{gtrack.var.set}},
#' \code{\link{gtrack.var.ls}}
#' @keywords ~variable
#' @examples
#' 
#' gdb.init_examples()
#' gtrack.var.set("sparse_track", "test_var1", 1 : 10)
#' gtrack.var.set("sparse_track", "test_var2", "v")
#' gtrack.var.ls("sparse_track")
#' gtrack.var.rm("sparse_track", "test_var1")
#' gtrack.var.rm("sparse_track", "test_var2")
#' gtrack.var.ls("sparse_track")
#' 
#' @export gtrack.var.rm
gtrack.var.rm <- function(track = NULL, var = NULL) {
	if (is.null(substitute(track)) || is.null(var))
		stop("Usage: gtrack.var.rm(track, var)", call. = F)
	.gcheckroot()

	trackname <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
	if (is.na(match(trackname, get("GTRACKS"))))
		stop(sprintf("Track %s does not exist", trackname), call. = F)

	trackdir <- sprintf("%s.track", paste(get("GWD"), gsub("\\.", "/", trackname), sep = "/"))
	filename <- paste(trackdir, "vars", var, sep = "/")
	invisible(file.remove(filename))
}



#' Assigns value to a track variable
#' 
#' Assigns value to a track variable.
#' 
#' This function creates a track variable and assigns 'value' to it. If the
#' track variable already exists its value is overwritten.
#' 
#' @param track track name
#' @param var track variable name
#' @param value value
#' @return None.
#' @seealso \code{\link{gtrack.var.get}}, \code{\link{gtrack.var.ls}},
#' \code{\link{gtrack.var.rm}}
#' @keywords ~variable
#' @examples
#' 
#' gdb.init_examples()
#' gtrack.var.set("sparse_track", "test_var", 1 : 10)
#' gtrack.var.get("sparse_track", "test_var")
#' gtrack.var.rm("sparse_track", "test_var")
#' 
#' @export gtrack.var.set
gtrack.var.set <- function(track = NULL, var = NULL, value = NULL) {
	if (is.null(substitute(track)) || is.null(var) || is.null(value))
		stop("Usage: gtrack.var.set(track, var, value)", call. = F)
	.gcheckroot()

	trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
	.gtrack.var.set(trackstr, var, value)
}



#' Defines rules for a single value calculation of a virtual 'Array' track
#' 
#' Defines how a single value within an interval is achieved for a virtual
#' track based on 'Array' track.
#' 
#' A track (regular or virtual) used in a track expression is expected to
#' return one value for each track interval. 'Array' tracks store multiple
#' values per interval (one for each 'column') and hence if used in a track
#' expression one must define the way of how a single value should be deduced
#' from several ones.
#' 
#' By default if an 'Array' track is used in a track expressions, its interval
#' value would be the average of all column values that are not NaN.
#' 'gvtrack.array.slice' allows to select specific columns and to specify the
#' function applied to their values.
#' 
#' 'slice' parameter allows to choose the columns. Columns can be indicated by
#' their names or their indices. If 'slice' is 'NULL' the non-NaN values of all
#' track columns are used.
#' 
#' 'func' parameter determines the function applied to the columns' values. Use
#' the following table for a reference of all valid functions and parameters
#' combinations:
#' 
#' \emph{func = "avg", params = NULL} \cr Average of columns' values.
#' 
#' \emph{func = "max", params = NULL} \cr Maximum of columns' values.
#' 
#' \emph{func = "min", params = NULL} \cr Minimum of columns' values.
#' 
#' \emph{func = "stdev", params = NULL} \cr Unbiased standard deviation of
#' columns' values.
#' 
#' \emph{func = "sum", params = NULL} \cr Sum of columns' values.
#' 
#' \emph{func = "quantile", params = [Percentile in the range of [0, 1]]} \cr
#' Quantile of columns' values.
#' 
#' @param vtrack virtual track name
#' @param slice a vector of column names or column indices or 'NULL'
#' @param func,params see below
#' @return None.
#' @seealso \code{\link{gvtrack.create}},
#' \code{\link{gtrack.array.get_colnames}}, \code{\link{gtrack.array.extract}}
#' @keywords ~virtual ~array
#' @examples
#' 
#' gdb.init_examples()
#' gvtrack.create("vtrack1", "array_track")
#' gvtrack.array.slice("vtrack1", c("col2", "col4"), "max")
#' gextract("vtrack1", gintervals(1, 0, 1000))
#' 
#' @export gvtrack.array.slice
gvtrack.array.slice <- function(vtrack = NULL, slice = NULL, func = "avg", params = NULL) {
	if (is.null(substitute(vtrack)))
		stop("Usage: gvtrack.array.slice(vtrack, slice = NULL, func = \"avg\", params = NULL)", call. = F);
	.gcheckroot()

	vtrackstr <- do.call(.gvtrack, list(substitute(vtrack)), envir = parent.frame())

	var <- .gvtrack.get(vtrackstr)
	.slice <- list()
	.slice$slice <- .gslice(var$src, slice)$slice
	.slice$func <- func
	.slice$params <- params
	var$slice <- .slice
	.gvtrack.set(vtrackstr, var)

	retv <- NULL
}



#' Creates a new virtual track
#' 
#' Creates a new virtual track.
#' 
#' This function creates a new virtual track named 'vtrack' with the given
#' source, function and parameters. 'src' can be either a track or intervals
#' (1D or 2D). Use the following table for a reference of all valid source,
#' function and parameters combinations:
#' 
#' \emph{src = [Track], func = "avg", params = NULL} \cr Average track value in
#' iterator interval.
#' 
#' \emph{src = [Track], func = "max", params = NULL} \cr Maximal track value in
#' iterator interval.
#' 
#' \emph{src = [Track], func = "min", params = NULL} \cr Minimal track value in
#' iterator interval.
#' 
#' \emph{src = ['Dense' / 'Sparse' / 'Array' track], func = "nearest", params =
#' NULL} \cr Mean track value in iterator interval. If there are no track
#' values covered by an iterator interator (can occur only in 'Sparse' track),
#' the nearest track value is returned.
#' 
#' \emph{src = ['Dense' / 'Sparse' / 'Array' track], func = "stddev", params =
#' NULL} \cr Unbiased standard deviation of track values in iterator interval.
#' 
#' \emph{src = ['Dense' / 'Sparse' / 'Array' track], func = "sum", params =
#' NULL} \cr Sum of track values in iterator interval.
#' 
#' \emph{src = ['Dense' / 'Sparse' / 'Array' track], func = "quantile", params
#' = [Percentile in the range of [0, 1]]} \cr Quantile of track values in
#' iterator interval.
#' 
#' \emph{src = ['Dense' track], func = "global.percentile", params = NULL} \cr
#' Percentile of an average track value in iterator interval relatively to all
#' values of the track.
#' 
#' \emph{src = ['Dense' track], func = "global.percentile.max", params = NULL}
#' \cr Percentile of a maximal track value in iterator interval relatively to
#' all values of the track.
#' 
#' \emph{src = ['Dense' track], func = "global.percentile.min", params = NULL}
#' \cr Percentile of a minimal track value in iterator interval relatively to
#' all values of the track.
#' 
#' \emph{src = [2D track], func = "area", params = NULL} \cr Area covered by
#' iterator interval.
#' 
#' \emph{src = [2D track], func = "weighted.sum", params = NULL} \cr Weighted
#' sum of values where each weight equals to the intersection area between the
#' iterator interval and the rectangle containing the value.
#' 
#' \emph{src = [1D intervals], func = "distance", params = [Minimal distance
#' from center (default: 0)]} \cr Given the center 'C' of the current iterator
#' interval returns 'DC * X/2', where 'DC' is the normalized distance to the
#' center of the interval that contains 'C', and 'X' is the value of the
#' parameter. If no interval contains 'C' the resulted value is 'D + XXX/2'
#' where 'D' is the distance between 'C' and the edge of the closest interval.
#' Distance can be positive or negative depending on the position of the
#' coordinate relative to the interval and the strand (-1 or 1) of the
#' interval. Distance is always positive if 'strand' is '0' or if 'strand'
#' column is missing. Distance is 'NA' if no intervals exist for the current
#' chromosome.
#' 
#' \emph{src = [1D intervals], func = "distance.center", params = NULL} \cr
#' Given the center 'C' of the current iterator interval returns 'NaN' if 'C'
#' is outside of the intervals, otherwise returns the distance between 'C' and
#' the center of the closest interval. Distance can be positive or negative
#' depending on the position of the coordinate relative to the interval and the
#' strand (-1 or 1) of the interval. Distance is always positive if 'strand' is
#' '0' or if 'strand' column is missing.
#' 
#' Once a virtual track is created one can modify its iterator behavior by
#' calling 'gvtrack.iterator' or 'gvtrack.iterator.2d'.
#' 
#' @param vtrack virtual track name
#' @param src source (track or intervals)
#' @param func,params see below
#' @return None.
#' @seealso \code{\link{gvtrack.info}}, \code{\link{gvtrack.iterator}},
#' \code{\link{gvtrack.iterator.2d}}, \code{\link{gvtrack.array.slice}},
#' \code{\link{gvtrack.ls}}, \code{\link{gvtrack.rm}}
#' @keywords ~virtual
#' @examples
#' 
#' gdb.init_examples()
#' 
#' gvtrack.create("vtrack1", "dense_track", "max")
#' gvtrack.create("vtrack2", "dense_track", "quantile", 0.5)
#' gextract("dense_track", "vtrack1", "vtrack2",
#'          gintervals(1, 0, 10000), iterator = 1000)
#' 
#' gvtrack.create("vtrack3", "dense_track", "global.percentile")
#' gvtrack.create("vtrack4", "annotations", "distance")
#' gdist("vtrack3", seq(0, 1, l = 10), "vtrack4", seq(-500, 500, 200))
#' 
#' @export gvtrack.create
gvtrack.create <- function(vtrack = NULL, src = NULL, func = NULL, params = NULL) {
	if (is.null(substitute(vtrack)) || is.null(substitute(src)))
		stop("Usage: gvtrack.create(vtrack, src, func = NULL, params = NULL)", call. = F)
	.gcheckroot()

	vtrackstr <- do.call(.gexpr2str, list(substitute(vtrack)), envir = parent.frame())
	srcstr <- do.call(.gexpr2str, list(substitute(src)), envir = parent.frame())

	if (!is.na(match(vtrackstr, get("GTRACKS"))))
		stop(sprintf("Cannot create virtual track: regular track named %s already exists", vtrackstr), call. = F)
	
	if (!is.na(match(vtrackstr, get("GINTERVS"))))
		stop(sprintf("Cannot create virtual track: intervals named %s already exists", vtrackstr), call. = F)
	
    if (vtrackstr != make.names(vtrackstr))
        stop(sprintf("\"%s\" is not a syntactically valid name for a variable", vtrackstr), call. = F)

	var <- list()
	if (is.character(srcstr) && !is.na(match(srcstr, get("GTRACKS"))))
		var$src <- srcstr
	else
		var$src <- src
	var$func <- func
	var$params <- params

    .gvtrack.set(vtrackstr, var)

	retv <- NULL
}



#' Returns the definition of a virtual track
#' 
#' Returns the definition of a virtual track.
#' 
#' This function returns the internal represenation of a virtual track.
#' 
#' @param vtrack virtual track name
#' @return Internal representation of a virtual track.
#' @seealso \code{\link{gvtrack.create}}
#' @keywords ~virtual
#' @examples
#' 
#' gdb.init_examples()
#' gvtrack.create("vtrack1", "dense_track", "max")
#' gvtrack.info("vtrack1")
#' 
#' @export gvtrack.info
gvtrack.info <- function(vtrack = NULL) {
	if (is.null(substitute(vtrack)))
		stop("Usage: gvtrack.info(vtrack)", call. = F);
	.gcheckroot()

	vtrackstr <- do.call(.gvtrack, list(substitute(vtrack)), envir = parent.frame())
	.gvtrack.get(vtrackstr)
}



#' Defines modification rules for a one-dimensional iterator in a virtual track
#' 
#' Defines modification rules for a one-dimensional iterator in a virtual
#' track.
#' 
#' This function defines modification rules for one-dimensional iterator
#' intervals in a virtual track.
#' 
#' 'dim' converts a 2D iterator interval (chrom1, start1, end1, chrom2, start2,
#' end2) to a 1D interval. If 'dim' is '1' the interval is converted to
#' (chrom1, start1, end1). If 'dim' is '2' the interval is converted to
#' (chrom2, start2, end2). If 1D iterator is used 'dim' must be set to 'NULL'
#' or '0' (meaning: no conversion is made).
#' 
#' Iterator interval's 'start' coordinate is modified by adding 'sshift'.
#' Similarly 'end' coordinate is altered by adding 'eshift'.
#' 
#' @param vtrack virtual track name
#' @param dim use 'NULL' or '0' for 1D iterators. '1' converts 2D iterator to
#' (chrom1, start1, end1) , '2' converts 2D iterator to (chrom2, start2, end2)
#' @param sshift shift of 'start' coordinate
#' @param eshift shift of 'end' coordinate
#' @return None.
#' @seealso \code{\link{gvtrack.create}}, \code{\link{gvtrack.iterator.2d}}
#' @keywords ~virtual
#' @examples
#' 
#' gdb.init_examples()
#' 
#' gvtrack.create("vtrack1", "dense_track")
#' gvtrack.iterator("vtrack1", sshift = 200, eshift = 200)
#' gextract("dense_track", "vtrack1", gintervals(1, 0, 500))
#' 
#' gvtrack.create("vtrack2", "dense_track")
#' gvtrack.iterator("vtrack2", dim = 1)
#' gextract("vtrack2", gintervals.2d(1, 0, 1000, 1, 0, -1),
#'          iterator = "rects_track")
#' 
#' @export gvtrack.iterator
gvtrack.iterator <- function(vtrack = NULL, dim = NULL, sshift = 0, eshift = 0) {
	if (is.null(substitute(vtrack)))
		stop("Usage: gvtrack.iterator(vtrack, dim = NULL, sshift = 0, eshift = 0)", call. = F);
	.gcheckroot()

	vtrackstr <- do.call(.gvtrack, list(substitute(vtrack)), envir = parent.frame())

	var <- .gvtrack.get(vtrackstr)
	itr <- list()
	itr$type <- "1d"
	itr$dim <- dim
	itr$sshift <- sshift
	itr$eshift <- eshift
	var$itr <- itr
	.gvtrack.set(vtrackstr, var)
	retv <- NULL
}



#' Defines modification rules for a two-dimensional iterator in a virtual track
#' 
#' Defines modification rules for a two-dimensional iterator in a virtual
#' track.
#' 
#' This function defines modification rules for one-dimensional iterator
#' intervals in a virtual track.
#' 
#' Iterator interval's 'start1' coordinate is modified by adding 'sshift1'.
#' Similarly 'end1', 'start2', 'end2' coordinates are altered by adding
#' 'eshift1', 'sshift2' and 'eshift2' accordingly.
#' 
#' @param vtrack virtual track name
#' @param sshift1 shift of 'start1' coordinate
#' @param eshift1 shift of 'end1' coordinate
#' @param sshift2 shift of 'start2' coordinate
#' @param eshift2 shift of 'end2' coordinate
#' @return None.
#' @seealso \code{\link{gvtrack.create}}, \code{\link{gvtrack.iterator}}
#' @keywords ~virtual
#' @examples
#' 
#' gdb.init_examples()
#' gvtrack.create("vtrack1", "rects_track")
#' gvtrack.iterator.2d("vtrack1", sshift1 = 1000, eshift1 = 2000)
#' gextract("rects_track", "vtrack1",
#'          gintervals.2d(1, 0, 5000, 2, 0, 5000))
#' 
#' @export gvtrack.iterator.2d
gvtrack.iterator.2d <- function(vtrack = NULL, sshift1 = 0, eshift1 = 0, sshift2 = 0, eshift2 = 0) {
	if (is.null(substitute(vtrack)))
		stop("Usage: gvtrack.iterator.2d(vtrack, sshift1 = 0, eshift1 = 0, sshift2 = 0, eshift2 = 0)", call. = F);
	.gcheckroot()

	vtrackstr <- do.call(.gvtrack, list(substitute(vtrack)), envir = parent.frame())

	var <- .gvtrack.get(vtrackstr)
	itr <- list()
	itr$type <- "2d"
	itr$sshift1 <- sshift1
	itr$eshift1 <- eshift1
	itr$sshift2 <- sshift2
	itr$eshift2 <- eshift2
	var$itr <- itr
	.gvtrack.set(vtrackstr, var)
	retv <- NULL
}



#' Returns a list of virtual track names
#' 
#' Returns a list of virtual track names.
#' 
#' This function returns a list of virtual tracks that exist in current R
#' environment that match the pattern (see 'grep'). If called without any
#' arguments all virtual tracks are returned.
#' 
#' @param pattern,ignore.case,perl,fixed,useBytes see 'grep'
#' @return An array that contains the names of virtual tracks.
#' @seealso \code{\link{grep}}, \code{\link{gvtrack.create}},
#' \code{\link{gvtrack.rm}}
#' @keywords ~virtual ~ls
#' @examples
#' 
#' gdb.init_examples()
#' gvtrack.create("vtrack1", "dense_track", "max")
#' gvtrack.create("vtrack2", "dense_track", "quantile", 0.5)
#' gvtrack.ls()
#' gvtrack.ls(pattern = "*2")
#' 
#' @export gvtrack.ls
gvtrack.ls <- function(pattern = "", ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE) {
	.gcheckroot()

    if (!exists("GVTRACKS"))
        return(NULL)

    gvtracks <- get("GVTRACKS")
    gwds <- names(gvtracks)
	if (!is.list(gvtracks) || (length(gvtracks) && !is.character(gwds)) || length(gvtracks) != length(gwds))
		stop("Invalid format of GVTRACKS variable.\nTo continue working with virtual tracks please remove this variable from the environment.", call. = F)

    gwd <- get("GWD")
    idx <- match(gwd, gwds)
    if (is.na(idx))
        return(NULL)

    vtracks <- gvtracks[[idx]]
    vtracknames <- names(vtracks)
	if (!is.list(vtracks) || (length(vtracks) && !is.character(vtracknames)) || length(vtracks) != length(vtracknames))
		stop("Invalid format of GVTRACKS variable.\nTo continue working with virtual tracks please remove this variable from the environment.", call. = F)

	if (!length(vtracks))
		return(NULL)

	if (pattern != "")
		grep(pattern, vtracknames, value = TRUE, ignore.case = ignore.case, perl = perl, fixed = fixed, useBytes = useBytes)
	else
		vtracknames
}



#' Deletes a virtual track
#' 
#' Deletes a virtual track.
#' 
#' This function deletes a virtual track from current R environment.
#' 
#' @param vtrack virtual track name
#' @return None.
#' @seealso \code{\link{gvtrack.create}}, \code{\link{gvtrack.ls}}
#' @keywords ~virtual
#' @examples
#' 
#' gdb.init_examples()
#' gvtrack.create("vtrack1", "dense_track", "max")
#' gvtrack.create("vtrack2", "dense_track", "quantile", 0.5)
#' gvtrack.ls()
#' gvtrack.rm("vtrack1")
#' gvtrack.ls()
#' 
#' @export gvtrack.rm
gvtrack.rm <- function(vtrack = NULL) {
	if (is.null(substitute(vtrack)))
		stop("Usage: gvtrack.rm(vtrack)", call. = F);
	.gcheckroot()

	vtrackstr <- do.call(.gvtrack, list(substitute(vtrack)), envir = parent.frame())
	.gvtrack.set(vtrackstr, NULL)
	retv <- NULL
}

#' @export
.gdb.convert_attrs <- function() {
	.gcheckroot()

	ro_attrs <- c("created.by", "created.date")
	.gcall_noninteractive(gdb.set_readonly_attrs, ro_attrs)
		
	for (track in GTRACKS) {
		for (attr in ro_attrs) {
			try({
				if (.gcall_noninteractive(.gtrack.var.exists, track, attr)) {
				   .gcall_noninteractive(.gtrack.attr.set, track, attr, as.character(.gtrack.var.get(track, attr))[1], T)
				   .gcall_noninteractive(gtrack.var.rm, track, attr)
				}
			}, silent = TRUE)
		}
		cat(sprintf("%s\n", track))
	}
}

#' @export
.gdb.convert_tracks <- function() {
	.gcheckroot()

	for (track in GTRACKS) {
		try({
			retv <- try(.gcall_noninteractive(gtrack.info, track), silent = TRUE)
			if (inherits(retv, "try-error") & length(grep("obsolete", retv)) > 0) {
				cat(sprintf("Converting track %s\n", track))
				.gcall_noninteractive(gtrack.convert, track)
			}
		}, silent = TRUE)
	}
}

