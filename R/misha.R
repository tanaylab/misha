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

.gintervals.is1d <- function(intervals) {
	if (is.character(intervals)) {
		if (.gintervals.is_bigset(intervals))
			return(.gintervals.big.is1d(intervals))
		intervals <- .gintervals.load(intervals)
	}
	all(colnames(intervals)[1 : 3] == c("chrom", "start", "end"))
}

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

gcompute_strands_autocorr <- function(file = NULL, chrom = NULL, binsize = NULL, maxread = 400, cols.order = c(9, 11, 13, 14), min.coord = 0, max.coord = 3e+8) {
	if (is.null(file) || is.null(chrom) || is.null(binsize))
		stop("Usage: gcompute_strands_autocorr(file, chrom, binsize, maxread = 400, cols.order = c(9, 11, 13, 14), min.coord = 0, max.coord = 3e+8)", call. = F)
	.gcheckroot()

	if (substr(chrom, 1, 3) != "chr")
		chrom <- paste("chr", chrom, sep = "")
	
	res <- .gcall("gcompute_strands_autocorr", file, chrom, binsize, maxread, cols.order, min.coord, max.coord, new.env(parent = parent.frame()))
	res
}

gdb.create <- function(groot = NULL, fasta = NULL, genes.file = NULL, annots.file = NULL, annots.names = NULL) {
	if (is.null(groot) || is.null(fasta))
		stop("Usage: gdb.create(groot, fasta, genes.file = NULL, annots.file = NULL, annots.names = NULL)", call. = F);

	if (file.exists(groot))
		stop(sprintf("Directory %s already exists", groot), call. = F)

	success <- FALSE
    allgenome.old <- NULL
    groot.old <- NULL
    if (exists("ALLGENOME"))
	    allgenome.old <- get("ALLGENOME")
    if (exists("GROOT"))
        groot.old <- get("GROOT")

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
    	assign("GROOT", groot, envir = .GlobalEnv)

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
		assign("GROOT", groot.old, envir = .GlobalEnv)
		if (!success)
			unlink(groot, recursive = TRUE)
	})
	retv <- 0 # suppress return value
}

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

gdb.init <- function(groot = NULL, dir = NULL, rescan = FALSE) {
	gsetroot(groot, dir, rescan)
}

gdb.init_examples <- function() {
	gsetroot(paste(.GLIBDIR, "trackdb/test", sep = "/"))
}

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

gdir.cwd <- function() {
	.gcheckroot()
	get("GWD")
}

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

gintervals <- function(chroms = NULL, starts = 0, ends = -1, strands = NULL) {
	if (is.null(chroms))
		stop("Usage: gintervals(chroms, starts = 0, ends = -1, strands = NULL)", call. = F)
	.gcheckroot()

	intervals <- .gintervals(chroms, starts, ends, strands)
	.gcall("gintervsort", intervals, new.env(parent = parent.frame()))
}

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

gintervals.2d.all <- function() {
	.gcheckroot()
	get("ALLGENOME")[[2]]
}

gintervals.all <- function() {
	.gcheckroot()
	get("ALLGENOME")[[1]]
}

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

gintervals.canonic <- function(intervals = NULL, unify_touching_intervals = TRUE) {
	if (is.null(intervals))
		stop("Usage: gintervals.canonic(intervals, unify_touching_intervals = TRUE)", call. = F)

	res <- .gcall("gintervcanonic", intervals, unify_touching_intervals, new.env(parent = parent.frame()))
	res
}

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

gintervals.exists <- function(intervals.set = NULL) {
	if (is.null(substitute(intervals.set)))
		stop("Usage: gintervals.exists(intervals.set)", call. = F)
	.gcheckroot()

	intervals.set <- do.call(.gexpr2str, list(substitute(intervals.set)), envir = parent.frame())
	!is.na(match(intervals.set, get("GINTERVS")))
}

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

gintervals.is.bigset <- function(intervals.set = NULL) {
	if (is.null(intervals.set))
		stop("Usage: gintervals.is.bigset(intervals.set)", call. = F)
	.gcheckroot()

	.gintervals.is_bigset(intervals.set) && !.gintervals.loadable(intervals.set)
}

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

gintervals.load <- function(intervals.set = NULL, chrom = NULL, chrom1 = NULL, chrom2 = NULL) {
	.gintervals.load_ext(intervals.set, chrom, chrom1, chrom2, TRUE)
}

gintervals.load_chain <- function(file = NULL) {
	if (is.null(file))
		stop("Usage: gintervals.load_chain(file)", call. = F)
	.gcall("gchain2interv", file, new.env(parent = parent.frame()))
}

gintervals.ls <- function(pattern = "", ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE) {
	.gcheckroot()
	grep(pattern, get("GINTERVS"), value = TRUE, ignore.case = ignore.case, perl = perl, fixed = fixed, useBytes = useBytes)
}

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

gintervals.save <- function(intervals.set.out = NULL, intervals = NULL) {
	if (is.null(substitute(intervals.set.out)) || is.null(intervals))
		stop("Usage: gintervals.save(intervals.set.out, intervals)", call. = F)
	.gcheckroot()
	
	intervals.set.out <- do.call(.gexpr2str, list(substitute(intervals.set.out)), envir = parent.frame())
	.gintervals.apply(gintervals.chrom_sizes(intervals), intervals, intervals.set.out, function(intervs, ...) { intervs[[1]] })
	retv <- NULL
}

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

gseq.extract <- function(intervals = NULL) {
	if (is.null(intervals))
		stop("Usage: gseq.extract(intervals)", call. = F);
	.gcheckroot()

	res <- .gcall("gseqread", intervals, new.env(parent = parent.frame()))
	res
}

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

gtrack.array.get_colnames <- function(track = NULL) {
	if (is.null(substitute(track)))
		stop("Usage: gtrack.array.get_colnames(track)", call. = F)

	trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
	names(.gtrack.array.get_colnames(trackstr))
}

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

gtrack.array.set_colnames <- function(track = NULL, names = NULL) {
	if (is.null(substitute(track)) || is.null(names))
		stop("Usage: gtrack.array.set_colnames(track, names)", call. = F)

	trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
	.gtrack.array.set_colnames(trackstr, names, TRUE)
}

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

gtrack.exists <- function(track = NULL) {
	if (is.null(substitute(track)))
		stop("Usage: gtrack.exists(track)", call. = F)
	.gcheckroot()

	trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
	!is.na(match(trackstr, get("GTRACKS")))
}

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

gtrack.attr.get <- function(track = NULL, attr = NULL) {
	if (is.null(substitute(track)) || is.null(attr))
	    stop("Usage: gtrack.attr.get(track, attr)", call. = F)
	.gcheckroot()

	trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
	res <- gtrack.attr.export(trackstr, attr)
	res[1, 1]
}

gtrack.attr.set <- function(track = NULL, attr = NULL, value = NULL) {
	if (is.null(substitute(track)) || is.null(attr) || is.null(value))
	    stop("Usage: gtrack.attr.set(track, attr, value)", call. = F)
	.gcheckroot()

	trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
	.gtrack.attr.set(trackstr, attr, value, FALSE)
	retv <- 0 # suppress return value
}

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

gtrack.info <- function(track = NULL) {
	if (is.null(substitute(track)))
		stop("Usage: gtrack.info(track)", call. = F)
	.gcheckroot()

	trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
	.gcall("gtrackinfo", trackstr, new.env(parent = parent.frame()))
}

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

gtrack.var.get <- function(track = NULL, var = NULL) {
	if (is.null(substitute(track)) || is.null(var))
		stop("Usage: gtrack.var.get(track, var)", call. = F)
	.gcheckroot()

	trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
	.gtrack.var.get(trackstr, var)
}

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

gtrack.var.set <- function(track = NULL, var = NULL, value = NULL) {
	if (is.null(substitute(track)) || is.null(var) || is.null(value))
		stop("Usage: gtrack.var.set(track, var, value)", call. = F)
	.gcheckroot()

	trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())
	.gtrack.var.set(trackstr, var, value)
}

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

gvtrack.info <- function(vtrack = NULL) {
	if (is.null(substitute(vtrack)))
		stop("Usage: gvtrack.info(vtrack)", call. = F);
	.gcheckroot()

	vtrackstr <- do.call(.gvtrack, list(substitute(vtrack)), envir = parent.frame())
	.gvtrack.get(vtrackstr)
}

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

gvtrack.rm <- function(vtrack = NULL) {
	if (is.null(substitute(vtrack)))
		stop("Usage: gvtrack.rm(vtrack)", call. = F);
	.gcheckroot()

	vtrackstr <- do.call(.gvtrack, list(substitute(vtrack)), envir = parent.frame())
	.gvtrack.set(vtrackstr, NULL)
	retv <- NULL
}

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

