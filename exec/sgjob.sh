#!/bin/bash

commandidx=$1
dirname=$2
R=$3

$R  --silent --vanilla --slave <<EOF

retv <- try({        
	# the caller's misha (its .GLIBDIR is saved with .misha), unless the caller ran a source tree
	# (devtools::load_all) or its library is not visible here
	local({
		saved <- new.env()
		load(paste("${dirname}", "misha", sep="/"), envir = saved)
		pkgdir <- get(".GLIBDIR", envir = get(".misha", envir = saved))
		lib <- if (file.exists(file.path(pkgdir, "Meta", "package.rds"))) dirname(pkgdir)
		library("misha", lib.loc = c(lib, .libPaths()))
	})
	if (!exists(".gcluster.restore_db", envir = asNamespace("misha")))
		stop("the job loaded misha ", packageVersion("misha"), " from ", find.package("misha"), "; gcluster.run needs misha >= 5.12.1 on the nodes")
    load(paste("${dirname}", "opts", sep="/"))
    options(opts)
	options(echo = FALSE)
    remove(opts)
    misha:::.gcluster.restore_db(paste("${dirname}", "misha", sep="/"))
    load(paste("${dirname}", "envir", sep="/"))
    load(paste("${dirname}", "commands", sep="/"))
    eval(.GSGECMD[[${commandidx}]])
})
save(retv, file = paste("${dirname}", "${commandidx}.retv", sep="/"))

EOF

