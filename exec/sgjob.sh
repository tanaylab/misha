#!/bin/bash

commandidx=$1
dirname=$2
R=$3

$R  --silent --vanilla --slave <<EOF

retv <- try({        
	library("misha")	
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

