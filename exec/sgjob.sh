#!/bin/bash

commandidx=$1
dirname=$2
R=$3

$R  --silent --vanilla --slave <<EOF

retv <- try({
    load(paste("${dirname}", "libdir", sep="/"))
	library("misha", lib.loc=dirname(.GLIBDIR))
	remove(.GLIBDIR)
    load(paste("${dirname}", "opts", sep="/"))
    options(opts)
	options(echo = FALSE)
    remove(opts)
    load(paste("${dirname}", "envir", sep="/"))
    load(paste("${dirname}", "commands", sep="/"))
    eval(.GSGECMD[[${commandidx}]])
})
save(retv, file = paste("${dirname}", "${commandidx}.retv", sep="/"))

EOF

