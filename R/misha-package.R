

#' Toolkit for analysis of genomic data
#'
#' 'misha' package is intended to help users to efficiently analyze genomic
#' data achieved from various experiments.
#'
#' For a complete list of help resources, use \code{library(help = "misha")}.
#'
#' The following options are available for the package. Use 'options' function
#' to alter the value of the options.
#'
#' \tabular{lll}{ NAME \tab DEFAULT \tab DESCRIPTION \cr gmax.data.size \tab
#' 10000000 \tab Maximal number of datums (intervals, ...) in large data sets
#' stored \cr \tab \tab in memory. Prevents excessive memory usage by various
#' functions \cr \tab \tab such as 'gextract', 'gscreen', etc. \cr
#' gbig.intervals.size \tab 1000000 \tab Minimal number of intervals in a big
#' intervals set format \cr gmax.mem.usage \tab 10000000 \tab Maximal memory
#' consumption of all child processes in KB before the limiting algorithm is
#' invoked. \cr gmax.processes \tab 16 \tab Maximal number of processes for
#' multitasking \cr gmax.processes2core \tab 2 \tab Maximal number of processes
#' per CPU core for multitasking \cr gmin.scope4process \tab 10000 \tab Minimal
#' scope range (for 2D: surface) assigned to a \cr \tab \tab process in
#' multitasking mode. \cr gbuf.size \tab 1000 \tab Size of track expression
#' values buffer. \cr gtrack.chunk.size \tab 100000 \tab Chunk size in bytes of
#' a 2D track. If '0' chunk size is unlimited. \cr gtrack.num.chunks \tab 0
#' \tab Maximal number of 2D track chunks simultaneously stored in \cr \tab
#' \tab memory.\cr }
#'
#' More information about the options can be found in 'User manual' of the
#' package.
#'
#'
#' @name misha-package
#' @useDynLib misha
#' @aliases misha-package misha
#' @docType package
#' @keywords package
NULL