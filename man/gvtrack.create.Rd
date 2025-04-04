% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vtrack.R
\name{gvtrack.create}
\alias{gvtrack.create}
\title{Creates a new virtual track}
\usage{
gvtrack.create(vtrack = NULL, src = NULL, func = NULL, params = NULL, ...)
}
\arguments{
\item{vtrack}{virtual track name}

\item{src}{source (track/intervals). NULL for PWM functions}

\item{func}{function name (see above)}

\item{params}{function parameters (see above)}

\item{...}{additional PWM parameters}
}
\value{
None.
}
\description{
Creates a new virtual track.
}
\details{
This function creates a new virtual track named 'vtrack' with the given
source, function and parameters. 'src' can be either a track or intervals
(1D or 2D). Use the following table for a reference of all valid source,
function and parameters combinations:

\emph{src = [Track], func = "avg", params = NULL} \cr Average track value in
iterator interval.

\emph{src = [Track], func = "max", params = NULL} \cr Maximal track value in
iterator interval.

\emph{src = [Track], func = "min", params = NULL} \cr Minimal track value in
iterator interval.

\emph{src = ['Dense' / 'Sparse' / 'Array' track], func = "nearest", params =
NULL} \cr Mean track value in iterator interval. If there are no track
values covered by an iterator interator (can occur only in 'Sparse' track),
the nearest track value is returned.

\emph{src = ['Dense' / 'Sparse' / 'Array' track], func = "stddev", params =
NULL} \cr Unbiased standard deviation of track values in iterator interval.

\emph{src = ['Dense' / 'Sparse' / 'Array' track], func = "sum", params =
NULL} \cr Sum of track values in iterator interval.

\emph{src = ['Dense' / 'Sparse' / 'Array' track], func = "quantile", params
= [Percentile in the range of [0, 1]]} \cr Quantile of track values in
iterator interval.

\emph{src = ['Dense' track], func = "global.percentile", params = NULL} \cr
Percentile of an average track value in iterator interval relatively to all
values of the track.

\emph{src = ['Dense' track], func = "global.percentile.max", params = NULL}
\cr Percentile of a maximal track value in iterator interval relatively to
all values of the track.

\emph{src = ['Dense' track], func = "global.percentile.min", params = NULL}
\cr Percentile of a minimal track value in iterator interval relatively to
all values of the track.

\emph{src = [2D track], func = "area", params = NULL} \cr Area covered by
iterator interval.

\emph{src = [2D track], func = "weighted.sum", params = NULL} \cr Weighted
sum of values where each weight equals to the intersection area between the
iterator interval and the rectangle containing the value.

\emph{src = [1D intervals], func = "distance", params = [Minimal distance
from center (default: 0)]} \cr Given the center 'C' of the current iterator
interval returns 'DC * X/2', where 'DC' is the normalized distance to the
center of the interval that contains 'C', and 'X' is the value of the
parameter. If no interval contains 'C' the resulted value is 'D + XXX/2'
where 'D' is the distance between 'C' and the edge of the closest interval.
Distance can be positive or negative depending on the position of the
coordinate relative to the interval and the strand (-1 or 1) of the
interval. Distance is always positive if 'strand' is '0' or if 'strand'
column is missing. Distance is 'NA' if no intervals exist for the current
chromosome.

\emph{src = [1D intervals], func = "distance.center", params = NULL} \cr
Given the center 'C' of the current iterator interval returns 'NaN' if 'C'
is outside of the intervals, otherwise returns the distance between 'C' and
the center of the closest interval. Distance can be positive or negative
depending on the position of the coordinate relative to the interval and the
strand (-1 or 1) of the interval. Distance is always positive if 'strand' is
'0' or if 'strand' column is missing.

\emph{src = [1D intervals], func = "coverage", params = NULL} \cr
For each iterator interval, calculates the fraction of its length that is covered by the
source intervals. Returns a value between 0 and 1. For example, if an iterator interval is [100,200]
and the source intervals cover positions 120-140 and 160-170, the coverage would be 0.3
((20 + 10) / 100 = 0.3). Overlapping source intervals are first unified.

\emph{func = "pwm", params = list(pssm = matrix, bidirect = TRUE,
prior = 0.01, extend = TRUE)} \cr
Calculates total log-likelihood score of DNA sequence against PSSM.
Uses log-sum-exp over all positions. For bidirect=TRUE, scans both
strands. Prior adds pseudocounts, extend=TRUE allows scoring at boundaries.

\emph{func = "pwm.max", params = list(pssm = matrix, bidirect = TRUE,
prior = 0.01, extend = TRUE)} \cr
Returns maximum log-likelihood score of best PSSM match. bidirect=TRUE
checks both strands. Prior adds pseudocounts, extend=TRUE allows boundary
scoring.

\emph{func = "pwm.max.pos", params = list(pssm = matrix, bidirect = TRUE,
prior = 0.01, extend = TRUE)} \cr
Returns 1-based position of best PSSM match.
If bidirect=TRUE, the position would be positive if the best hit was at the
forward strand, and negative if it was at the reverse strand. When strand is
-1 the position is still according to the forward strand, but the hit is at
the end of the match.
Prior adds pseudocounts, extend=TRUE allows boundary scoring.

For all PWM functions:
\itemize{
  \item pssm: Position-specific scoring matrix (A,C,G,T frequencies)
  \item bidirect: If TRUE, scans both strands; if FALSE, forward only
  \item prior: Pseudocount for frequencies (default: 0.01)
  \item extend: If TRUE, computes boundary scores
  \item strand: If 1, scans forward strand; if -1, scans reverse strand.
For strand == 1, the energy (and position of the best match) would be at
the beginning of the match, for strand == -1, the energy (and position of
the best match) would be at the end of the match.
}

PWM parameters are accepted as list or individual parameters (see examples).

\emph{func = "kmer.count", params = list(kmer = "ACGT", extend = TRUE, strand = 0)} \cr
Counts occurrences of the specified kmer in each interval. The extend=TRUE
parameter (default) allows counting kmers that span interval boundaries.
The strand parameter can be 1 (forward strand), -1 (reverse strand), or 0 (both strands).

\emph{func = "kmer.frac", params = list(kmer = "ACGT", extend = TRUE, strand = 0)} \cr
Calculates the fraction of possible positions in each interval that contain
the specified kmer. The extend=TRUE parameter (default) allows counting kmers
that span interval boundaries. The strand parameter can be 1 (forward strand), -1
(reverse strand), or 0 (both strands).

For kmer functions:
\itemize{
  \item kmer: The DNA sequence to count (case-insensitive)
  \item extend: If TRUE, counts kmers that span interval boundaries
  \item strand: If 1, counts kmers on forward strand; if -1, counts kmers on reverse strand. If
 0, counts kmers on both strands. Default is 0.
}

Kmer parameters are accepted as list or individual parameters (see examples).
Note that for palindromic kmers, setting strand to 1 or -1 is recommended to avoid double counting.

Modify iterator behavior with 'gvtrack.iterator' or 'gvtrack.iterator.2d'.
}
\examples{
\dontshow{
options(gmax.processes = 2)
}

gdb.init_examples()

gvtrack.create("vtrack1", "dense_track", "max")
gvtrack.create("vtrack2", "dense_track", "quantile", 0.5)
gextract("dense_track", "vtrack1", "vtrack2",
    gintervals(1, 0, 10000),
    iterator = 1000
)

gvtrack.create("vtrack3", "dense_track", "global.percentile")
gvtrack.create("vtrack4", "annotations", "distance")
gdist(
    "vtrack3", seq(0, 1, l = 10), "vtrack4",
    seq(-500, 500, 200)
)

gvtrack.create("cov", "annotations", "coverage")
gextract("cov", gintervals(1, 0, 1000), iterator = 100)

pssm <- matrix(
    c(
        0.7, 0.1, 0.1, 0.1, # Example PSSM
        0.1, 0.7, 0.1, 0.1,
        0.1, 0.1, 0.7, 0.1,
        0.1, 0.1, 0.7, 0.1,
        0.1, 0.1, 0.7, 0.1,
        0.1, 0.1, 0.7, 0.1
    ),
    ncol = 4, byrow = TRUE
)
colnames(pssm) <- c("A", "C", "G", "T")
gvtrack.create(
    "motif_score", NULL, "pwm",
    list(pssm = pssm, bidirect = TRUE, prior = 0.01)
)
gvtrack.create("max_motif_score", NULL, "pwm.max",
    pssm = pssm, bidirect = TRUE, prior = 0.01
)
gvtrack.create("max_motif_pos", NULL, "pwm.max.pos",
    pssm = pssm
)
gextract(
    c(
        "dense_track", "motif_score", "max_motif_score",
        "max_motif_pos"
    ),
    gintervals(1, 0, 10000),
    iterator = 500
)

# Kmer counting examples
gvtrack.create("cg_count", NULL, "kmer.count", kmer = "CG", strand = 1)
gvtrack.create("cg_frac", NULL, "kmer.frac", kmer = "CG", strand = 1)
gextract(c("cg_count", "cg_frac"), gintervals(1, 0, 10000), iterator = 1000)

gvtrack.create("at_pos", NULL, "kmer.count", kmer = "AT", strand = 1)
gvtrack.create("at_neg", NULL, "kmer.count", kmer = "AT", strand = -1)
gvtrack.create("at_both", NULL, "kmer.count", kmer = "AT", strand = 0)
gextract(c("at_pos", "at_neg", "at_both"), gintervals(1, 0, 10000), iterator = 1000)

# GC content
gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")
gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C")
gextract("g_frac + c_frac", gintervals(1, 0, 10000),
    iterator = 1000,
    colnames = "gc_content"
)
}
\seealso{
\code{\link{gvtrack.info}}, \code{\link{gvtrack.iterator}},
\code{\link{gvtrack.iterator.2d}}, \code{\link{gvtrack.array.slice}},
\code{\link{gvtrack.ls}}, \code{\link{gvtrack.rm}}
}
\keyword{~virtual}
