\name{gbins.quantiles}
\alias{gbins.quantiles}
\title{
Calculates quantiles of a track expression for bins
}
\description{
Calculates quantiles of a track expression for bins.
}
\usage{
gbins.quantiles([bin_expr, breaks]+, expr, percentiles = 0.5,
                intervals = ALLGENOME, include.lowest = FALSE,
                iterator = NULL, band = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bin_expr}{a track expression that determines the bin}
  \item{breaks}{breaks that define the bins}
  \item{expr}{track expression for which quantiles are calculated}
  \item{percentiles}{an array of percentiles of quantiles in [0, 1] range}
  \item{intervals}{genomic scope for which the function is applied.}
  \item{include.lowest}{if 'TRUE', the lowest value of the range determined by breaks is included}
  \item{iterator}{track expression iterator. If 'NULL' iterator is determined implicitly based on track expressions.}
  \item{band}{track expression band. If 'NULL' no band is used.}
}
\details{
This function is a binned version of 'gquantiles'. For each iterator interval the value of 'bin_expr' is calculated and assigned to the corresponding bin determined by 'breaks'. The quantiles of 'expr' are calculated then separatedly for each bin.

The bins can be multi-dimensional depending on the number of 'bin_expr'-'breaks' pairs.

The range of bins is determined by 'breaks' argument. For example: 'breaks=c(x1, x2, x3, x4)' represents three different intervals (bins): (x1, x2], (x2, x3], (x3, x4].

If 'include.lowest' is 'TRUE' the the lowest value will be included in the first interval, i.e. in [x1, x2].
}
\value{
Multi-dimensional array representing quantiles for each percentile and bin.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{gquantiles}}, \code{\link{gintervals.quantiles}}, \code{\link{gdist}}
}
\examples{
gdb.init_examples()
gbins.quantiles("dense_track", c(0, 0.2, 0.4, 2), "sparse_track",
                percentiles = c(0.2, 0.5),
                intervals = gintervals(1),
                iterator = "dense_track")
}
\keyword{ ~quantiles }
\keyword{ ~percentiles }
