\name{gintervals.quantiles}
\alias{gintervals.quantiles}
\title{
Calculates quantiles of a track expression for intervals
}
\description{
Calculates quantiles of a track expression for intervals.
}
\usage{
gintervals.quantiles(expr, percentiles = 0.5, intervals,
                     iterator = NULL, band = NULL,
                     intervals.set.out = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expr}{track expression for which quantiles are calculated}
  \item{percentiles}{an array of percentiles of quantiles in [0, 1] range}
  \item{intervals}{set of intervals}
  \item{iterator}{track expression iterator. If 'NULL' iterator is determined implicitly based on track expressions.}
  \item{band}{track expression band. If 'NULL' no band is used.}
  \item{intervals.set.out}{intervals set name where the function result is optionally outputed}
}
\details{
This function calculates quantiles of 'expr' for each interval in 'intervals'.

If 'intervals.set.out' is not 'NULL' the result is saved as an intervals set. Use this parameter if the result size exceeds the limits of the physical memory.
}
\value{
If 'intervals.set.out' is 'NULL' a set of intervals with additional columns representing quantiles for each percentile.
}

\seealso{
\code{\link{gquantiles}}, \code{\link{gbins.quantiles}}
}
\examples{
gdb.init_examples()
intervs <- gintervals(c(1, 2), 0, 5000)
gintervals.quantiles("dense_track",
                     percentiles = c(0.5, 0.3, 0.9), intervs)
}
\keyword{ ~quantiles }
\keyword{ ~percentiles }
