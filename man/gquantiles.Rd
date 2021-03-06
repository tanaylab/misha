\name{gquantiles}
\alias{gquantiles}
\title{
Calculates quantiles of a track expression
}
\description{
Calculates the quantiles of a track expression for the given percentiles.
}
\usage{
gquantiles(expr, percentiles = 0.5, intervals = ALLGENOME,
           iterator = NULL, band = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expr}{track expression}
  \item{percentiles}{an array of percentiles of quantiles in [0, 1] range}
  \item{intervals}{genomic scope for which the function is applied}
  \item{iterator}{track expression iterator. If 'NULL' iterator is determined implicitly based on track expression.}
  \item{band}{track expression band. If 'NULL' no band is used.}
}
\details{
This function calculates the quantiles for the given percentiles.

If data size exceeds the limit (see: 'getOption(gmax.data.size)'), the data is randomly sampled to fit the limit. A warning message is generated. The seed of the pseudo-random generator can be controled through 'grnd.seed' option.

Note: this function is capable to run in multitasking mode. Sampling may vary according to the extent of multitasking. Since multitasking depends on the number of available CPU cores, running the function on two different machines might give different results. Please switch off multitasking if you want to achieve identical results on any machine. For more information regarding multitasking please refer "User Manual".
}
\value{
An array that represent quantiles.
}

\seealso{
\code{\link{gbins.quantiles}}, \code{\link{gintervals.quantiles}}, \code{\link{gdist}}
}
\examples{
gdb.init_examples()
gquantiles("dense_track", c(0.1, 0.6, 0.8), gintervals(c(1, 2)))
}
\keyword{ ~quantiles }
\keyword{ ~percentiles }
