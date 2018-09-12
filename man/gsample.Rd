\name{gsample}
\alias{gsample}
\title{
Returns samples from the values of track expression
}
\description{
Returns a sample of the specified size from the values of track expression.
}
\usage{
gsample(expr, n, intervals = ALLGENOME, iterator = NULL,
        band = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expr}{track expression}
  \item{n}{a number of items to choose}
  \item{intervals}{genomic scope for which the function is applied}
  \item{iterator}{track expression iterator. If 'NULL' iterator is determined implicitly based on track expression.}
  \item{band}{track expression band. If 'NULL' no band is used.}
}
\details{
This function returns a sample of the specified size from the values of track expression. If 'n' is less than the total number of values, the data is randomally sampled. The seed of the pseudo-random generator can be controled through 'grnd.seed' option.

If 'n' is higher than the total number of values, all values are returned (yet reshuffled).
}
\value{
An array that represent quantiles.
}

\seealso{
\code{\link{gextract}}
}
\examples{
gdb.init_examples()
gsample("sparse_track", 10)
}
\keyword{ ~sample }
