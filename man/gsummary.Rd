\name{gsummary}
\alias{gsummary}
\title{
Calculates summary statistics of track expression
}
\description{
Calculates summary statistics of track expression.
}
\usage{
gsummary(expr, intervals = ALLGENOME, iterator = NULL,
         band = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expr}{track expression}
  \item{intervals}{genomic scope for which the function is applied}
  \item{iterator}{track expression iterator. If 'NULL' iterator is determined implicitly based on track expression.}
  \item{band}{track expression band. If 'NULL' no band is used.}
}
\details{
This function returns summary statistics of a track expression: total number of bins, total number of bins whose value is NaN, min, max, sum, mean and standard deviation of the values.
}
\value{
An array that represents summary statistics.
}

\seealso{
\code{\link{gintervals.summary}}, \code{\link{gbins.summary}}
}
\examples{
gdb.init_examples()
gsummary("rects_track")
}
\keyword{ ~summary }
\keyword{ ~statistics }
