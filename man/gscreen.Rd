\name{gscreen}
\alias{gscreen}
\title{
Finds intervals that match track expression
}
\description{
Finds all intervals where track expression is 'TRUE'.
}
\usage{
gscreen(expr, intervals = ALLGENOME, iterator = NULL, band = NULL,
        intervals.set.out = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expr}{logical track expression}
  \item{intervals}{genomic scope for which the function is applied}
  \item{iterator}{track expression iterator. If 'NULL' iterator is determined implicitly based on track expression.}
  \item{band}{track expression band. If 'NULL' no band is used.}
  \item{intervals.set.out}{intervals set name where the function result is optionally outputed}
}
\details{
This function finds all intervals where track expression's value is 'TRUE'.

If 'intervals.set.out' is not 'NULL' the result is saved as an intervals set. Use this parameter if the result size exceeds the limits of the physical memory.
}
\value{
If 'intervals.set.out' is 'NULL' a set of intervals that match track expression.
}

\seealso{
\code{\link{gsegment}}, \code{\link{gextract}}
}
\examples{
gdb.init_examples()
gscreen("dense_track > 0.2 & sparse_track < 0.4",
        iterator = "dense_track")
}
\keyword{ ~screen }
\keyword{ ~interval }
\keyword{ ~intervals }
