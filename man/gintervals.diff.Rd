\name{gintervals.diff}
\alias{gintervals.diff}
\title{
Calculates difference of two intervals sets
}
\description{
Returns difference of two sets of intervals.
}
\usage{
gintervals.diff(intervals1, intervals2, intervals.set.out = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{intervals1, intervals2}{set of one-dimensional intervals}
  \item{intervals.set.out}{intervals set name where the function result is optionally outputed}
}
\details{
This function returns a genomic space that is covered by 'intervals1' but not covered by 'intervals2'.

If 'intervals.set.out' is not 'NULL' the result is saved as an intervals set. Use this parameter if the result size exceeds the limits of the physical memory.
}
\value{
If 'intervals.set.out' is 'NULL' a data frame representing the intervals.
}

\seealso{
\code{\link{gintervals}}, \code{\link{gintervals.intersect}}, \code{\link{gintervals.union}}
}
\examples{
gdb.init_examples()

intervs1 <- gscreen("dense_track > 0.15")
intervs2 <- gscreen("dense_track < 0.2")

## 'res3' equals to 'res4'
res3 <- gintervals.diff(intervs1, intervs2)
res4 <- gscreen("dense_track >= 0.2")
}
\keyword{ ~diff }
