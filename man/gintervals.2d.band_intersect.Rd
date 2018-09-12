\name{gintervals.2d.band_intersect}
\alias{gintervals.2d.band_intersect}
\title{
Intersects two-dimenstional intervals with a band
}
\description{
Intersects two-dimenstional intervals with a band.
}
\usage{
gintervals.2d.band_intersect(intervals, band = NULL,
                             intervals.set.out = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{intervals}{two-dimensional intervals}
  \item{band}{track expression band. If 'NULL' no band is used.}
  \item{intervals.set.out}{intervals set name where the function result is optionally outputed}
}
\details{
This function intersects each two-dimensional interval from 'intervals' with 'band'. If the intersection is not empty, the interval is shrunk to the minimal rectangle that contains the band and added to the return value.

If 'intervals.set.out' is not 'NULL' the result is saved as an intervals set. Use this parameter if the result size exceeds the limits of the physical memory.
}
\value{
If 'intervals.set.out' is 'NULL' a data frame representing the intervals.
}

\seealso{
\code{\link{gintervals.2d}}, \code{\link{gintervals.intersect}}
}
\examples{
gdb.init_examples()
gintervals.2d.band_intersect(gintervals.2d(1), c(10000, 20000))
}
\keyword{ ~band }
\keyword{ ~intersect }
