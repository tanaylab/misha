\name{gintervals.force_range}
\alias{gintervals.force_range}
\title{
Limits intervals to chromosomal range
}
\description{
Limits intervals to chromosomal range.
}
\usage{
gintervals.force_range(intervals, intervals.set.out = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{intervals}{intervals}
  \item{intervals.set.out}{intervals set name where the function result is optionally outputed}
}
\details{
This function enforces the intervals to be within the chromosomal range [0, chrom length) by altering the intervals' boundaries. Intervals that lay entirely outside of the chromosomal range are eliminated. The new intervals are returned.

If 'intervals.set.out' is not 'NULL' the result is saved as an intervals set. Use this parameter if the result size exceeds the limits of the physical memory.
}
\value{
If 'intervals.set.out' is 'NULL' a data frame representing the intervals.
}

\seealso{
\code{\link{gintervals}}, \code{\link{gintervals.2d}}, \code{\link{gintervals.canonic}}
}
\examples{
gdb.init_examples()
intervs <- data.frame(chrom = "chr1",
                      start = c(11000, -100, 10000, 10500),
                      end = c(12000, 200, 13000000, 10600))
gintervals.force_range(intervs)
}
\keyword{ ~intervals }
