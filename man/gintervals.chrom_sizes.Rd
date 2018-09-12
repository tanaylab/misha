\name{gintervals.chrom_sizes}
\alias{gintervals.chrom_sizes}
\title{
Returns number of intervals per chromosome
}
\description{
Returns number of intervals per chromosome (or chromosome pair).
}
\usage{
gintervals.chrom_sizes(intervals.set)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{intervals}{intervals set}
}
\details{
This function returns number of intervals per chromosome (for 1D intervals) or chromosome pair (for 2D intervals).
}
\value{
Data frame representing number of intervals per chromosome (for 1D intervals) or chromosome pair (for 2D intervals).
}

\seealso{
\code{\link{gintervals.load}}, \code{\link{gintervals.save}}, \code{\link{gintervals.exists}}, \code{\link{gintervals.ls}}, \code{\link{gintervals}}, \code{\link{gintervals.2d}}
}
\examples{
gdb.init_examples()
gintervals.chrom_sizes("annotations")
}
\keyword{ ~intervals }

