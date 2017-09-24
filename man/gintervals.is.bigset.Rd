\name{gintervals.is.bigset}
\alias{gintervals.is.bigset}
\title{
Tests for big intervals set
}
\description{
Tests for big intervals set.
}
\usage{
gintervals.is.bigset(intervals.set)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{intervals.set}{name of an intervals set}
}
\details{
This function tests whether 'intervals.set' is a big intervals set. Intervals set is big if it is stored in big intervals set format and given the current limits it cannot be fully loaded into memory.

Memory limit is controlled by 'gmax.data.size' option (see: 'getOption("gmax.data.size")').
}
\value{
'TRUE' if intervals set is big, otherwise 'FALSE'.
}

\seealso{
\code{\link{gintervals.load}}, \code{\link{gintervals.save}}, \code{\link{gintervals.exists}}, \code{\link{gintervals.ls}}
}
\examples{
gdb.init_examples()
gintervals.is.bigset("annotations")
}
\keyword{ ~intervals }
