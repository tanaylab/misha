\name{gintervals.liftover}
\alias{gintervals.liftover}
\title{
Converts intervals from another assembly
}
\description{
Converts intervals from another assembly to the current one.
}
\usage{
gintervals.liftover(intervals, chain)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{intervals}{intervals from another assembly}
  \item{chain}{name of chain file or data frame as returned by 'gintervals.load_chain'}
}
\details{
This function converts 'intervals' from another assembly to the current one. Chain file instructs how the conversion of coordinates should be done. It can be either a name of a chain file or a data frame in the same format as returned by 'gintervals.load_chain' function.

The converted intervals are returned. An additional column named 'intervalID' is added to the resulted data frame. For each interval in the resulted intervals it indicates the index of the original interval.}
\value{
A data frame representing the converted intervals.
}

\seealso{
\code{\link{gintervals.load_chain}}, \code{\link{gtrack.liftover}}, \code{\link{gintervals}}
}
\examples{
gdb.init_examples()
chainfile <- paste(GROOT, "data/test.chain", sep = "/")
intervs <- data.frame(chrom = "chr25", start = c(0, 7000),
                      end = c(6000, 20000))
gintervals.liftover(intervs, chainfile)
}
\keyword{ ~intervals }
\keyword{ ~liftover }
\keyword{ ~chain }
