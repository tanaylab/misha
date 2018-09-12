\name{gintervals.save}
\alias{gintervals.save}
\title{
Creates a named intervals set
}
\description{
Saves intervals to a named intervals set.
}
\usage{
gintervals.save(intervals.set, intervals)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{intervals.set}{name of an intervals set}
  \item{intervals}{intervals}
}
\details{
This function saves 'intervals' as a named intervals set.
}
\value{
None.
}

\seealso{
\code{\link{gintervals.rm}}, \code{\link{gintervals.load}}, \code{\link{gintervals.exists}}, \code{\link{gintervals.ls}}, \code{\link{gintervals}}, \code{\link{gintervals.2d}}
}
\examples{
gdb.init_examples()
intervs <- gintervals(c(1, 2))
gintervals.save("testintervs", intervs)
gintervals.ls()
gintervals.rm("testintervs", force = TRUE)
}
\keyword{ ~intervals }
