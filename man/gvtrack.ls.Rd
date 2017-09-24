\name{gvtrack.ls}
\alias{gvtrack.ls}
\title{
Returns a list of virtual track names
}
\description{
Returns a list of virtual track names.
}
\usage{
gvtrack.ls(pattern = "", ignore.case = FALSE, perl = FALSE,
           fixed = FALSE, useBytes = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pattern, ignore.case, perl, fixed, useBytes}{see 'grep'}
}
\details{
This function returns a list of virtual tracks that exist in current R environment that match the pattern (see 'grep'). If called without any arguments all virtual tracks are returned.
}
\value{
An array that contains the names of virtual tracks.
}

\seealso{
\code{\link{grep}}, \code{\link{gvtrack.create}}, \code{\link{gvtrack.rm}}
}
\examples{
gdb.init_examples()
gvtrack.create("vtrack1", "dense_track", "max")
gvtrack.create("vtrack2", "dense_track", "quantile", 0.5)
gvtrack.ls()
gvtrack.ls(pattern = "*2")
}
\keyword{ ~virtual }
\keyword{ ~ls }
