\name{gvtrack.rm}
\alias{gvtrack.rm}
\title{
Deletes a virtual track
}
\description{
Deletes a virtual track.
}
\usage{
gvtrack.rm(vtrack)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vtrack}{virtual track name}
}
\details{
This function deletes a virtual track from current R environment.
}
\value{
None.
}

\seealso{
\code{\link{gvtrack.create}}, \code{\link{gvtrack.ls}}
}
\examples{
gdb.init_examples()
gvtrack.create("vtrack1", "dense_track", "max")
gvtrack.create("vtrack2", "dense_track", "quantile", 0.5)
gvtrack.ls()
gvtrack.rm("vtrack1")
gvtrack.ls()
}
\keyword{ ~virtual }
