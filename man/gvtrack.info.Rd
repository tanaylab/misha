\name{gvtrack.info}
\alias{gvtrack.info}
\title{
Returns the definition of a virtual track
}
\description{
Returns the definition of a virtual track.
}
\usage{
gvtrack.info(vtrack)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vtrack}{virtual track name}
}
\details{
This function returns the internal represenation of a virtual track.
}
\value{
Internal representation of a virtual track.
}

\seealso{
\code{\link{gvtrack.create}}
}
\examples{
gdb.init_examples()
gvtrack.create("vtrack1", "dense_track", "max")
gvtrack.info("vtrack1")
}
\keyword{ ~virtual }

