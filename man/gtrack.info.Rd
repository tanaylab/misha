\name{gtrack.info}
\alias{gtrack.info}
\title{
Returns information about a track
}
\description{
Returns information about a track.
}
\usage{
gtrack.info(track)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{track}{track name}
}
\details{
Returns information about the track (type, dimensions, size in bytes, etc.). The fields in the returned value vary depending on the type of the track.
}
\value{
A list that contains track properties
}

\seealso{
\code{\link{gtrack.exists}}, \code{\link{gtrack.ls}}
}
\examples{
gdb.init_examples()
gtrack.info("dense_track")
gtrack.info("rects_track")
}
\keyword{ ~track }
\keyword{ ~info }
\keyword{ ~property }
