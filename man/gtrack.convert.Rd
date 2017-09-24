\name{gtrack.convert}
\alias{gtrack.convert}
\title{
Converts a track to the most current format
}
\description{
Converts a track (if needed) to the most current format.
}
\usage{
gtrack.convert(src.track, tgt.track = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{src.track}{source track name}
  \item{tgt.track}{target track name. If 'NULL' the source track is overwritten.}
}
\details{
This function converts a track to the most current format. It should be used if a track created by an old version of the library cannot be read anymore by the newer version. The old track is given by 'src.track'. After conversion a new track 'tgt.track' is created. If 'tgt.track' is 'NULL' the source track is overwritten.
}
\value{
None
}

\seealso{
\code{\link{gtrack.create}}, \code{\link{gtrack.2d.create}}, \code{\link{gtrack.create_sparse}}
}
\keyword{ ~track }
\keyword{ ~convert }
