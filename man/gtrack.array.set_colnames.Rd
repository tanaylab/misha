\name{gtrack.array.set_colnames}
\alias{gtrack.array.set_colnames}
\title{
Sets column names of array track
}
\description{
Sets column names of array track.
}
\usage{
gtrack.array.get_colnames(track, names)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{track}{track name}
  \item{track}{vector of column names}
}
\details{
This sets the column names of an array track.
}
\value{
None. 
}

\seealso{
\code{\link{gtrack.array.get_colnames}}, \code{\link{gtrack.array.extract}}, \code{\link{gvtrack.array.slice}}, \code{\link{gtrack.info}}
}
\examples{
old.names <- gtrack.array.get_colnames("array_track")
new.names <- paste("modified", old.colnames, sep = "_")
gtrack.array.set_colnames("array_track", new.names)
gtrack.array.get_colnames("array_track")
gtrack.array.set_colnames("array_track", old.names)
gtrack.array.get_colnames("array_track")
}
\keyword{ ~array }
\keyword{ ~columns }
