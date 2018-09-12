\name{gtrack.array.get_colnames}
\alias{gtrack.array.get_colnames}
\title{
Returns column names of array track
}
\description{
Returns column names of array track.
}
\usage{
gtrack.array.get_colnames(track)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{track}{track name}
}
\details{
This function returns the column names of an array track.
}
\value{
A character vector with column names. 
}

\seealso{
\code{\link{gtrack.array.set_colnames}}, \code{\link{gtrack.array.extract}}, \code{\link{gvtrack.array.slice}}, \code{\link{gtrack.info}}
}
\examples{
gtrack.array.get_colnames("array_track")
}
\keyword{ ~array }
\keyword{ ~columns }
