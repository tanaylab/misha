\name{gtrack.exists}
\alias{gtrack.exists}
\title{
Tests for a track existence
}
\description{
Tests for a track existence.
}
\usage{
gtrack.exists(track)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{track}{track name}
}
\details{
This function returns 'TRUE' if a track exists in Genomic Database.
}
\value{
'TRUE' if a track exists. Otherwise 'FALSE'.
}

\seealso{
\code{\link{gtrack.ls}}, \code{\link{gtrack.info}}, \code{\link{gtrack.create}}, \code{\link{gtrack.rm}}
}
\examples{
gdb.init_examples()
gtrack.exists("dense_track")
}
\keyword{ ~track }
