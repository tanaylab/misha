\name{gtrack.var.ls}
\alias{gtrack.var.ls}
\title{
Returns a list of track variables for a track
}
\description{
Returns a list of track variables for a track.
}
\usage{
gtrack.var.ls(track, pattern = "", ignore.case = FALSE,
              perl = FALSE, fixed = FALSE, useBytes = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{track}{track name}
  \item{pattern, ignore.case, perl, fixed, useBytes}{see 'grep'}
}
\details{
This function returns a list of track variables of a track that match the pattern (see 'grep'). If called without any arguments all track variables of a track are returned.
}
\value{
An array that contains the names of track variables.
}

\seealso{
\code{\link{grep}}, \code{\link{gtrack.var.get}}, \code{\link{gtrack.var.set}}, \code{\link{gtrack.var.rm}}
}
\examples{
gdb.init_examples()
gtrack.var.ls("sparse_track")
gtrack.var.set("sparse_track", "test_var1", 1 : 10)
gtrack.var.set("sparse_track", "test_var2", "v")
gtrack.var.ls("sparse_track")
gtrack.var.ls("sparse_track", pattern = "2")
gtrack.var.rm("sparse_track", "test_var1")
gtrack.var.rm("sparse_track", "test_var2")
}
\keyword{ ~variable }
\keyword{ ~ls }

