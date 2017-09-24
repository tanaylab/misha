\name{gtrack.rm}
\alias{gtrack.rm}
\title{
Deletes a track
}
\description{
Deletes a track.
}
\usage{
gtrack.rm(track, force = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{track}{track name}
  \item{force}{if 'TRUE', supresses user confirmation of a named track removal}
}
\details{
This function deletes a track from the Genomic Database. By default 'gtrack.rm' requires the user to interactively confirm the deletion. Set 'force' to 'TRUE' to suppress the user prompt.
}
\value{
None.
}

\seealso{
\code{\link{gtrack.exists}}, \code{\link{gtrack.ls}}, \code{\link{gtrack.create}}, \code{\link{gtrack.2d.create}}, \code{\link{gtrack.create_sparse}}, \code{\link{gtrack.smooth}}
}
\examples{
gdb.init_examples()
gtrack.create("new_track", "Test track", "2 * dense_track")
gtrack.exists("new_track")
gtrack.rm("new_track", force = T)
gtrack.exists("new_track")
}
\keyword{ ~track }
