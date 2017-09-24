\name{gtrack.var.rm}
\alias{gtrack.var.rm}
\title{
Deletes a track variable
}
\description{
Deletes a track variable.
}
\usage{
gtrack.var.rm(track, var)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{track}{track name}
  \item{var}{track variable name}
}
\details{
This function deletes a track variable.
}
\value{
None.
}

\seealso{
\code{\link{gtrack.var.get}}, \code{\link{gtrack.var.set}}, \code{\link{gtrack.var.ls}}
}
\examples{
gdb.init_examples()
gtrack.var.set("sparse_track", "test_var1", 1 : 10)
gtrack.var.set("sparse_track", "test_var2", "v")
gtrack.var.ls("sparse_track")
gtrack.var.rm("sparse_track", "test_var1")
gtrack.var.rm("sparse_track", "test_var2")
gtrack.var.ls("sparse_track")
}
\keyword{ ~variable }

