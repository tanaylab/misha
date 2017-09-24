\name{gtrack.var.set}
\alias{gtrack.var.set}
\title{
Assigns value to a track variable
}
\description{
Assigns value to a track variable.
}
\usage{
gtrack.var.set(track, var, value)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{track}{track name}
  \item{var}{track variable name}
  \item{value}{value}
}
\details{
This function creates a track variable and assigns 'value' to it. If the track variable already exists its value is overwritten.
}
\value{
None.
}

\seealso{
\code{\link{gtrack.var.get}}, \code{\link{gtrack.var.ls}}, \code{\link{gtrack.var.rm}}
}
\examples{
gdb.init_examples()
gtrack.var.set("sparse_track", "test_var", 1 : 10)
gtrack.var.get("sparse_track", "test_var")
gtrack.var.rm("sparse_track", "test_var")
}
\keyword{ ~variable }

