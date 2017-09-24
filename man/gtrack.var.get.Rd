\name{gtrack.var.get}
\alias{gtrack.var.get}
\title{
Returns value of a track variable
}
\description{
Returns value of a track variable.
}
\usage{
gtrack.var.get(track, var)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{track}{track name}
  \item{var}{track variable name}
}
\details{
This function returns the value of a track variable. If the variable does not exist an error is reported.
}
\value{
Track variable value.
}

\seealso{
\code{\link{gtrack.var.set}}, \code{\link{gtrack.var.ls}}, \code{\link{gtrack.var.rm}}
}
\examples{
gdb.init_examples()
gtrack.var.set("sparse_track", "test_var", 1 : 10)
gtrack.var.get("sparse_track", "test_var")
gtrack.var.rm("sparse_track", "test_var")
}
\keyword{ ~variable }

