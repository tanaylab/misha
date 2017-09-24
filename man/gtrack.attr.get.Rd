\name{gtrack.attr.get}
\alias{gtrack.attr.get}
\title{
Returns value of a track attribute
}
\description{
Returns value of a track attribute.
}
\usage{
gtrack.attr.get(track, attr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{track}{track name}
  \item{attr}{attribute name}
}
\details{
This function returns the value of a track attribute. If the attribute does not exist an empty sting is returned.
}
\value{
Track attribute value.
}

\seealso{
\code{\link{gtrack.attr.import}}, \code{\link{gtrack.attr.set}}
}
\examples{
gdb.init_examples()
gtrack.attr.set("sparse_track", "test_attr", "value")
gtrack.attr.get("sparse_track", "test_var")
gtrack.attr.set("sparse_track", "test_attr", "")
}
\keyword{ ~attr }
\keyword{ ~attribute }

