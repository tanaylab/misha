\name{gtrack.attr.set}
\alias{gtrack.attr.set}
\title{
Assigns value to a track attribute
}
\description{
Assigns value to a track attribute.
}
\usage{
gtrack.attr.set(track, attr, value)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{track}{track name}
  \item{attr}{attribute name}
  \item{value}{value}
}
\details{
This function creates a track attribute and assigns 'value' to it. If the attribute already exists its value is overwritten.

If 'value' is an empty string the attribute is removed.

Error is reported on an attempt to modify a value of a read-only attribute.
}
\value{
None.
}

\seealso{
\code{\link{gtrack.attr.get}}, \code{\link{gtrack.attr.import}}, \code{\link{gtrack.var.set}}, \code{\link{gdb.get_readonly_attrs}}
}
\examples{
gdb.init_examples()
gtrack.attr.set("sparse_track", "test_attr", "value")
gtrack.attr.get("sparse_track", "test_var")
gtrack.attr.set("sparse_track", "test_attr", "")
}
\keyword{ ~attr }
\keyword{ ~attribute }

