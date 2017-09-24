\name{gdb.get_readonly_attrs}
\alias{gdb.get_readonly_attrs}
\title{
Returns a list of read-only track attributes
}
\description{
Returns a list of read-only track attributes.
}
\usage{
gdb.get_readonly_attrs()
}
\details{
This function returns a list of read-only track attributes. These attributes are not allowed to be modified or deleted.

If no attributes are marked as read-only a 'NULL' is returned.
}
\value{
A list of read-only track attributes.
}

\seealso{
\code{\link{gdb.set_readonly_attrs}}, \code{\link{gtrack.attr.get}}, \code{\link{gtrack.attr.set}}
}
\keyword{ ~attr }
\keyword{ ~attribute }

