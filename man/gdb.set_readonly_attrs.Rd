\name{gdb.set_readonly_attrs}
\alias{gdb.set_readonly_attrs}
\title{
Sets read-only track attributes
}
\description{
Sets read-only track attributes.
}
\usage{
gdb.set_readonly_attrs()
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{attrs}{a vector of read-only attributes names or 'NULL'}
}
\details{
This function sets the list of read-only track attributes. The specified attributes may or may not already exist in the tracks.

If 'attrs' is 'NULL' the list of read-only attributes is emptied.
}
\value{
None.
}

\seealso{
\code{\link{gdb.get_readonly_attrs}}, \code{\link{gtrack.attr.get}}, \code{\link{gtrack.attr.set}}
}
\keyword{ ~attr }
\keyword{ ~attribute }

