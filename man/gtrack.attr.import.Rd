\name{gtrack.attr.import}
\alias{gtrack.attr.import}
\title{
Imports track attributes values
}
\description{
Imports track attributes values.
}
\usage{
gtrack.attr.import(table, remove.others = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{table}{a data frame containing attribute values}
  \item{remove.others}{specifies what to do with the attributes that are not in the table}
}
\details{
This function makes imports attribute values contained in a data frame 'table'. The format of a table is similar to the one returned by 'gtrack.attr.export'. The values of the table must be character strings. Column names of the table should specify the attribute names, while row names should contain the track names.

The specified attributes of the specified tracks are modified. If an attribute value is an empty string this attribute is removed from the track.

If 'remove.others' is 'TRUE' all non-readonly attributes that do not appear in the table are removed, otherwise they are preserved unchanged.

Error is reported on an attempt to modify a value of a read-only attribute.
}
\value{
None.
}

\seealso{
\code{\link{gtrack.attr.import}}, \code{\link{gtrack.attr.set}}, \code{\link{gtrack.attr.get}}, \code{\link{gdb.get_readonly_attrs}}
}
\examples{
gdb.init_examples()
t <- gtrack.attr.export()
t$newattr <- as.character(1 : dim(t)[1])
gtrack.attr.import(t)
gtrack.attr.export(attrs = "newattr")

# roll-back the changes
t$newattr <- ""
gtrack.attr.import(t)
}
\keyword{ ~attr }
\keyword{ ~attribute }

