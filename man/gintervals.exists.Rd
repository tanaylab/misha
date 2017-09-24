\name{gintervals.exists}
\alias{gintervals.exists}
\title{
Tests for a named intervals set existence
}
\description{
Tests for a named intervals set existence.
}
\usage{
gintervals.exists(intervals.set)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{intervals.set}{name of an intervals set}
}
\details{
This function returns 'TRUE' if a named intervals set exists in Genomic Database.
}
\value{
'TRUE' if a named intervals set exists. Otherwise 'FALSE'.
}

\seealso{
\code{\link{gintervals.ls}}, \code{\link{gintervals.load}}, \code{\link{gintervals.rm}}, \code{\link{gintervals.save}}, \code{\link{gintervals}}, \code{\link{gintervals.2d}}
}
\examples{
gdb.init_examples()
gintervals.exists("annotations")
}
\keyword{ ~intervals }
