\name{gintervals.ls}
\alias{gintervals.ls}
\title{
Returns a list of named intervals sets
}
\description{
Returns a list of named intervals sets in Genomic Database.
}
\usage{
gintervals.ls(pattern = "", ignore.case = FALSE, perl = FALSE,
              fixed = FALSE, useBytes = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pattern, ignore.case, perl, fixed, useBytes}{see 'grep'}
}
\details{
This function returns a list of named intervals sets that match the pattern (see 'grep'). If called without any arguments all named intervals sets are returned.
}
\value{
An array that contains the names of intervals sets.
}

\seealso{
\code{\link{grep}}, \code{\link{gintervals.exists}}, \code{\link{gintervals.load}}, \code{\link{gintervals.save}}, \code{\link{gintervals.rm}}, \code{\link{gintervals}}, \code{\link{gintervals.2d}}
}
\examples{
gdb.init_examples()
gintervals.ls()
gintervals.ls(pattern = "annot*")
}
\keyword{ ~intervals }
\keyword{ ~ls }
