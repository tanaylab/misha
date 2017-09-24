\name{gintervals.rm}
\alias{gintervals.rm}
\title{
Deletes a named intervals set
}
\description{
Deletes a named intervals set.
}
\usage{
gintervals.rm(intervals.set, force = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{intervals.set}{name of an intervals set}
  \item{force}{if 'TRUE', supresses user confirmation of a named intervals set removal}
}
\details{
This function deletes a named intervals set from the Genomic Database. By default 'gintervals.rm' requires the user to interactively confirm the deletion. Set 'force' to 'TRUE' to suppress the user prompt.
}
\value{
None.
}

\seealso{
\code{\link{gintervals.save}}, \code{\link{gintervals.exists}}, \code{\link{gintervals.ls}}, \code{\link{gintervals}}, \code{\link{gintervals.2d}}
}
\examples{
gdb.init_examples()
intervs <- gintervals(c(1, 2))
gintervals.save("testintervs", intervs)
gintervals.ls()
gintervals.rm("testintervs", force = TRUE)
gintervals.ls()
}
\keyword{ ~intervals }
