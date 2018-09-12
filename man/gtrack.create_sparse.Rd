\name{gtrack.create_sparse}
\alias{gtrack.create_sparse}
\title{
Creates a 'Sparse' track from intervals and values
}
\description{
Creates a 'Sparse' track from intervals and values.
}
\usage{
gtrack.create_sparse(track, description, intervals, values)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{track}{track name}
  \item{description}{a character string description}
  \item{intervals}{a set of one-dimensional intervals}
  \item{values}{an array of numeric values - one for each interval}
}
\details{
This function creates a new 'Sparse' track with values at given intervals. 'description' is added as a track attribute.
}
\value{
None.
}

\seealso{
\code{\link{gtrack.create}}, \code{\link{gtrack.2d.create}}, \code{\link{gtrack.smooth}}, \code{\link{gtrack.modify}}, \code{\link{gtrack.rm}}, \code{\link{gtrack.info}}, \code{\link{gdir.create}}
}
\examples{
gdb.init_examples()
intervs <- gintervals.load("annotations")
gtrack.create_sparse("test_sparse", "Test track", intervs,
                     1 : dim(intervs)[1])
gextract("test_sparse", ALLGENOME)
gtrack.rm("test_sparse", force = TRUE)
}
\keyword{ ~create }
\keyword{ ~sparse }
\keyword{ ~track }
