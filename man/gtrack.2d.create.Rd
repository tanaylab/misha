\name{gtrack.2d.create}
\alias{gtrack.2d.create}
\title{
Creates a 'Rectangles' track from intervals and values
}
\description{
Creates a 'Rectangles' track from intervals and values.
}
\usage{
gtrack.2d.create(track, description, intervals, values)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{track}{track name}
  \item{description}{a character string description}
  \item{intervals}{a set of two-dimensional intervals}
  \item{values}{an array of numeric values - one for each interval}
}
\details{
This function creates a new 'Rectangles' (two-dimensional) track with values at given intervals. 'description' is added as a track attribute.
}
\value{
None.
}

\seealso{
\code{\link{gtrack.create}}, \code{\link{gtrack.create_sparse}}, \code{\link{gtrack.smooth}}, \code{\link{gtrack.modify}}, \code{\link{gtrack.rm}}, \code{\link{gtrack.info}}, \code{\link{gdir.create}}, \code{\link{gtrack.attr.get}}
}
\examples{
gdb.init_examples()
intervs1 <- gintervals.2d(1, (1 : 5) * 200, (1 : 5) * 200 + 100,
                          1, (1 : 5) * 300, (1 : 5) * 300 + 200)
intervs2 <- gintervals.2d("X", (7 : 10) * 100, (7 : 10) * 100 + 50,
                          2, (1 : 4) * 200, (1 : 4) * 200 + 130)
intervs <- rbind(intervs1, intervs2)
gtrack.2d.create("test_rects", "Test 2d track", intervs,
                 runif(dim(intervs)[1], 1, 100))
gextract("test_rects", ALLGENOME)
gtrack.rm("test_rects", force = TRUE)
}
\keyword{ ~create }
\keyword{ ~track }
