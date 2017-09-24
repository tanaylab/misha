\name{gtrack.create}
\alias{gtrack.create}
\title{
Creates a track from a track expression
}
\description{
Creates a track from a track expression.
}
\usage{
gtrack.create(track, description, expr, iterator = NULL,
              band = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{track}{track name}
  \item{description}{a character string description}
  \item{expr}{track expression}
  \item{iterator}{track expression iterator. If 'NULL' iterator is determined implicitly based on track expression.}
  \item{band}{track expression band. If 'NULL' no band is used.}
}
\details{
This function creates a new track named track. The values of the track are determined by evaluation of 'expr' - a numeric track expression. The type of the new track is determined by the type of the iterator. 'Fixed bin', 'Sparse' or 'Rectangles' track can be created accordingly. 'description' is added as a track attribute.
}
\value{
None.
}

\seealso{
\code{\link{gtrack.2d.create}}, \code{\link{gtrack.create_sparse}}, \code{\link{gtrack.smooth}}, \code{\link{gtrack.modify}}, \code{\link{gtrack.rm}}, \code{\link{gtrack.info}}, \code{\link{gdir.create}}
}
\examples{
gdb.init_examples()

## Creates a new track that is a sum of values from 'dense' and
## 2 * non-nan values of 'sparse' track. The new track type is
## Dense with a bin size that equals to '70'.
gtrack.create("mixed_track", "Test track",
              "dense_track +
              replace(sparse_track, is.nan(sparse_track), 0) * 2",
              iterator = 70)
gtrack.info("mixed_track")
gtrack.rm("mixed_track", force = T)
}
\keyword{ ~create }
\keyword{ ~track }
