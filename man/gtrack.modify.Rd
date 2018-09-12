\name{gtrack.modify}
\alias{gtrack.modify}
\title{
Modifies track contents
}
\description{
Modifies 'Dense' track contents.
}
\usage{
gtrack.modify(track, expr, intervals = ALLGENOME)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{track}{track name}
  \item{expr}{track expression}
  \item{intervals}{genomic scope for which track is modified}
}
\details{
This function modifies the contents of a 'Dense' track by the values of 'expr'. 'intervals' argument controls which portion of the track is modified. The iterator policy is set internally to the bin size of the track.
}
\value{
None.
}

\seealso{
\code{\link{gtrack.create}}, \code{\link{gtrack.rm}}
}
\examples{
gdb.init_examples()
intervs <- gintervals(1, 300, 800)
gextract("dense_track", intervs)
gtrack.modify("dense_track", "dense_track * 2", intervs)
gextract("dense_track", intervs)
gtrack.modify("dense_track", "dense_track / 2", intervs)
}
\keyword{ ~modify }
\keyword{ ~track }
