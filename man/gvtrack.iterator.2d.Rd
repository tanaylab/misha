\name{gvtrack.iterator.2d}
\alias{gvtrack.iterator.2d}
\title{
Defines modification rules for a two-dimensional iterator in a virtual track
}
\description{
Defines modification rules for a two-dimensional iterator in a virtual track.
}
\usage{
gvtrack.iterator.2d(vtrack, sshift1 = 0, eshift1 = 0, sshift2 = 0,
                    eshift2 = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vtrack}{virtual track name}
  \item{sshift1}{shift of 'start1' coordinate}
  \item{eshift1}{shift of 'end1' coordinate}
  \item{sshift2}{shift of 'start2' coordinate}
  \item{eshift2}{shift of 'end2' coordinate}
}
\details{
This function defines modification rules for one-dimensional iterator intervals in a virtual track. 

Iterator interval's 'start1' coordinate is modified by adding 'sshift1'. Similarly 'end1', 'start2', 'end2' coordinates are altered by adding 'eshift1', 'sshift2' and 'eshift2' accordingly.
}
\value{
None.
}

\seealso{
\code{\link{gvtrack.create}}, \code{\link{gvtrack.iterator}}
}
\examples{
gdb.init_examples()
gvtrack.create("vtrack1", "rects_track")
gvtrack.iterator.2d("vtrack1", sshift1 = 1000, eshift1 = 2000)
gextract("rects_track", "vtrack1",
         gintervals.2d(1, 0, 5000, 2, 0, 5000))
}
\keyword{ ~virtual }
