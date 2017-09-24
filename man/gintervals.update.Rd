\name{gintervals.update}
\alias{gintervals.update}
\title{
Updates a named intervals set
}
\description{
Updates a named intervals set.
}
\usage{
gintervals.update(intervals.set, intervals, chrom = NULL,
                  chrom1 = NULL, chrom2 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{intervals.set}{name of an intervals set}
  \item{intervals}{intervals or 'NULL'}
  \item{chrom}{chromosome for 1D intervals set}
  \item{chrom1}{first chromosome for 2D intervals set}
  \item{chrom2}{second chromosome for 2D intervals set}
}
\details{
This function replaces all intervals of given chromosome (or chromosome pair) within 'intervals.set' with 'intervals'. Chromosome is specified by 'chrom' for 1D intervals set or 'chrom1', 'chrom2' for 2D intervals set.

If 'intervals' is 'NULL' all intervals of given chromosome are removed from 'intervals.set'.
}
\value{
None.
}

\seealso{
\code{\link{gintervals.save}}, \code{\link{gintervals.load}}, \code{\link{gintervals.exists}}, \code{\link{gintervals.ls}}
}
\examples{
gdb.init_examples()
intervs <- gscreen("sparse_track > 0.2",
                   gintervals(c(1,2), 0, 10000))
gintervals.save("testintervs", intervs)
gintervals.load("testintervs")
gintervals.update("testintervs", intervs[intervs$chrom == "chr2",][1:5,], chrom = 2)
gintervals.load("testintervs")
gintervals.update("testintervs", NULL, chrom = 2)
gintervals.load("testintervs")
gintervals.rm("testintervs", force = TRUE)
}
\keyword{ ~intervals }

