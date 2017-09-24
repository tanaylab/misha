\name{gintervals.load}
\alias{gintervals.load}
\title{
Loads a named intervals set
}
\description{
Loads a named intervals set.
}
\usage{
gintervals.load(intervals.set, chrom = NULL, chrom1 = NULL,
                chrom2 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{intervals.set}{name of an intervals set}
  \item{chrom}{chromosome for 1D intervals set}
  \item{chrom1}{first chromosome for 2D intervals set}
  \item{chrom2}{second chromosome for 2D intervals set}
}
\details{
This function loads and returns intervals stored in a named intervals set.

If intervals set contains 1D intervals and 'chrom' is not 'NULL' only the intervals of the given chromosome are returned.

Likewise if intervals set contains 2D intervals and 'chrom1', 'chrom2' are not 'NULL' only the intervals of the given pair of chromosomes are returned.

For big intervals sets 'chrom' parameter (1D case) / 'chrom1', 'chrom2' parameters (2D case) must be specified. In other words: big intervals sets can be loaded only by chromosome or chromosome pair.
}
\value{
A data frame representing the intervals.
}

\seealso{
\code{\link{gintervals.save}}, \code{\link{gintervals.is.bigset}}, \code{\link{gintervals.exists}}, \code{\link{gintervals.ls}}, \code{\link{gintervals}}, \code{\link{gintervals.2d}}
}
\examples{
gdb.init_examples()
gintervals.load("annotations")
}
\keyword{ ~intervals }
