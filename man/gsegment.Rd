\name{gsegment}
\alias{gsegment}
\title{
Divides track expression into segments
}
\description{
Divides the values of track expression into segments by using Wilcoxon test.
}
\usage{
gsegment(expr, minsegment, maxpval = 0.05, onetailed = TRUE,
         intervals = ALLGENOME, iterator = NULL,
         intervals.set.out = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expr}{track expression}
  \item{minsegment}{minimal segment size}
  \item{maxpval}{maximal P-value that separates two adjacent segments}
  \item{onetailed}{if 'TRUE', Wilcoxon test is performed one tailed, otherwise two tailed}
  \item{intervals}{genomic scope for which the function is applied}
  \item{iterator}{track expression iterator of "fixed bin" type. If 'NULL' iterator is determined implicitly based on track expression.}
  \item{intervals.set.out}{intervals set name where the function result is optionally outputed}
}
\details{
This function divides the values of track expression into segments, where each segment size is at least of 'minsegment' size and the P-value of comparing the segment with the first 'minsegment' values from the next segment is at most 'maxpval'. Comparison is done using Wilcoxon (also known as Mann-Whitney) test.

If 'intervals.set.out' is not 'NULL' the result is saved as an intervals set. Use this parameter if the result size exceeds the limits of the physical memory.
}
\value{
If 'intervals.set.out' is 'NULL' a set of intervals where each interval represents a segment.
}

\seealso{
\code{\link{gscreen}}, \code{\link{gwilcox}}
}
\examples{
gdb.init_examples()
gsegment("dense_track", 5000, 0.0001)
}
\keyword{ ~segment }
\keyword{ ~wilcoxon }
\keyword{ ~Mann-Whitney }
