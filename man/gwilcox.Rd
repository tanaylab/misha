\name{gwilcox}
\alias{gwilcox}
\title{
Calculates Wilcoxon test on sliding windows over track expression
}
\description{
Calculates Wilcoxon test on sliding windows over the values of track expression.
}
\usage{
gwilcox(expr, winsize1, winsize2, maxpval = 0.05,
        onetailed = TRUE, what2find = 1, intervals = ALLGENOME,
        iterator = NULL, intervals.set.out = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expr}{track expression}
  \item{winsize1}{number of values in the first sliding window}
  \item{winsize2}{number of values in the second sliding window}
  \item{maxpval}{maximal P-value}
  \item{onetailed}{if 'TRUE', Wilcoxon test is performed one tailed, otherwise two tailed}
  \item{what2find}{if '-1', lows are searched. If '1', peaks are searched. If '0', both peaks and lows are searched}
  \item{intervals}{genomic scope for which the function is applied}
  \item{iterator}{track expression iterator of "fixed bin" type. If 'NULL' iterator is determined implicitly based on track expression.}
  \item{intervals.set.out}{intervals set name where the function result is optionally outputed}
}
\details{
This function runs a Wilcoxon test (also known as a Mann-Whitney test) over the values of track expression in the two sliding windows having an identical center. The sizes of the windows are specified by 'winsize1' and 'winsize2'. 'gwilcox' returns intervals where the smaller window tested against a larger window gives a P-value below 'maxpval'. The test can be one or two tailed.

'what2find' argument controls what should be searched: peaks, lows or both.

If 'intervals.set.out' is not 'NULL' the result is saved as an intervals set. Use this parameter if the result size exceeds the limits of the physical memory.
}
\value{
If 'intervals.set.out' is 'NULL' a data frame representing the intervals with an additional 'pval' column where P-value is below 'maxpval'.
}

\seealso{
\code{\link{gscreen}}, \code{\link{gsegment}}
}
\examples{
gdb.init_examples()
gwilcox("dense_track", 100000, 1000, maxpval = 0.01,
        what2find = 1)
}
\keyword{ ~wilcoxon }
\keyword{ ~Mann-Whitney }
