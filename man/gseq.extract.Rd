\name{gseq.extract}
\alias{gseq.extract}
\title{
Returns DNA sequences
}
\description{
Returns DNA sequences for given intervals
}
\usage{
gseq.extract(intervals)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{intervals}{intervals for which DNA sequence is returned}
}
\details{
This function returns an array of sequence strings for each interval from 'intervals'. If  intervals contain an additional 'strand' column and its value is '-1', the reverse-complementary sequence is returned.
}
\value{
An array of character strings representing DNA sequence.
}

\seealso{
\code{\link{gextract}}
}
\examples{
gdb.init_examples()
intervs <- gintervals(c(1, 2), 10000, 10020)
gseq.extract(intervs)
}
\keyword{ ~extract }
\keyword{ ~DNA }
\keyword{ ~sequence }
