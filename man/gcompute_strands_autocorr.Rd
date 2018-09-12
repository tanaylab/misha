\name{gcompute_strands_autocorr}
\alias{gcompute_strands_autocorr}
\title{
Computes auto-correlation between the strands for a file of mapped sequences
}
\description{
  Calculates auto-correlation between plus and minus strands for the given chromosome in a file of mapped sequences.
}
\usage{
gcompute_strands_autocorr(file, chrom, binsize, maxread = 400,
                          cols.order = c(9, 11, 13, 14),
                          min.coord = 0, max.coord = 3e+8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{the name of the file containing mapped sequences}
  \item{chrom}{chromosome for which the auto-correlation is computed}
  \item{binsize}{calculate the auto-correlation for bins in the range of [-maxread, maxread]}
  \item{maxread}{maximal length of the sequence used for statistics}
  \item{cols.order}{order of sequence, chromosome, coordinate and strand columns in file}
  \item{min.coord}{minimal coordinate used for statistics}
  \item{max.coord}{maximal coordinate used for statistics}
}
\details{
 This function calculates auto-correlation between plus and minus strands for the given chromosome in a file of mapped sequences. Each line in the file describes one read. Each column is separated by a TAB character.

The following columns must be presented in the file: sequence, chromosome, coordinate and strand. The position of these columns are controlled by 'cols.order' argument accordingly. The default value of 'cols.order' is a vector (9,11,13,14) meaning that sequence is expected to be found at column number 9, chromosome - at column 11, coordinate - at column 13 and strand - at column 14. The first column should be referenced by 1 and not by 0.

Coordinates that are not in [min.coord, max.coord] range are ignored.

gcompute_strands_autocorr outputs the total statistics and the auto-correlation given by bins. The size of the bin is indicated by 'binsize' parameter. Statistics is calculated for bins in the range of [-maxread, maxread].
}
\value{
Statistics for each strand and auto-correlation by given bins.
}

\examples{
gdb.init_examples()
gcompute_strands_autocorr(paste(GROOT, "reads", sep = "/"),
                          "chr1", 50, maxread = 300)
}
\keyword{ ~gcompute_strands_autocorr }
\keyword{ ~auto-correlation }
\keyword{ ~autocorrelation }
\keyword{ ~correlation }
