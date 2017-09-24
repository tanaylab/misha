\name{gintervals.import_genes}
\alias{gintervals.import_genes}
\title{
Imports genes and annotations from files
}
\description{
Imports genes and annotations from files.
}
\usage{
gintervals.import_genes(genes.file, annots.file = NULL, annots.names = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genes.file}{name or URL of file that contains genes}
  \item{annots.file}{name of URL file that contains annotations. If 'NULL' no annotations are imported}
  \item{annots.names}{annotations names}
}
\details{
This function reads a definition of genes from 'genes.file' and returns four sets of intervals: TSS, exons, 3utr and 5utr. In addition to the regular intervals columns 'strand' column is added. It contains '1' values for '+' strands and '-1' values for '-' strands. 

If annotation file 'annots.file' is given then annotations are attached too to the intervals. The names of the annotations as they would appear in the return value must be specified in 'annots.names' argument.

Both 'genes.file' and 'annots.file' can be either a file path or URL in a form of 'ftp://[address]/[file]'. Files that these arguments point to can be zipped or unzipped.

Examples of 'genes.file' and 'annots.file' can be found here:

\code{ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz}
\code{ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/kgXref.txt.gz}

If a few intervals overlap (for example: two TSS regions) they are all unified to an interval that covers the whole overlapping region. 'strand' value is set to '0' if two or more of the overlapping intervals have different strands. The annotations of the overlapping intervals are concatenated to a single character string separated by semicolons. Identical values of overlapping intervals' annotation are eliminated.
}
\value{
A list of four intervals sets named 'tss', 'exons', 'utr3' and 'utr5'. 'strand' column and annotations are attached to the intevals.
}

\seealso{
\code{\link{gintervals}}, \code{\link{gdb.create}}
}
\keyword{ ~intervals }
\keyword{ ~import }
\keyword{ ~genes }
