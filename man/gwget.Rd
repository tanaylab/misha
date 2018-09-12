\name{gwget}
\alias{gwget}
\title{
Downloads files from FTP server
}
\description{
Downloads multiple files from FTP server
}
\usage{
gwget(url, path = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{url}{URL of FTP server}
  \item{path}{directory path where the downloaded files are stored}
}
\details{
This function downloads files from FTP server given by 'url'. The address in 'url' can contain wildcards to download more than one file at once. Files are downloaded to a directory given by 'path' argument.  If 'path' is 'NULL', file are downloaded into 'GROOT/downloads'.
}
\value{
An array of file names that have been downloaded.
}

\seealso{
\code{\link{gtrack.import_set}}
}
\examples{
gdb.init_examples()
gwget("ftp://hgdownload.cse.ucsc.edu/logs/2012/access_log.2012071*")
}
\keyword{ ~ftp }
