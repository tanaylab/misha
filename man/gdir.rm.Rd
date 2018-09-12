\name{gdir.rm}
\alias{gdir.rm}
\title{
Deletes a directory from Genomic Database
}
\description{
Deletes a directory from Genomic Database.
}
\usage{
gdir.rm(dir, recursive = FALSE, force = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dir}{directory path}
  \item{recursive}{if 'TRUE', the directory is deleted recursively}
  \item{force}{if 'TRUE', supresses user confirmation of tracks/intervals removal}
}
\details{
This function deletes a directory from Genomic Database. If 'recursive' is 'TRUE', the directory is deleted with all the files/directories it contains. If the directory contains tracks or intervals, the user is prompted to confirm the deletion. Set 'force' to 'TRUE' to suppress the prompt.
}
\value{
None.
}

\seealso{
\code{\link{gdb.init}}, \code{\link{gdir.create}}, \code{\link{gdir.cd}}, \code{\link{gdir.cwd}}
}
\keyword{ ~db }
\keyword{ ~data }
\keyword{ ~database }
\keyword{ ~dir }
\keyword{ ~directory }
\keyword{ ~folder }
\keyword{ ~rm }
