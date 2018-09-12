\name{gdir.cd}
\alias{gdir.cd}
\title{
Changes current working directory in Genomic Database
}
\description{
Changes current working directory in Genomic Database.
}
\usage{
gdir.cd(dir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dir}{directory path}
}
\details{
This function changes the current working directory in Genomic Database (not to be confused with shell's current working directory). The list of database objects - tracks, intervals, track variables - is rescanned recursively under 'dir'. Object names are updated with the respect to the new current working directory. Example: a track named 'subdir.dense' will be referred as 'dense' once current working directory is set to 'subdir'. All virtual tracks are removed.
}
\value{
None.
}

\seealso{
\code{\link{gdb.init}}, \code{\link{gdir.cwd}}, \code{\link{gdir.create}}, \code{\link{gdir.rm}}, \code{\link{gset_input_mode}}
}
\examples{
gdb.init_examples()
gdir.cd("subdir")
gtrack.ls()
gdir.cd("..")
gtrack.ls()
}
\keyword{ ~db }
\keyword{ ~data }
\keyword{ ~database }
\keyword{ ~cd }
\keyword{ ~dir }
\keyword{ ~directory }
\keyword{ ~folder }
