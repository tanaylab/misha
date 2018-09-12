\name{gdb.reload}
\alias{gdb.reload}
\title{
Reloads database from the disk
}
\description{
  Reloads database from disk: list of tracks, intervals, etc.
}
\usage{
gdb.reload(rescan = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rescan}{indicates whether the file structure should be rescanned}
}
\details{
Reloads Genomic Database from disk: list of tracks, intervals, etc. Use this function if you manually add tracks or if for any reason the database becomes corrupted. If 'rescan' is 'TRUE', the list of tracks and intervals is achieved by rescanning directory structure under the current current working directory. Otherwise 'gdb.reload' attempts to use the cached list that resides in 'GROOT/.db.cache' file.
}

\seealso{
\code{\link{gdb.init}}, \code{\link{gdb.create}}, \code{\link{gdir.cd}}, \code{\link{gset_input_mode}}
}
\keyword{ ~db }
