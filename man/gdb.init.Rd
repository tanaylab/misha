\name{gdb.init}
\alias{gdb.init}
\alias{gdb.init.examples}
\title{
Initializes connection with Genomic Database
}
\description{
Initializes connection with Genomic Database: loads the list of tracks, intervals, etc.
}
\usage{
gdb.init(groot, dir = NULL, rescan = FALSE)
gdb.init_examples()
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{groot}{the root directory of the Genomic Database}
  \item{dir}{the current working directory inside the Genomic Database}
  \item{rescan}{indicates whether the file structure should be rescanned}
}
\details{
'gdb.init' initializes the connection with the Genomic Database. It is typically called first prior to any other function. When the package is attached it internally calls to 'gdb.init.examples' which opens the connection with the database located at 'PKGDIR/trackdb/test' directory, where 'PKGDIR' is the directory where the package is installed.

The current working directory inside the Genomic Database is set to 'dir'. If 'dir' is 'NULL', the current working directory is set to 'GROOT/tracks'.

If 'rescan' is 'TRUE', the list of tracks and intervals is achieved by rescanning directory structure under the current current working directory. Otherwise 'gdb.init' attempts to use the cached list that resides in 'groot/.db.cache' file.

Upon completion the connection is established with the database. If auto-completion mode is switched on (see 'gset_input_method') the list of tracks and intervals sets is loaded and added as variables to the global environment allowing auto-completion of object names with <TAB> key. Also a few global variables are defined. These variables should not be modified by user.

\tabular{ll}{
    GROOT \tab Root directory of Genomic Database\cr
    GWD \tab Current working directory inside Genomic Database\cr
    GTRACKS \tab List of all available tracks\cr
    GINTERVS \tab List of all available intervals\cr
    GVTRACKS \tab List of all available virtual tracks\cr
    ALLGENOME \tab List of all chromosomes and their sizes\cr
    GITERATOR.INTERVALS \tab A set of iterator intervals for which the track expression is evaluated\cr
}
}
\value{
None.
}

\seealso{
\code{\link{gdb.reload}}, \code{\link{gdb.create}}, \code{\link{gdir.cd}}, \code{\link{gtrack.ls}}, \code{\link{gintervals.ls}}, \code{\link{gvtrack.ls}}, \code{\link{gset_input_mode}}
}
\keyword{ ~db }
\keyword{ ~data }
\keyword{ ~database }
