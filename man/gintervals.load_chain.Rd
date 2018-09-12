\name{gintervals.load_chain}
\alias{gintervals.load_chain}
\title{
Loads assembly conversion table from a chain file
}
\description{
Loads assembly conversion table from a chain file.
}
\usage{
gintervals.load_chain(file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{name of chain file}
}
\details{
This function reads a file in 'chain' format and returns assembly conversion table that can be used in 'gtrack.liftover' and 'gintervals.liftover'.

Note: chain file might map a few different source intervals into a single target one. These ambigous mappings are not presented in the data frame returned by 'gintervals.load_chain'.
}
\value{
A data frame representing assembly conversion table.
}

\seealso{
\code{\link{gintervals.liftover}}, \code{\link{gtrack.liftover}}
}
\examples{
gdb.init_examples()
chainfile <- paste(GROOT, "data/test.chain", sep = "/")
gintervals.load_chain(chainfile)
}
\keyword{ ~intervals }
\keyword{ ~liftover }
\keyword{ ~chain }
