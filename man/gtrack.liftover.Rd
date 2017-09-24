\name{gtrack.liftover}
\alias{gtrack.liftover}
\title{
Imports a track from another assembly
}
\description{
Imports a track from another assembly.
}
\usage{
gtrack.liftover(track, description, src.track.dir, chain)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{track}{name of a created track}
  \item{description}{a character string description}
  \item{src.track.dir}{path to the directory of the source track}
  \item{chain}{name of chain file or data frame as returned by 'gintervals.load_chain'}
}
\details{
This function imports a track located in 'src.track.dir' of another assembly to the current database. Chain file instructs how the conversion of coordinates should be done. It can be either a name of a chain file or a data frame in the same format as returned by 'gintervals.load_chain' function. The name of the newly created track is specified by 'track' argument and 'description' is added as a track attribute.
}
\value{
None.
}

\seealso{
\code{\link{gintervals.load_chain}}, \code{\link{gintervals.liftover}}
}
\keyword{ ~track }
\keyword{ ~liftover }
\keyword{ ~chain }
