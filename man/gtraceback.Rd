\name{gtraceback}
\alias{gtraceback}
\title{
Prints call stack of the last uncaught error
}
\description{
Prints call stack of the last uncaught error in a friendly way.
}
\usage{
gtraceback(x = NULL, max.lines = getOption("deparse.max.lines"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{see 'traceback'}
  \item{max.lines}{see 'traceback'}
}
\details{
Similarly to 'traceback' this function prints the call stack of the last uncaught error. Yet 'gtraceback' does it in a more friendly way by omitting the calls that occurred inside the library.
}
\value{
See 'traceback'.
}

\seealso{
\code{\link{traceback}}
}
\examples{
gdb.init_examples()
f <- function() { gscreen("blablabla") }
f()
gtraceback()
}
\keyword{ ~trace }
\keyword{ ~error }
\keyword{ ~exception }
