\name{gtrack.ls}
\alias{gtrack.ls}
\title{
Returns a list of track names
}
\description{
Returns a list of track names in Genomic Database.
}
\usage{
gtrack.ls(..., ignore.case = FALSE, perl = FALSE, fixed = FALSE,
          useBytes = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{these arguments are of either form 'pattern' or 'attribute = pattern'}
  \item{ignore.case, perl, fixed, useBytes}{see 'grep'}
}
\details{
This function returns a list of tracks whose name or track attribute value match a  pattern (see 'grep'). If called without any arguments all tracks are returned.

If pattern is specified without a track attribute (i.e. in the form of 'pattern') then filtering is applied to the track names. If pattern is supplied with a track attribute (i.e. in the form of 'name = pattern') then track attribute is matched against the pattern.

Multiple patterns are applied one after another. The resulted list of tracks should match all the patterns.
}
\value{
An array that contains the names of tracks that match the supplied patterns.
}

\seealso{
\code{\link{grep}}, \code{\link{gtrack.exists}}, \code{\link{gtrack.create}}, \code{\link{gtrack.rm}}
}
\examples{
gdb.init_examples()

# get all track names
gtrack.ls()

# get track names that match the pattern "den*"
gtrack.ls("den*")

# get track names whose "created.by" attribute match the pattern
# "create_sparse"
gtrack.ls(created.by = "create_sparse")

# get track names whose names match the pattern "den*" and whose
# "created.by" attribute match the pattern "track"
gtrack.ls("den*", created.by = "track")
}
\keyword{ ~intervals }
\keyword{ ~ls }
