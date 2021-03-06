\name{gset_input_mode}
\alias{gset_input_mode}
\title{
Sets input mode and auto-completion for track expressions, track names, virtual tracks and interval sets
}
\description{
Sets input mode and auto-completion for track expressions, track names, virtual tracks and interval sets.
}
\usage{
gset_input_mode(autocompletion = FALSE, interactive = FALSE)
}
\details{
This function enables / disables auto-completion for track names, virtual tracks and interval sets. It also controls whether these objects together with track expressions should be passed as strings or "as is" to the various package functions.

If 'autocompletion' is 'TRUE' all the track names, virtual tracks and intervals sets are defined as R variables (auxiliary variables) which allows them to be auto-completed by TAB key. The values of these variables are meaningless for the user and they should not be altered.

If 'interactive' is 'TRUE' track names, virtual tracks, interval sets and track expressions are passed to the package functions "as is", i.e. unquoted.

'autocompletion' is required to be switched on if 'interactive' mode is on too.

Please beware of the consequences of using interactive mode as it creates a bunch of new variables in R environment. Though collision with the existing variables is checked at the time of the call to 'gset_input_mode', yet nothing prevents the user to modify the value of the auxiliary variables later. This might cause unexpected behaviour in some of the package functions. Also the auxiliary variables are automatically undefined once the interactive mode is switched off. User who mistakenly uses auxiliary variables to store the data might therefore accidentially loose it. 
}
\value{
None.
}
\examples{
gdb.init_examples()
gset_input_mode(interactive = FALSE)
gsummary("dense_track + 10")
gset_input_mode(autocompletion = TRUE, interactive = TRUE)
gsummary(dense_track + 10)

# roll-back to default input mode
gset_input_mode()
}
\keyword{ ~autocompletion }
\keyword{ ~interactive }
\keyword{ ~quoted }

