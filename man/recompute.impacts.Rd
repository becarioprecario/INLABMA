\name{recompute.impacts}
\alias{recompute.impacts}

\title{
Recompute the impact summaries from the marginals
}
\description{
This functions recomputes the impacts summaries using the (approximated)
marginals rather than by weighting on the different summaries.
}
\usage{
recompute.impacts(obj, impacts = c("total", "direct", "indirect"))
}

\arguments{
  \item{obj}{
Object with a resulting model obtained by Bayesian Model Averaging with INLABMA.
}
  \item{impacts}{
Types of impacts to recompute.
}
}
\details{
This function uses the impacts marginals to compute some summary statistics.
By default, the summary of the impacts is obtained by weighting the different
summaries used in Bayesian MOdel Averaging with function \code{INLABMA}.
}
\value{
Original object with the updated summary statistics of the impacts.
}
\references{
Bivand et al. (2013)
}
\author{
Virgilio Gómez-Rubio <virgilio.gomez@uclm.es>
}
%%\note{
%%  ~~further notes~~
%%}


\seealso{
\code{\link{INLABMA}}
}
%%\examples{
%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
