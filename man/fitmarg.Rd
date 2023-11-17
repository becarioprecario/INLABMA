\name{fitmarg}
\alias{fitmarg}
\title{
Fit posterior marginal distributions to points
}
\description{
Compute (and re-scale, if necessary) the marginal from a set of
points \code{x} and values of log-likelihood \code{logy} and 
log-prior density \code{logp}.
}
\usage{
fitmarg(x, logy, logp = 0, usenormal = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Values of the random variable.
}
  \item{logy}{
Log-likelihood.
}
  \item{logp}{
Log-prior density.
}
  \item{usenormal}{
Whether use a Normal distribution for the fitted marginal.
}
}
\details{
Fits a marginal at a set of points \code{x} from their log-likelihood
and log-prior. The fitted marginal is re-scaled to integrate one if
necessary. If \code{usenormal=TRUE} then the fitted marginal is supposed
to be Normal, which is computed using the posterior mean and standard
deviation of \code{x}.
}
\value{
A function with the fitted marginal is returned.
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Virgilio Gómez-Rubio <virgilio.gomez@uclm.es>
}
%%\note{
%%  ~~further notes~~
%%}

\seealso{
\code{\link{fitmargBMA}}, \code{\link{fitmargBMA2}},\code{\link{mysplinefun}}
}
%%\examples{
%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{smooth}
