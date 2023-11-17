\name{BMArho}
\alias{BMArho}
\title{
Compute BMA of \code{fitted.values} from a list of INLA objects
}
\description{
This functions performs a weighted average of the component
\code{fitted.values} from a list of INLA objects.
}
\usage{
BMArho(models, rho, logrhoprior = rep(1, length(rho)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{models}{
List of INLA models fitted for different values of \code{rho}
}
  \item{rho}{
Vector fo values of \code{rho} used to compute the list in \code{models}.
}
  \item{logrhoprior}{
Log-prior density for each value of \code{rho}.
}
}
\details{
The different \code{fitted.values} are weighted using the values of the
marginal likelihood of the fitted models and the prior of parameter
\code{rho}. \code{rho} is a parameter that is fixed when computing
\code{models} and that have a log-prior density  defined in \code{pogrhoprior}.
}
\value{
Vector of averaged fitted values.
}
%5\references{
%% ~put references to the literature/web site here ~
%%}
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

\references{
  Roger S. Bivand, Virgilio Gómez-Rubio, Håvard Rue (2014). Approximate 
  Bayesian inference for spatial econometrics models. Spatial Statistics, 
  Volume 9, 146-165.

  Roger S. Bivand, Virgilio Gómez-Rubio, Håvard Rue (2015). Spatial
  Data Analysis with R-INLA with Some Extensions. Journal of
  Statistical Software, 63(20), 1-31. URL http://www.jstatsoft.org/v63/i20/.
}


\keyword{smooth}
