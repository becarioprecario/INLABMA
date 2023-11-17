\name{INLABMA}
\alias{INLABMA}
\title{
Perform complete Bayesian Model Averaging on some Spatial Econometrics models
}
\description{
This function performs Bayesian Model Averaging on a list of 
different Spatial Econometrics models. These models have been computed
under different values of the spatial autocorrelation parameter \code{rho}.
}
\usage{
INLABMA(models, rho, logrhoprior = rep(1, length(rho)), impacts = FALSE, 
 usenormal = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{models}{
List of INLA models, computed for different values of \code{rho}.
}
  \item{rho}{
A vector with the values of \code{rho} used to compute \code{models}.
}
  \item{logrhoprior}{
Vector with the values of the log-prior density of \code{rho}.
}
  \item{impacts}{
Logical. Whether impacts should be computed.
}
  \item{usenormal}{
Logical. Whether the posterior marginal of \code{rho} is assumed to be
Gaussian.
}
}
\details{
This functions perfomrs BMA on most of the compponents of an INLA model
using the marginal likelihoods of the models and the provided
log-prior density of \code{rho}. 
}
\value{
A list with the averaged components. Another component called
\code{rho} is added, with its posterior marginal and some other
summary information.
}
%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Virgilio Gómez-Rubio <virgilio.gomez@uclm.es>
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{sem.inla}}, \code{\link{slm.inla}}, \code{\link{sdm.inla}} 
}
%%\examples{
%%}

\references{
  Roger S. Bivand, Virgilio Gómez-Rubio, Håvard Rue (2014). Approximate 
  Bayesian inference for spatial econometrics models. Spatial Statistics, 
  Volume 9, 146-165.

  Roger S. Bivand, Virgilio Gómez-Rubio, Håvard Rue (2015). Spatial
  Data Analysis with R-INLA with Some Extensions. Journal of
  Statistical Software, 63(20), 1-31. URL http://www.jstatsoft.org/v63/i20/.
}

\keyword{smooth}
