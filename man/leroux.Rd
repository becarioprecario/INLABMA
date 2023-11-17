\name{leroux.inla}
\alias{leroux.inla}

\title{
Fit Leroux et al's spatial model.
}
\description{
This function fits the model by Leroux et al. for a given value
of the parameter \code{lambda}, i.e., the mixture parameter that
appears in the variance..
}
\usage{
leroux.inla(formula, d, W, lambda, improve = TRUE, fhyper = NULL, ...)
}

\arguments{
  \item{formula}{
Formula of the fixed effects.
}
  \item{d}{
A data.frame with the data to be used.
}
  \item{W}{
Adjacency matrix.
}
  \item{lambda}{
Parameter used in the mixture of the two precission matrices.
}
  \item{improve}{
Logical. Whether to improve the fitted models to obtain better
estimates of the marginal likelihoods.
}
  \item{fhyper}{
Extra arguments passed to the definition of the hyperparameters.
}
  \item{\dots}{
Extra arguments passed to function \code{inla}.
}
}
\details{
This function fits the model proposed by Leroux et al. (1999)
for a given value of parameter \code{lambda}. This parameter
controls the mixture between a diagonal precission (\code{lambda}=1) 
and an intrinsic CAR precission (\code{lambda}=0).

The marginal log-likelihood is corrected to add half the log-determinant of the precission matrix.
}
\value{
An INLA object.
}
\author{
Virgilio Gómez-Rubio <virgilio.gomez@uclm.es>
}
%%\note{
%%  ~~further notes~~
%%}

\seealso{
\code{\link{sem.inla}},\code{\link{slm.inla}},\code{\link{sdm.inla}}
}
%%\examples{
%%}

\references{
  Leroux B, Lei X, Breslow N (1999). Estimation of Disease Rates in Small 
  Areas: A New Mixed Model for Spatial Dependence. In M Halloran, D Berry 
  (eds.), Statistical Models in Epidemiology, the Environment and Clinical
  Trials, pp. 135-178. Springer-Verlag, New York.

  Roger S. Bivand, Virgilio Gómez-Rubio, Håvard Rue (2014). Approximate 
  Bayesian inference for spatial econometrics models. Spatial Statistics, 
  Volume 9, 146-165.

  Roger S. Bivand, Virgilio Gómez-Rubio, Håvard Rue (2015). Spatial
  Data Analysis with R-INLA with Some Extensions. Journal of
  Statistical Software, 63(20), 1-31. URL http://www.jstatsoft.org/v63/i20/.
}
\keyword{models}
