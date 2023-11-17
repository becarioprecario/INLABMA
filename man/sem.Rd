\name{sem.inla}
\alias{sem.inla}
\alias{slm.inla}
\alias{sdm.inla}
\alias{sac.inla}
\title{
Fit spatial econometrics models with INLA
}
\description{
These functions fit some spatial econometrics models for a given
value of \code{rho} (the spatial autocorrelation parameter).
\code{sem.inla} fits a spatial error model, \code{slm} fits a spatial lag model
and \code{sdm.inla} fits a spatial Durbin model.
}
\usage{

sem.inla(formula, d, W, rho, improve = TRUE, impacts = FALSE, fhyper = NULL, 
   probit = FALSE, ...)
slm.inla(formula, d, W, rho, mmatrix = NULL, improve = TRUE, impacts = FALSE, 
   fhyper = NULL, probit = FALSE, ...)
sdm.inla(formula, d, W, rho, mmatrix = NULL, intercept = TRUE, impacts = FALSE,
   improve = TRUE, fhyper = NULL, probit = FALSE, ...)
sac.inla(formula, d, W.rho, W.lambda, rho, lambda, mmatrix = NULL, 
   improve = TRUE, impacts = FALSE, fhyper = NULL, probit = FALSE, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
Formula with the response variable, the fixed effects and, possibly, other
non-linear effects.
}
  \item{d}{
Data.frame with the data.
}
  \item{W}{
Adjacency matrix.
}
  \item{rho}{
Value of the spatial autocorrelation parameter. For the SAC model, spatial
autocorrelation term on the response.
}
  \item{W.rho}{
For the SAC model, adjacency matrix associated to the autocorrelation on the
response.
}
  \item{W.lambda}{
For the SAC model, adjacency matrix associated to the autocorrelation on the
error term.
}
\item{lambda}{
For the SAC model, spatial autocorrelation of the error term.
}

  \item{mmatrix}{
Design matrix of fixed effects.
}
  \item{intercept}{
Logical. Whether an intercept has been included in the model.
}
  \item{improve}{
Logical. Whether improve model fitting (this may require more computing time).
}
  \item{impacts}{
Logical. Whether impacts are computed.
}
  \item{fhyper}{
Options to be passed to the definition of the hyper-parameters in the
spatial effects.
}
  \item{probit}{
Logical. Whether a probit model is used. Note this is only used
when computing the impacts and that argument \code{family} must
be set accordingly.
}
  \item{\dots}{
Other arguments passed to function \code{inla}.
}
}
\details{
These functions fit a spatial econometrics model with a fixed value
of the spatial autocorrelation parameter \code{rho}. 

In addition, the marginal -log-likelihood is corrected to account for
the variance-covariance matrix of the error term or random effects.

}
\value{
An \code{inla} object.
}
\references{
  Roger S. Bivand, Virgilio Gómez-Rubio, Håvard Rue (2014). Approximate 
  Bayesian inference for spatial econometrics models. Spatial Statistics, 
  Volume 9, 146-165.

  Roger S. Bivand, Virgilio Gómez-Rubio, Håvard Rue (2015). Spatial
  Data Analysis with R-INLA with Some Extensions. Journal of
  Statistical Software, 63(20), 1-31. URL http://www.jstatsoft.org/v63/i20/.

  Virgilio Gómez-Rubio and Francisco-Palmí Perales (2016).  Spatial Models with the Integrated Nested Laplace Approximation within Markov Chain Monte Carlo. Submitted.

}
\author{
Virgilio Gómez-Rubio <virgilio.gomez@uclm.es>
}
%%\note{
%%  ~~further notes~~
%%}

\seealso{
\code{\link{leroux.inla}}
}

%Example taken from the spdep::columbus manual page
\examples{

\dontrun{

if(requireNamespace("INLA", quietly = TRUE)) {
require(INLA)
require(spdep)

data(columbus)

lw <- nb2listw(col.gal.nb, style="W")

#Maximum Likelihood (ML) estimation
colsemml <- errorsarlm(CRIME ~ INC + HOVAL, data=columbus, lw, method="eigen", 
	quiet=FALSE)
colslmml <- lagsarlm(CRIME ~ INC + HOVAL, data=columbus, lw, method="eigen", 
	type="lag", quiet=FALSE)
colsdmml <- lagsarlm(CRIME ~ INC + HOVAL, data=columbus, lw, method="eigen", 
	type="mixed", quiet=FALSE)

#Define grid on rho
rrho<-seq(-1, .95, length.out=40)

#Adjacency matrix
W <- as(as_dgRMatrix_listw(nb2listw(col.gal.nb)), "CsparseMatrix")
#Index for spatial random effects
columbus$idx<-1:nrow(columbus)

#Formula
form<- CRIME ~ INC + HOVAL

zero.variance = list(prec=list(initial = 25, fixed=TRUE))



seminla<-mclapply(rrho, function(rho){

                sem.inla(form, d=columbus, W=W, rho=rho,
                        family = "gaussian", impacts=FALSE,
                        control.family = list(hyper = zero.variance),
                        control.predictor=list(compute=TRUE),
                        control.compute=list(dic=TRUE, cpo=TRUE),
                        control.inla=list(print.joint.hyper=TRUE), 
				#tolerance=1e-20, h=1e-6),
			verbose=FALSE
                )

})



slminla<-mclapply(rrho, function(rho){

                slm.inla(form, d=columbus, W=W, rho=rho,
                        family = "gaussian", impacts=FALSE,
                        control.family = list(hyper = zero.variance),
                        control.predictor=list(compute=TRUE),
                        control.compute=list(dic=TRUE, cpo=TRUE),
                        control.inla=list(print.joint.hyper=TRUE), 
				#tolerance=1e-20, h=1e-6),
			verbose=FALSE
                )
})


sdminla<-mclapply(rrho, function(rho){

                sdm.inla(form, d=columbus, W=W, rho=rho,
                        family = "gaussian", impacts=FALSE,
                        control.family = list(hyper = zero.variance),
                        control.predictor=list(compute=TRUE),
                        control.compute=list(dic=TRUE, cpo=TRUE),
                        control.inla=list(print.joint.hyper=TRUE), 
				#tolerance=1e-20, h=1e-6),
			verbose=FALSE
                )
})

#BMA using a uniform prior (in the log-scale) and using a Gaussian 
#approximation to the marginal
sembma<-INLABMA(seminla, rrho, 0, usenormal=TRUE)
slmbma<-INLABMA(slminla, rrho, 0, usenormal=TRUE)
sdmbma<-INLABMA(sdminla, rrho, 0, usenormal=TRUE)

#Display results
plot(sembma$rho$marginal, type="l", ylim=c(0,5))
lines(slmbma$rho$marginal, lty=2)
lines(sdmbma$rho$marginal, lty=3)
#Add ML estimates
abline(v=colsemml$lambda, col="red")
abline(v=colslmml$rho, col="red", lty=2)
abline(v=colsdmml$rho, col="red", lty=3)
#Legend
legend(-1,5, c("SEM", "SLM", "SDM"), lty=1:3)
}

}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
