\name{logprrho}
\alias{logprrho}

\title{
Log-prior density for the spatial autocorrelation parameter \code{rho}
}
\description{
Compute log-prior density for \code{rho}
}
\usage{
logprrho(rho)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rho}{
The value to compute the log-density.
}
}
\details{
This function computes the log-density of the prior for \code{rho}
according to logit(rho) ~ N(0, prec=.1). THis is one of the default
priors in \pkg{R-INLA} for spatial autocorrelation parameters.
}
\value{
Numerical.
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

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{

rrho<-seq(.01, .99, length.out=100)
plot(rrho, exp(logprrho(rrho)))
}


\keyword{distribution}
