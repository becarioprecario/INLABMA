\name{trIrhoWinv}
\alias{trIrhoWinv}

\title{
Compute trace of (I-rho*W)^{-1} matrix
}
\description{
This function computes (or estimates) the trace of matrix (I-rho*W)^{-1},
which is often needed when computing impacts in some spatial econometrics 
models.
}
\usage{
trIrhoWinv(W, rho, offset = 0, order = 20, direct = TRUE, Df = Matrix::Diagonal(nrow(W)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{W}{
Adjacency matrix. Usually, it is row-standardised.
}
  \item{rho}{
Value of spatial autocorrelation parameter \code{rho}.
}
  \item{offset}{
Number of times (I-rho*W)^{-1} is multiplied by W (for sdm model).
}
  \item{order}{
Order of Taylor expansion used in the approximation of the trace.
}
  \item{direct}{
Use direct method, i.e., matrix multiplication, etc.
}
  \item{Df}{
Diagonal matrix used to compute the impacts in the Probit model
only used if direct=TRUE.
}
}
\details{
This function computes the trace of (I-rho*W)^{-1}, which is later used
to computed the impacts. This is an internal function.
}
\value{
Numerica value.
}
\references{
LeSage and Page (2008)
Bivand et al. (2013)
}
\author{
Virgilio Gómez-Rubio <virgilio.gomez@uclm.es>
}
%%\note{
%%  ~~further notes~~
%%}

\seealso{
\code{\link{sem.inla}}, \code{\link{slm.inla}}, \code{\link{sdm.inla}}
}
%%\examples{
%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{math}
