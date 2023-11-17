\name{rescalemarg}
\alias{rescalemarg}

\title{
Re-scale marginal distribution to compute the distribution of w*x
}
\description{
This function takes a marginal distribution (represetend by a 2-column matrix)
and computes the marginal distribution of w*x.
}
\usage{
rescalemarg(xx, w)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xx}{
2-column matrix with x and y-values.
}
  \item{w}{
Weight to re-scale the y-values.
}
}
\details{
This function simply re-scales 
}
\value{
A 2-column matrix with the new values of w*x and their associated
probability densities. This is also an object of classes \code{inla.marginal}.
}
\references{
INLA
}
\author{
Virgilio Gómez-Rubio <virgilio.gomez@uclm.es>
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{inla.tmarginal}
}
\examples{

if(requireNamespace("INLA", quietly = TRUE)) {
require(INLA)
x<-seq(-3,3, by=.01)
xx<-cbind(x, dnorm(x))

xx2<-rescalemarg(xx, 3)

plot(xx, type="l", xlim=c(-9,9))
lines(xx2, col="red")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
