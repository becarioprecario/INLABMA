\name{fitmargBMA}
\alias{fitmargBMA}
\alias{fitmargBMA2}
\alias{fitmatrixBMA}
\alias{fitlistBMA}
\title{
Compute marginals using Bayesian Model Averaging
}
\description{
\code{fitmargBMA} takes a list of marginal distributions and weights 
(presumably, based on some marginal likelihoods) and computes a 
final distribution by weighting.

\code{fitmargBMA2} takes a list of INLA models and computes Bayesian
Model Averaging on some of their components.

\code{fitmatrixBMA} performs  averaging on a list of matrices.

\code{fitlistBMA} performs averaging of elements in lists.
}
\usage{
fitmargBMA(margs, ws, len = 100)
fitmargBMA2(models, ws, item)
fitmatrixBMA(models, ws, item)
fitlistBMA(models, ws, item)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{margs}{
List of 2-column matrices with the values of the (marginal) distributions.
}
  \item{models}{
List of INLA models to be averaged.
}
  \item{ws}{
Vector of weights. They do not need to sum up to one.
}
  \item{len}{
Length of the x-vector to compute the weighted distribution.
}
  \item{item}{
Name of the elements of an INLA object to be used in the Model
Averaging.
}
}
\details{
For \code{fitmargBMA}, distributions provided are averaging according to the
weights provided.  A new probability distribution is obtained.  

\code{fitmargBMA2} uses a list of INLA models to compute Model Averaging
on some of their components (for example, the fitted values).

\code{fitmatrixBMA} performs  averaging on a list of matrices.

\code{fitlistBMA} performs  averaging of a list of a list of matrices.
}
\value{
\code{fitmargBMA} returns a 2-column matrix with the weighted marginal 
distribution.

\code{fitmargBMA2} returns a list of weighted components.

\code{fitmatrixBMA} returns a matrix.

\code{fitlistBMA} returns a list.
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
%%\examples{
%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{smooth}
