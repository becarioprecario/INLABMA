\name{mysplinefun}
\alias{mysplinefun}
\title{
Compute spline function
}
\description{
This function is similar to \code{splinefun} but it returns 0
outside the range of \code{x}.
}
\usage{
mysplinefun(x, y = NULL, method = c("fmm", "periodic", "natural", "monoH.FC")[1],
   ties = mean)
}
\arguments{
  \item{x}{
x-values to use in the interpolation.
}
  \item{y}{
y-values to use in the interpolation (optional).
}
  \item{method}{
	Method used to compute the spline. See \code{\link{splinefun}} for details.
}
  \item{ties}{
Handling of tied 'x' values. See \code{\link{splinefun}} for details.
}
}
\details{
This function calls \code{\link{splinefun}} and returns a function
with the fitted spline. The main difference is that this new function
returns 0 outside the range of 0.
}
\value{
Returns a function with \code{x} and \code{deriv} arguments.  See \code{\link{splinefun}} for details.
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
\code{\link{splinefun}}
}
%%\examples{
%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{smooth}
