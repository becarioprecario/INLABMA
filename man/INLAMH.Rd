\name{INLAMH}
\alias{INLAMH}

\title{
Perform INLA with MCMC.
}
\description{
This function implements the Metropolis-Hastings algorithm using repeated
calls to R-INLA to fint conditional model on the current state of the
MCMC simulations.
}
\usage{
INLAMH(d, fit.inla, b.init, rq, dq, prior, n.sim = 200, n.burnin = 100,
  n.thin = 1, n.errors = 20, verbose = FALSE)
}

\arguments{
  \item{d}{Data.frame with the data used to fit the model with R-INLA.
}
  \item{fit.inla}{A function used to fit the model with R-INLA. It should
take at least two arguments: a data.frame (first) and an object with
the actual value of the sampled parameters. This function must return
a vector of two components: model.sim (an 'inla' object with the fitted model)
and 'mlik' (the marginal likelihood as returned by INLA in model.sim$mlik).
}
  \item{b.init}{Initial values of the model parameters for the 
Metropolis-Hastings algorithm.
}
  \item{rq}{Sampling from the proposal distribution. It must take one
argument: the current state of the Markov chain.
}
  \item{dq}{Density of the proposal distribution. It takes two arguments:
current state and proposed new state.
}
  \item{prior}{Prior distribution of the model parameters.
}
  \item{n.sim}{Total of simulations to be done.
}
  \item{n.burnin}{Number of burn-in simulation (thinning is ignored here).
}
  \item{n.thin}{Thinning to be applied to the simulations after burn-in.
}
  \item{n.errors}{This is the number of errores allowed when calling inla().
}
  \item{verbose}{Whether to show some running information or not (defaut to
FALSE).
}
}
\details{
This function implements the Metropolis-Hastings algorithm using INLA
(i.e., INLA within MCMC) at every step. In practice, only a few of the model
parameters are sampled in the MCMC steps and the posterior marginal of the
remainder of parameters is obtained by Bayesian model averaging of the 
conditional marginals returned by R-INLA at each step of the Metropolis-Hastings
algorithm.

}
\value{
A list with three components:
\item{acc.sim}{A vector of logical values (of length 'n.sim') showing whether a given proposal
has been accepted or not. This is useful to compute the acceptance rate.}
\item{model.sim}{A list with the models fitted, as returned by fit.inla().}
\item{b.sim}{List of all sampled values of the models parameters. It is
a list beacuse the sampled values can be vectors.}
}
\references{
Virgilio Gómez-Rubio and Håvard Rue (2017).  Markov Chain Monte Carlo with the Integrated Nested Laplace Approximation. \url{https://arxiv.org/abs/1210.4908}.
}
\author{
Virgilio Gómez-Rubio.
}

%%\examples{
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.
%%}

\keyword{misc}
