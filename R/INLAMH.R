#d: data to be passed to inla() in a suitable format
#fit.inla: function that fits a model given d and a. Takes two arguments: 'data'
#  for the data and 'b' for the  
#b.init: Initial values of the model parameters
#rq: sample from the proposal distribution q.
#dq: density of the proposal distribution
#prior: prior on the model parameters
#n.sim: Total number of simulations
#n.burin: Simulations for burn in
#n.thin: Thinning AFTER burnin
#n.errors: Number of total INLA errors allowed
#verbose: Whether to plot some information messages (default to FALSE).

INLAMH <- function(d, fit.inla, b.init, rq, dq, prior, n.sim = 200, 
  n.burnin = 100, n.thin = 1, n.errors = 20, verbose = FALSE) {


b.sim <- as.list(rep(NA, n.sim))
model.sim <- as.list(rep(NA, n.sim))
acc.sim <- rep(NA, n.sim)

#Initial values
#b.sim[[1]] <- b.init
#model.sim[[1]] <- fit.inla(d, b.sim[[1]])

b.cur <- b.init
model.cur <- fit.inla(d, b.cur)

#Index for saved iterations
save.idx <- 0

#Total number of simulations
n.sim.tot <- n.burnin + n.thin * n.sim

#Set number of errors to zero
n.err.idx <- 0

i <- 1
while (i <= n.sim.tot) {

   #Sample new proposal
   b.new <- rq(b.cur)
   #Fit model using try to handle possible errors in INLA
   model.new <- try(fit.inla(d, b.new))

   if(class(model.new) == "try-error") {
     #Update number of errors
     n.err.idx <- n.err.idx + 1

     if(verbose) {
        cat("INLA error number ", n.err.idx, " at iteration ", i, ".\n")
     }

     if(n.err.idx > n.errors) {
       stop("INLA failed too many times to fit a model.")
     }


     #Set current to random previous state (it is from posterior)
     if(i > n.burnin & save.idx > 100) {
       idx.aux <- sample(save.idx - 1, 1) #i - 1 to avoid previous state
       b.cur <- b.sim[[idx.aux]]
       model.cur <- model.sim[[idx.aux]]
     } else {
       stop("INLA failed early in the sampling process.")
     }

     #DO NOT incrase 'i'

  } else {

  
   #Compute alpha (log-scale of densities)
   alpha <- model.new$mlik + prior(b.new) + dq(b.new, b.cur)
   alpha <- alpha - ( model.cur$mlik + prior(b.cur) + dq(b.cur, b.new) )
  
   alpha <- min(log(1), alpha)
   alpha <- exp(alpha)

   acc.sim[i] <- runif(1) < alpha

   if(acc.sim[i]) { #Accept
      b.cur <- b.new
      model.cur <- model.new
   }

   #Save results (if needed)
   if(i > n.burnin & (i - n.burnin + 1 ) %% n.thin == 0) {
     save.idx <- save.idx + 1
     b.sim[[save.idx]] <- b.cur
     model.sim[[save.idx]] <- model.cur
   }

   if(verbose) {
    if(i %% 100 ==0) {cat("Iteration ", i, "completed. Time:", 
      as.character(Sys.time()), "\n")}
   }

   #INCREASE 'i'
   i <- i + 1 
  }#Try-error

}

  res <- list(acc.sim = acc.sim, model.sim = model.sim,
    b.sim = b.sim)

  return(res)
}
