#  File ergm/R/mvmodel.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
mvmodel <- function(object, ...)
UseMethod("mvmodel")



mvmodel.default <- function(object,...)
{
  stop("Either a network, an ergm object or a formula argument must be given")
}


#######################################################################
# The <mvmodel.X> functions each perform a simulation study, which
# simulates a given number of networks according to a provided formula
# or ergm X and summarizes the given statistics 
########################################################################

mvmodel.formula <- function (formula, ..., init, nsim=100,
                             burnin=10000, interval=1000,
                             constraints=NULL,
                             control=control.simulate.formula(),
                             seed=NULL, 
                             statistic=NULL
		      ) {
  trms <- ergm.getterms(formula)
  if(length(trms)>2){
    g <- eval(trms[[2]], sys.parent())
  }else{
     stop("A network object on the RHS of the formula must be given")
  }
  g <- as.network(g)
  if(!is.network(g)){
    stop("A network object on the RHS of the formula must be given")
  }

  m <- ergm.getmodel(formula, g)
  Clist <- ergm.Cprepare(g, m)

  if(missing(init)){
      init <- rep(0,Clist$nstats)
      warning("No parameter values given, using 0\n\t")
  }

  n <- network.size(g)
  if(is.null(seed)){seed <- sample(10000000, size=1)}

  probabilites <- FALSE
  if(!is.function(statistic)){
   if(is.character(statistic) && statistic=="density"){
    statistic <- function(x){summary(x ~ density)}
   }else{
    statistic <- function(x){as.sociomatrix(x)}
    probabilites <- TRUE
   }
  }

  # Simulate an exponential family random network model

  SimGraphSeriesObj <- simulate(formula, nsim=nsim, seed=seed,
                                constraints=constraints,
                                control=control,
                                coef=init)
  
# cat("\nCollating simulations\n")

  # Set up the output arrays

  simcentrality <- statistic(SimGraphSeriesObj$networks[[1]])

  if(!probabilites){
    sim.mvmodel <-array(0,dim=c(nsim,length(simcentrality)))
    dimnames(sim.mvmodel) <- list(paste(c(1:nsim)),names(simcentrality))
    sim.mvmodel[1,] <- simcentrality
  }else{
    sim.mvmodel <- simcentrality
  }

  for (i in 2:nsim)
  { 
    simcentrality <- statistic(SimGraphSeriesObj$networks[[i]])
    if(!probabilites){
      sim.mvmodel[i,] <- simcentrality
    }else{
      sim.mvmodel <- sim.mvmodel + simcentrality
    }
  }

  if(!probabilites){
    simcentrality <- apply(sim.mvmodel,2,mean)
    return(list(mean=simcentrality,
                   sim=sim.mvmodel))
  }else{
    return(list(mean=sim.mvmodel / nsim))
  }
}




mvmodel.ergm <- function (object, ..., nsim=100,
                          burnin=10000, interval=1000,
                          constraints=NULL,
                          seed=NULL,
                          control=control.simulate.ergm(),
                          statistic=NULL) {

# trms <- ergm.getterms(object$formula)
# g <- as.network(eval(trms[[2]], sys.parent()))
  g <- object$network
  g <- as.network(g)
  if(!is.network(g)){
    stop("A network object must be given")
  }
  n <- network.size(g)
  if(is.null(seed)){seed <- sample(10000000, size=1)}
  if(is.null(constraints)) {constraints <- object$constraints}

  probabilites <- FALSE
  if(!is.function(statistic)){
   if(is.character(statistic) && statistic=="density"){
    statistic <- function(x){summary(x ~ density)}
   }else{
    statistic <- function(x){as.sociomatrix(x)}
    probabilites <- TRUE
   }
  }

  # Simulate an exponential family random network model

  SimGraphSeriesObj <- simulate(object, nsim=nsim, seed=seed,
                                constraints=constraints,
                                control=control)

  # cat("\nCollating simulations\n")

  # Set up the output arrays

  simcentrality <- statistic(SimGraphSeriesObj$networks[[1]])

  if(!probabilites){
    sim.mvmodel <-array(0,dim=c(nsim,length(simcentrality)))
    dimnames(sim.mvmodel) <- list(paste(c(1:nsim)),names(simcentrality))
    sim.mvmodel[1,] <- simcentrality
  }else{
    sim.mvmodel <- simcentrality
  }

  for (i in 2:nsim)
  { 
    simcentrality <- statistic(SimGraphSeriesObj$networks[[i]])
    if(!probabilites){
      sim.mvmodel[i,] <- simcentrality
    }else{
      sim.mvmodel <- sim.mvmodel + simcentrality
    }
  }

  if(!probabilites){
    simcentrality <- apply(sim.mvmodel,2,mean)
    return(list(mean=simcentrality,
                   sim=sim.mvmodel))
  }else{
    return(sim.mvmodel / nsim)
  }
}
