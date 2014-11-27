#  File R/mvmodel.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################
#=============================================================
# This file contains the following 4 files for gathering
# summary statistics from simulated networks:
#      <mvmodel>           <mvmodel.formula>
#      <mvmodel.default>   <mvmodel.ergm>
#=============================================================



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
#
# --PARAMETERS--
#   formula/object: X, either a formula or ergm
#   ...           : any parameters passed via ... are ignored
#   init        : the vector of initial theta values to use when
#                   simulating networks
#   nsim          : the number of simulations to draw
#   burnin        : the number of proposals to disregard for the initial
#                   burn-in period
#   interval      : the number of proposals to disregard in between
#                   the drawn simulations
#   constraints   : a one-sided formula specifying the constraints on the
#                   support of the distribution of networks being simulated;
#                   default=NULL
#   control       : a list of control parameters for algorithm tuning, as
#                   returned by <control.simulate.ergm> or
#                   <control.simulate.formula>; default=<control.simulate.X>
#   seed          : an integer at which to set the random number generator
#   statistic     : this parameter may have one of two forms - either
#                   a function that accepts a network 'nw' as input and
#                   as output gives a summary statistic, or this may be
#                   the character string "density", in which case, the
#                   'statistic' function is defined for you; the
#                   default is to return the sociomatrix of 'nw'
#
# --RETURNED--
#   if 'statistic' takes on its default value, a "sociomatrix" is
#      returned, where the ij entry is the percent of simulations in which
#      edge ij was present;
#   otherwise, the simulation summary is returned as a list containing:
#     mean:  the vector of mean 'statistic's over the set of simulations
#     sim :  the matrix of summary statistics where entry ij is the
#            jth statistic returned by 'statistic' for the ith simulation
#
########################################################################

mvmodel.formula <- function (formula, ..., init, nsim=100,
                             burnin=10000, interval=1000,
                             constraints=NULL,
                             control=control.simulate.formula(),
                             seed=NULL, 
                             statistic=NULL
		      ) {
  check.control.class("simulate.formula", "ERGM mvmodel.formula")
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

  simcentrality <- statistic(SimGraphSeriesObj[[1]])

  if(!probabilites){
    sim.mvmodel <-array(0,dim=c(nsim,length(simcentrality)))
    dimnames(sim.mvmodel) <- list(paste(c(1:nsim)),names(simcentrality))
    sim.mvmodel[1,] <- simcentrality
  }else{
    sim.mvmodel <- simcentrality
  }

  for (i in 2:nsim)
  { 
    simcentrality <- statistic(SimGraphSeriesObj[[i]])
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




#######################################################################
#
#  (see the <mvmodel.formula> function for header details)
#
########################################################################

mvmodel.ergm <- function (object, ..., nsim=100,
                          burnin=10000, interval=1000,
                          constraints=NULL,
                          seed=NULL,
                          control=control.simulate.ergm(),
                          statistic=NULL) {

  check.control.class("simulate.ergm")
  
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

  simcentrality <- statistic(SimGraphSeriesObj[[1]])

  if(!probabilites){
    sim.mvmodel <-array(0,dim=c(nsim,length(simcentrality)))
    dimnames(sim.mvmodel) <- list(paste(c(1:nsim)),names(simcentrality))
    sim.mvmodel[1,] <- simcentrality
  }else{
    sim.mvmodel <- simcentrality
  }

  for (i in 2:nsim)
  { 
    simcentrality <- statistic(SimGraphSeriesObj[[i]])
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
