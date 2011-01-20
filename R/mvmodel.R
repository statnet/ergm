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
#   theta0        : the vector of initial theta values to use when
#                   simulating networks
#   nsim          : the number of simulations to draw
#   burnin        : the number of proposals to disregard for the initial
#                   burn-in period
#   interval      : the number of proposals to disregard in between
#                   the drawn simulations
#   constraints   : a one-sided formula specifying the constraints on the
#                   support of the distribution of networks being simulated;
#                   default=NULL
#   prop.weights  : specifies the method used to allocate probabilities
#                   if being proposed to dyads; options are "TNT",
#                   "random", "nonobserved" and "default"; default=
#                   NULL if X is an ergm (which then uses the weights
#                   that the ergm was fit by); default="default" if
#                   X is a formula (which picks a reasonable default
#                   considering any constraints)
#   prop.args     : an alternative, direct way of specifying additional
#                   arguments to proposal
#   seed          : an integer at which to set the random number generator
#   drop          : whether degenerate terms should be dropped from the
#                   fit (T or F); default=TRUE                             
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

mvmodel.formula <- function (formula, ..., theta0, nsim=100,
                             burnin=100, interval=100,
                             constraints=~., prop.weights="default", prop.args=NULL,
                             seed=NULL,  drop=FALSE,
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

  m <- ergm.getmodel(formula, g, drop=drop)
  Clist <- ergm.Cprepare(g, m)

  if(missing(theta0)){
      theta0 <- rep(0,Clist$nstats)
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

  SimGraphSeriesObj <- simulate(formula, burnin=burnin, interval=interval,
                                constraints=constraints,
                                control=control.simulate.ergm(prop.args=prop.args,
                                  prop.weights=prop.weights,
                                  drop=drop),
                                theta0=theta0,
                                n=nsim, seed=seed)
  
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




#######################################################################
#
#  (see the <mvmodel.formula> function for header details)
#
########################################################################

mvmodel.ergm <- function (object, ..., nsim=100,
                          burnin=100, interval=100,
                          constraints=NULL, prop.weights="default", prop.args=NULL,
                          seed=NULL, drop=FALSE,
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
  constraints <-
    if(is.null(constraints)) object
    else constraints

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

  SimGraphSeriesObj <- simulate(object, burnin=burnin, interval=interval,
                                constraints=constraints,
                                control=control.simulate.ergm(prop.args=prop.args,
                                  prop.weights=prop.weights,
                                  drop=drop),
                             n=nsim, seed=seed)

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
