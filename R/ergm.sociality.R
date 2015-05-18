#  File R/ergm.sociality.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2015 Statnet Commons
#######################################################################
#==================================================================
# This file contains the 5 following functions for ??
#    <sociality>            <sociality.network>
#    <sociality.default>    <sociality.formula>
#    <sociality.ergm>
#==================================================================



sociality <- function(object, ...)
UseMethod("sociality")



sociality.default <- function(object,...)
{
  stop("Either a network, an ergm object or a formula argument must be given")
}




#################################################################
# The <sociality.network> function ??
#
# --PARAMETERS--
#   object    : a network object
#   ...       : any parameters passed via ... are ignored
#   statistics: is a character vector naming the columns of the
#               returned matrix; default="centrality (eigen)" for
#               undirected networks and c("prestige (eigen)",
#               centrality (eigen)") for undirected networks
#
# --RETURNED--
#   stats:  a matrix whose i,j entry gives the ?? of the ith
#           statistic of 'statistics' for the jth node of 'object'
#
#################################################################

sociality.network <- function (object, ..., 
   statistics=NULL){
  .Deprecated(msg="the sociality.network function will not be supported in the future. see summary.formula and the ergm term 'sociality' for an alternate")
  requireNamespace('sna', quietly=TRUE, warn.conflicts=FALSE)
  if(!is.directed(object)){
    gmode <- "network"
    if(is.null(statistics)){
     statistics <- c("eigenvector")
    }
    statnames <- statistics
    statnames[statnames=="eigenvector"] <- "centrality (eigen)"
  }else{
    gmode <- "dinetwork"
    if(is.null(statistics)){
     statistics <- c("eigenvector","eigenvector.sym")
    }
    statnames <- statistics
    statnames[statnames=="eigenvector"] <- "prestige (eigen)"
    statnames[statnames=="eigenvector.sym"] <- "centrality (eigen)"
  }
  smatrix <- as.sociomatrix(object)
  symmatrix <- smatrix 
  symmatrix[t(smatrix) > 0] <- 1
  degreecent <- sna::prestige(smatrix,gmode=gmode,cmode="indegree")
  odeg <- order(-degreecent)[c(1,network.size(object))]
  stats <- matrix(0,ncol=length(statistics),nrow=network.size(object))
  for(i in seq(along=statistics)){
    if(statistics[i] %in% "eigenvector.sym"){
      stats[,i] <- sna::prestige(symmatrix,gmode=gmode,cmode="eigenvector")
    }else{
      stats[,i] <- Re(sna::prestige(smatrix,gmode=gmode,cmode=statistics[i]))
    }
    if(diff(stats[odeg,i])>0){stats[,i] <- -stats[,i]}
  }
  colnames(stats) <- statnames
  rownames(stats) <- network.vertex.names(object)
  stats
}




###########################################################################
# The <sociality.formula> function ??
#
# --PARAMETERS--
#   formula     :  a formula of the form 'nw ~ model term(s)'
#   ...         : any parameters passed via ... are ignored
#   init      : the vector of initial theta values
#   nsim        : the number of simulations to gather for the
#                 returned 'sim.sociality' vector
#   burnin      : the number of proposals to ignore before MCMC sampling
#                 begins; default=10,000
#   interval    : the number of proposals to disregard between sampled 
#                 statistics; default=100
#   constraints : a one-sided formula of the constraint terms; options are
#                      bd        degrees        nodegrees
#                      edges     degreedist     idegreedist
#                      observed  odegreedist
#                 default="~ ."
#   prop.weights: the method to allocate probabilities of being proposed
#                 to dyads as "TNT", "random", "nonobserved", or "default"
#                 default="default", which is based upon the ergm constraints
#   prop.args   : an alternative, direct way of specifying additional
#                 arguments to proposal              
#   seed        :  an integer starting value for the random number generator;
#                  default=NULL
#   drop        : whether degenerate terms should be dropped from the fit 
#                 (T or F); default=TRUE  
#   statistics  : is a character vector naming the statistics to ??
#
# --RETURNED--
#   the ?? as an invisible list containing:
#     sociality    :
#     sim.sociality:
#
###############################################################################

sociality.formula <- function (formula, ..., init, nsim=100,
                               burnin=100, interval=100,
                               constraints=~.,
                               prop.weights="default",
                               prop.args=list(),
                               seed=NULL,  drop=FALSE,
                               statistics=NULL
                               ) {
  .Deprecated(msg="the sociality.formula function will not be supported in the future. see summary.formula and the ergm term 'sociality' for an alternate")
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

  m <- ergm.getmodel(trms, g, drop=drop)
  Clist <- ergm.Cprepare(g, m)

  if(missing(init)){
      init <- rep(0,Clist$nstats)
      warning("No parameter values given, using 0\n\t")
  }

  n <- network.size(g)
  if(is.null(seed)){seed <- sample(10000000, size=1)}

  if(!is.directed(g)){
    if(is.null(statistics)){
     dimcentrality <- 1
    }else{
     dimcentrality <- length(statistics)
    }
  }else{
    if(is.null(statistics)){
     dimcentrality <- 2
    }else{
     dimcentrality <- length(statistics)
    }
  }

  # Set up the output arrays

  sim.sociality <-array(0,dim=c(nsim,n,dimcentrality))

  dimnames(sim.sociality) <- list(paste(c(1:nsim)),paste(1:n),rep("",dimcentrality))

  # Simulate an exponential family random network model

  SimGraphSeriesObj <- simulate(formula, nsim=nsim, seed=seed,
                                burnin=burnin, interval=interval,
                                constraints=constraints, coef=init,
                                control=control.simulate.formula(
                                  MCMC.burnin=burnin,
                                  MCMC.interval=interval,
                                  MCMC.prop.args=prop.args,
                                  MCMC.prop.weights=prop.weights))

# cat("\nCollating simulations\n")

  for (i in 1:nsim)
  { 
    simcentrality <- sociality(SimGraphSeriesObj$networks[[i]],
                                     statistics=statistics)
    sim.sociality[i,,] <- simcentrality
  }

  dimnames(sim.sociality)[[2]] <- rownames(simcentrality)
  dimnames(sim.sociality)[[3]] <- colnames(simcentrality)
  
  simcentrality <- apply(sim.sociality,c(2,3),mean)
  invisible(list(sociality=simcentrality,
                 sim.sociality=sim.sociality))
  }




###########################################################################
# The <sociality.ergm> function ??
#
# --PARAMETERS--
#   object      :  a formula of the form 'nw ~ model term(s)'
#   ...         : any parameters passed via ... are ignored
#   nsim        : the number of simulations to gather for the
#                 returned 'sim.sociality' vector
#   burnin      : the number of proposals to ignore before MCMC sampling
#                 begins; default=10,000
#   interval    : the number of proposals to disregard between sampled 
#                 statistics; default=100
#   constraints : a one-sided formula of the constraint terms; options are
#                      bd        degrees        nodegrees
#                      edges     degreedist     idegreedist
#                      observed  odegreedist
#                 default="~ ."
#   prop.weights: the method to allocate probabilities of being proposed
#                 to dyads as "TNT", "random", "nonobserved", or "default"
#                 default="default", which is based upon the ergm constraints
#   prop.args   : an alternative, direct way of specifying additional
#                 arguments to proposal              
#   seed        :  an integer starting value for the random number generator;
#                  default=NULL
#   drop        : whether degenerate terms should be dropped from the fit 
#                 (T or F); default=TRUE  
#   statistics  : is a character vector naming the statistics to ??
#
# --RETURNED--
#   the ?? as an invisible list containing:
#     sociality    :
#     sim.sociality:
#
###############################################################################

sociality.ergm <- function (object, ..., nsim=100,
                            burnin=100, interval=100,
                            constraints=NULL, prop.weights="default", prop.args =list(),
                            seed=NULL, drop=FALSE,
                            statistics=NULL) {
  .Deprecated(msg="the sociality.ergm function will not be supported in the future. see summary.formula and the ergm term 'sociality' for an alternate")

# trms <- ergm.getterms(object$formula)
# g <- as.network(eval(trms[[2]], sys.parent()))
  g <- object$network
  g <- as.network(g)
  if(!is.network(g)){
    stop("A network object must be given")
  }
  n <- network.size(g)
  if(is.null(seed)){seed <- sample(10000000, size=1)}
     
  if(!is.directed(g)){
    if(is.null(statistics)){
     dimcentrality <- 1
    }else{
     dimcentrality <- length(statistics)
    }
  }else{
    if(is.null(statistics)){
     dimcentrality <- 2
    }else{
     dimcentrality <- length(statistics)
    }
  }

  # Set up the output arrays

  sim.sociality <-array(0,dim=c(nsim,n,dimcentrality))

  dimnames(sim.sociality) <- list(paste(c(1:nsim)),paste(1:n),rep("",dimcentrality))

  # Simulate an exponential family random network model

  SimGraphSeriesObj <- simulate(object, nsim=nsim, seed=seed,
                                constraints=constraints,
                                control=control.simulate.ergm(
                                   MCMC.burnin=burnin,
                                   MCMC.interval=interval,
                                   MCMC.prop.weights=prop.weights,
                                   MCMC.prop.args=prop.args))

# cat("\nCollating simulations\n")

  for (i in 1:nsim)
  { 
    simcentrality <- sociality(SimGraphSeriesObj$networks[[i]],
                                     statistics=statistics)
    sim.sociality[i,,] <- simcentrality
  }

  dimnames(sim.sociality)[[2]] <- rownames(simcentrality)
  dimnames(sim.sociality)[[3]] <- colnames(simcentrality)
  
  simcentrality <- apply(sim.sociality,c(2,3),mean)
  invisible(list(sociality=simcentrality,
                 sim.sociality=sim.sociality))
}
