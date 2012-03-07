#  File ergm/R/ergm.sociality.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
sociality <- function(object, ...)
UseMethod("sociality")



sociality.default <- function(object,...)
{
  stop("Either a network, an ergm object or a formula argument must be given")
}



sociality.network <- function (object, ..., 
   statistics=NULL){
  require(sna, quietly=TRUE, warn.conflicts=FALSE)
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
  degreecent <- prestige(smatrix,gmode=gmode,cmode="indegree")
  odeg <- order(-degreecent)[c(1,network.size(object))]
  stats <- matrix(0,ncol=length(statistics),nrow=network.size(object))
  for(i in seq(along=statistics)){
    if(statistics[i] %in% "eigenvector.sym"){
      stats[,i] <- prestige(symmatrix,gmode=gmode,cmode="eigenvector")
    }else{
      stats[,i] <- Re(prestige(smatrix,gmode=gmode,cmode=statistics[i]))
    }
    if(diff(stats[odeg,i])>0){stats[,i] <- -stats[,i]}
  }
  colnames(stats) <- statnames
  rownames(stats) <- network.vertex.names(object)
  stats
}




sociality.formula <- function (formula, ..., init, nsim=100,
                               burnin=100, interval=100,
                               constraints=~.,
                               prop.weights="default",
                               prop.args=list(),
                               seed=NULL,  drop=FALSE,
                               statistics=NULL
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
                                control=control.simulate.ergm(
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



sociality.ergm <- function (object, ..., nsim=100,
                            burnin=100, interval=100,
                            constraints=NULL, prop.weights="default", prop.args =list(),
                            seed=NULL, drop=FALSE,
                            statistics=NULL) {

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
