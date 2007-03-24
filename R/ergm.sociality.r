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

sociality.formula <- function (formula, ..., theta0, nsim=100,
                      burnin=100, interval=100,
                      proposaltype="randomtoggle", multiplicity=1,
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

  if(missing(theta0)){
      theta0 <- rep(0,Clist$nparam)
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

  SimGraphSeriesObj <- simulate(formula, burnin=burnin, interval=interval,
                             proposaltype=proposaltype,
                             multiplicity=multiplicity,
                             theta0=theta0,
                             drop=drop,
                             n=nsim, seed=seed)

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
                      proposaltype="randomtoggle", multiplicity=1,
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
  if(missing(proposaltype) & !is.null(object$proposaltype)){
    proposaltype <- object$proposaltype
  }

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

  SimGraphSeriesObj <- simulate(object, burnin=burnin, interval=interval,
                             proposaltype=proposaltype,
                             multiplicity=multiplicity,
                             drop=drop,
                             n=nsim, seed=seed)

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
