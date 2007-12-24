simulate.formula <- function(object, nsim=1, seed=NULL, ...,theta0,
                             burnin=1000, interval=1000,
                             basis=NULL,
                             sequential=TRUE,
                             constraints=~.,
                             control=simulate.formula.control(),
                             verbose=FALSE) {
  out.list <- list()
  out.mat <- numeric(0)
  formula <- object
  
  if(!is.null(seed)) set.seed(as.integer(seed))
  if(!is.null(basis)) {
    nw <- basis
#   formula <- as.formula(paste(c("nw",as.character(formula)),
#                               collapse=" "))
    formula <- update(formula, nw ~ .)
    object <- formula
  } else {
    nw <- ergm.getnetwork(formula)
  }
  if(class(nw) =="network.series"){
    nw <- nw$networks[[1]]
  }
  nw <- as.network(nw)
  if(!is.network(nw)){
    stop("A network object on the LHS of the formula or via",
         " the 'basis' argument must be given")
  }

# m <- ergm.getmodel(formula, nw, drop=control$drop)
  m <- ergm.getmodel(formula, nw, drop=FALSE)
  MHproposal <- MHproposal(constraints,control$prop.args, nw, m, weights=control$prop.weights,class="c")

  Clist <- ergm.Cprepare(nw, m)
  
  MCMCsamplesize <- 1
  verb <- match(verbose,
                c("FALSE","TRUE", "very"), nomatch=1)-1
  if(missing(theta0)) {
   if(is.bipartite(nw) && MHproposal$name=="CondDegree" && require("networksis", quietly = TRUE)) {
     control$packagenames <- c(control$packagenames,"networksis")
     return(sis.simulate(nw, formula, m, Clist, nsim=nsim,
            control=control, verbose=verbose, ...))
   }else{
    theta0 <- rep(0,Clist$nparam)
    warning("No parameter values given, using Bernouli network\n\t")
   }
  }
  theta0.bdd <- theta0
  theta0.bdd[is.infinite(theta0.bdd)] <- -10000
  
  if(!is.null(seed)) set.seed(as.integer(seed))

  curstats<-summary.statistics.network(object)

  if (verb) {
    cat("Starting",nsim,"MCMC iterations of",burnin+interval*MCMCsamplesize,
        "steps each:\n")
  }
  for(i in 1:nsim){
    Clist <- ergm.Cprepare(nw, m)
    maxedges <- max(2000, Clist$nedges)
    if(i==1 | !sequential){
      use.burnin <- burnin
    }else{
      use.burnin <- interval
    }
#
#   Check for truncation of the returned edge list
#
    z <- list(newnwheads=maxedges+1)
    while(z$newnwheads[1] > maxedges){
     maxedges <- 10*maxedges
     if (verb) {
       cat(paste("# ", i, " of ", nsim, ": ", sep=""))
     }
     z <- .C("MCMC_wrapper",
             as.integer(Clist$heads), as.integer(Clist$tails), 
             as.integer(Clist$nedges), as.integer(Clist$n),
             as.integer(Clist$dir), as.integer(Clist$bipartite),
             as.integer(Clist$nterms), 
             as.character(Clist$fnamestring),
             as.character(Clist$snamestring), 
             as.character(MHproposal$name),
             as.character(MHproposal$package),
#  Add:  as.double(length(MHproposal$args)), as.double(MHproposal$args), 
             as.double(Clist$inputs),
             as.double(theta0.bdd),
             as.integer(MCMCsamplesize),
             s = as.double(rep(curstats,MCMCsamplesize)),
             as.integer(use.burnin), as.integer(interval), 
             newnwheads = integer(maxedges),
             newnwtails = integer(maxedges), 
             as.integer(verb),
             as.integer(MHproposal$bd$attribs), 
             as.integer(MHproposal$bd$maxout), as.integer(MHproposal$bd$maxin), as.integer(MHproposal$bd$minout), 
             as.integer(MHproposal$bd$minin), as.integer(MHproposal$bd$condAllDegExact),
             as.integer(length(MHproposal$bd$attribs)), 
             as.integer(maxedges), 
             as.integer(0.0), as.integer(0.0), 
             as.integer(0.0),
             PACKAGE="ergm")
    }
#
#   Next update the network to be the final (possibly conditionally)
#   simulated one
#
    out.list[[i]] <- newnw.extract(nw,z)
    curstats <- z$s
    out.mat <- rbind(out.mat,curstats)
    if(sequential){
      nw <-  out.list[[i]]
    }
  }
  if(nsim > 1){
    out.list <- list(formula = formula, networks = out.list, 
                     stats = out.mat, coef=theta0)
    class(out.list) <- "network.series"
  }else{
    out.list <- out.list[[1]]
  }
  return(out.list)
}


simulate.ergm <- function(object, nsim=1, seed=NULL, ..., theta0=NULL,
                          burnin=1000, interval=1000, 
                          sequential=TRUE, 
                          constraints=NULL,
                          control=simulate.ergm.control(),
                          verbose=FALSE) {
  out.list <- vector("list", nsim)
  out.mat <- numeric(0)
  
#  if(missing(multiplicity) & !is.null(object$multiplicity)){
#    multiplicity <- object$multiplicity
#  }
  if(!is.null(seed)) set.seed(as.integer(seed))
  
  nw <- object$network  
  
# m <- ergm.getmodel(object$formula, nw, drop=control$drop)
  m <- ergm.getmodel(object$formula, nw, drop=FALSE)
  ## By default, all arguments but the first are NULL, and all the information is borrowed from the fit.
  MHproposal <- MHproposal(object,constraints=constraints,arguments=control$prop.args, nw=nw, model=m, weights=control$prop.weights)
  MCMCsamplesize <- 1
  verb <- match(verbose,
                c("FALSE","TRUE", "very"), nomatch=1)-1
# multiplicity.constrained <- 1  
  if(missing(theta0))
    theta0 <- object$coef
  eta0 <- ergm.eta(theta0, m$etamap)
  if (verb) {
    cat("Starting",nsim,"MCMC iterations of",burnin+interval*MCMCsamplesize,
        "steps each:\n")
  }
  for(i in 1:nsim){
    Clist <- ergm.Cprepare(nw, m)
    maxedges <- max(5000, Clist$nedges)
    if(i==1 | !sequential){
      use.burnin <- burnin
    }else{
      use.burnin <- interval
    }
    
#
#   Check for truncation of the returned edge list
#
    z <- list(newnwheads=maxedges+1)
    while(z$newnwheads[1] > maxedges){
     maxedges <- 10*maxedges
     if (verb) {
       cat(paste("#", i, " of ", nsim, ": ", sep=""))
     }
     z <- .C("MCMC_wrapper",
            as.integer(Clist$heads), as.integer(Clist$tails), 
            as.integer(Clist$nedges), as.integer(Clist$n),
            as.integer(Clist$dir), as.integer(Clist$bipartite),
            as.integer(Clist$nterms), 
            as.character(Clist$fnamestring),
            as.character(Clist$snamestring), 
            as.character(MHproposal$name),
            as.character(MHproposal$package),
#  Add:  as.double(length(MHproposal$args)), as.double(MHproposal$args), 
            as.double(Clist$inputs),
            as.double(eta0),
            as.integer(MCMCsamplesize),
            s = double(MCMCsamplesize * Clist$nparam),
            as.integer(use.burnin), as.integer(interval), 
             newnwheads = integer(maxedges),
             newnwtails = integer(maxedges), 
            as.integer(verb),
            as.integer(MHproposal$bd$attribs), 
            as.integer(MHproposal$bd$maxout), as.integer(MHproposal$bd$maxin),
            as.integer(MHproposal$bd$minout), 
            as.integer(MHproposal$bd$minin), as.integer(MHproposal$bd$condAllDegExact),
            as.integer(length(MHproposal$bd$attribs)), 
            as.integer(maxedges), 
            as.integer(0.0), as.integer(0.0), 
            as.integer(0.0),
            PACKAGE="ergm")
    }
    #
    #   summarize stats
    if(control$summarizestats){
      class(Clist) <- "networkClist"
      if(i==1){
        globalstatsmatrix <- summary(Clist)
        statsmatrix <- matrix(z$s, MCMCsamplesize, Clist$nparam, byrow = TRUE)
        colnames(statsmatrix) <- m$coef.names
      }else{
        globalstatsmatrix <- rbind(globalstatsmatrix, summary(Clist))
        statsmatrix <- rbind(statsmatrix,
                             matrix(z$s, MCMCsamplesize,
                                    Clist$nparam, byrow = TRUE))
      }
    }
    #
    #   Next update the network to be the final (possibly conditionally)
    #   simulated one

    out.list[[i]] <- newnw.extract(nw, z)
    out.mat <- rbind(out.mat,z$s[(1):(Clist$nparam)])
    if(sequential){
      nw <-  out.list[[i]]
    }
  }
  if(nsim > 1){
    out.list <- list(formula = object$formula, networks = out.list, 
                     stats = out.mat, coef=theta0)
    class(out.list) <- "network.series"
  }else{
    out.list <- out.list[[1]]
  }
  if(control$summarizestats){
    colnames(globalstatsmatrix) <- colnames(statsmatrix)
    print(globalstatsmatrix)
    print(apply(globalstatsmatrix,2,summary.statsmatrix.ergm),scipen=6)
    print(apply(statsmatrix,2,summary.statsmatrix.ergm),scipen=6)
  }
  return(out.list)
}

