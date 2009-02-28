#  File ergm/R/simulate.ergm.R
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
# Copyright 2003 Mark S. Handcock, University of Washington
#                David R. Hunter, Penn State University
#                Carter T. Butts, University of California - Irvine
#                Steven M. Goodreau, University of Washington
#                Martina Morris, University of Washington
# Copyright 2007 The statnet Development Team
######################################################################
simulate.formula <- function(object, nsim=1, seed=NULL, ...,theta0,
                             burnin=1000, interval=1000,
                             basis=NULL,
                             statsonly=FALSE,
                             sequential=TRUE,
                             constraints=~.,
                             control=control.simulate.formula(),
                             verbose=FALSE) {
  if (!statsonly) {
    nw.list <- list()
  }
  out.mat <- numeric(0)
  formula <- object
  
  if(!is.null(seed)) set.seed(as.integer(seed))
  if(!is.null(basis)) {
    nw <- basis
#   formula <- as.formula(paste(c("nw",as.character(formula)),
#                               collapse=" "))
    formula <- safeupdate.formula(formula, nw ~ .)
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
  MHproposal <- MHproposal(constraints,arguments=control$prop.args,
                           nw=nw, model=m, weights=control$prop.weights, class="c")
  
  Clist <- ergm.Cprepare(nw, m)
  
  verb <- match(verbose,
                c("FALSE","TRUE", "very"), nomatch=1)-1
  if(missing(theta0)) {
    theta0 <- rep(0,Clist$nstats)
    warning("No parameter values given, using Bernouli network\n\t")
  }
  eta0 <- theta0
  eta0[is.infinite(eta0)] <- -10000
  
  if(!is.null(seed)) set.seed(as.integer(seed))
    
# curstats<-summary.statistics.network(object)
  curstats<-summary(safeupdate.formula(object,nw ~ .))
  MCMCparams <- list(samplesize=1,
                     stats=curstats,
                     burnin=burnin,
                     interval=interval,
                     parallel=control$parallel,
                     packagenames=control$packagenames,
                     Clist.miss=ergm.design(nw, m, verbose=verbose))
  
  if (verb) {
    cat(paste("Starting ",nsim," MCMC iteration", ifelse(nsim>1,"s",""),
              " of ", burnin+interval*(MCMCparams$samplesize-1), 
              " steps", ifelse(nsim>1, " each", ""), ".\n", sep=""))
  }
  if(sequential){
    for(i in 1:nsim){
      Clist <- ergm.Cprepare(nw, m)
      maxedges <- max(2000, Clist$nedges)
      if(i==1 | !sequential){
        MCMCparams$burnin <- burnin
      }else{
        MCMCparams$burnin <- interval
      }
      #
      #   Check for truncation of the returned edge list
      #
      z <- list(newnwheads=maxedges+1)
      while(z$newnwheads[1] > maxedges){
        maxedges <- 10*maxedges
        z <- ergm.mcmcslave(Clist,MHproposal,eta0,MCMCparams,maxedges,verb) 
      }
      #
      #   Next update the network to be the final (possibly conditionally)
      #   simulated one
      #
      if (!statsonly) {
        nw.list[[i]] <- newnw.extract(nw,z)
      }
      curstats <- z$s
      names(curstats) <- m$coef.names
      out.mat <- rbind(out.mat,curstats)
      if (sequential){
        if (!statsonly)
          nw <-  nw.list[[i]]
        else 
          nw <- newnw.extract(nw, z)
        MCMCparams$stats<-curstats
      }
    }
  }else{
    stop("Parallelization not currently enabled.")
  }
  if(nsim > 1){
    rownames(out.mat) <- NULL
    if (statsonly) {
      out.list <- as.matrix(out.mat)
    } else {
      out.list <- list(formula = formula, networks = nw.list, 
                       stats = out.mat, coef=theta0)
      class(out.list) <- "network.series"
    }
  } else if (statsonly) {
    out.list <- as.vector(out.mat)
    names(out.list) <- colnames(out.mat)
  } else {
    out.list <- nw.list[[1]]
  }
  return(out.list)
}


simulate.ergm <- function(object, nsim=1, seed=NULL, ..., theta0=NULL,
                          burnin=1000, interval=1000, 
                          statsonly=FALSE,
                          sequential=TRUE, 
                          constraints=NULL,
                          control=control.simulate.ergm(),
                          verbose=FALSE) {
  if (!statsonly) {
    nw.list <- vector("list", nsim)
  }
  out.mat <- numeric(0)
  
#  if(missing(multiplicity) & !is.null(object$multiplicity)){
#    multiplicity <- object$multiplicity
#  }
  if(!is.null(seed)) set.seed(as.integer(seed))
  
  nw <- object$newnetwork  
  
# m <- ergm.getmodel(object$formula, nw, drop=control$drop)
  m <- ergm.getmodel(object$formula, nw, drop=FALSE)
  ## By default, all arguments but the first are NULL, and all the information is borrowed from the fit.
  MHproposal <- MHproposal(object,
      constraints=constraints,
      arguments=control$prop.args, 
      nw=nw, model=m, weights=control$prop.weights)
  verb <- match(verbose,
                c("FALSE","TRUE", "very"), nomatch=1)-1
# multiplicity.constrained <- 1  
  if(missing(theta0)){theta0 <- object$coef}
  eta0 <- ergm.eta(theta0, m$etamap)

  MCMCparams <- list(samplesize=1,
      stats=summary(safeupdate.formula(object$formula,nw ~ .)),
      burnin=burnin,
      interval=interval,
      parallel=control$parallel,
      packagenames=control$packagenames,
      Clist.miss=ergm.design(nw, m, verbose=verbose))

  if (verb) {
    cat(paste("Starting ",nsim," MCMC iteration", ifelse(nsim>1,"s",""),
        " of ", burnin+interval*(MCMCparams$samplesize-1), 
        " steps", ifelse(nsim>1, " each", ""), ".\n", sep=""))
  }
  if(sequential){
    for(i in 1:nsim){
      Clist <- ergm.Cprepare(nw, m)
      maxedges <- max(5000, Clist$nedges)
      if(i==1 | !sequential){
        MCMCparams$burnin <- burnin
      }else{
        MCMCparams$burnin <- interval
      }
#
#   Check for truncation of the returned edge list
#
      z <- list(newnwheads=maxedges+1)
      while(z$newnwheads[1] > maxedges){
        maxedges <- 10*maxedges
        if (verb) {
          cat("   ")
          #cat(paste("  #", i, " of ", nsim, ": ", sep=""))
        }
        z <- ergm.mcmcslave(Clist,MHproposal,eta0,MCMCparams,maxedges,verb) 
      }
      #   summarize stats
      if(control$summarizestats){
        class(Clist) <- "networkClist"
        if(i==1){
          globalstatsmatrix <- summary(Clist)
          statsmatrix <- matrix(z$s, MCMCparams$samplesize, Clist$nstats, byrow = TRUE)
          colnames(statsmatrix) <- m$coef.names
        }else{
          globalstatsmatrix <- rbind(globalstatsmatrix, summary(Clist))
          statsmatrix <- rbind(statsmatrix,
                               matrix(z$s, MCMCparams$samplesize,
                                      Clist$nstats, byrow = TRUE))
        }
      }
      #
      #   Next update the network to be the final (possibly conditionally)
      #   simulated one
      if (!statsonly) {
        nw.list[[i]] <- newnw.extract(nw, z)
      }
      curstats <- z$s[(1):(Clist$nstats)]
      names(curstats) <- m$coef.names
      out.mat <- rbind(out.mat,curstats)
      if(sequential){
        if (!statsonly) 
          nw <-  nw.list[[i]]
        else 
          nw <- newnw.extract(nw, z)
        MCMCparams$stats<-curstats
      }
    }
  }else{
    stop("Parallelization not currently enabled.")   
  }
  if(nsim > 1){
    rownames(out.mat) <- NULL
    if (statsonly) {
      out.list <- as.matrix(out.mat)
    } else {
      out.list <- list(formula = object$formula, networks = nw.list, 
                       stats = out.mat, coef=theta0)
      class(out.list) <- "network.series"
    }
  } else if (statsonly) {
    out.list <- as.vector(out.mat)
    names(out.list) <- colnames(out.mat)
  } else {
    out.list <- nw.list[[1]]
  }
  return(out.list)
}
