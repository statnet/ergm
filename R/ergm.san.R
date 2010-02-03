san <- function(object, ...){
 UseMethod("san")
}

san.default <- function(object,...)
{
  stop("Either a ergm object or a formula argument must be given")
}

san.formula <- function(object, nsim=1, seed=NULL, ...,theta0=NULL,
                        tau=1, invcov=NULL,
                        burnin=10000, interval=10000,
                        meanstats=NULL,
                        sequential=TRUE,
                        constraints=~.,
                        basis=NULL,
                        control=control.san(),
                        verbose=FALSE) {
  out.list <- list()
  out.mat <- numeric(0)
  formula <- object

  if(!is.null(seed)) set.seed(as.integer(seed))
  if(!is.null(basis)) {
    nw <- basis
    formula <- safeupdate.formula(formula, nw ~ .)
    object <- formula
  } else {
    nw <- ergm.getnetwork(formula)
  }
  if(class(nw) =="network.series"){
    nw <- nw$networks[[1]]
  }
  nw <- as.network(nw)
  if(is.null(meanstats)){
    stop("You need to specify target statistic via",
         " the 'meanstats' argument")
  }
  if(!is.network(nw)){
    stop("A network object on the LHS of the formula ",
         "must be given")
  }

# model <- ergm.getmodel(formula, nw, drop=control$drop)
  model <- ergm.getmodel(formula, nw, drop=FALSE)
  Clist <- ergm.Cprepare(nw, model)
  Clist.miss <- ergm.design(nw, model, verbose=verbose)
  
  MCMCsamplesize <- 1
  verb <- match(verbose,
                c("FALSE","TRUE", "very"), nomatch=1)-1
  MHproposal<-MHproposal(constraints,control$prop.args,nw,model,weights=control$prop.weights)
# if(is.null(theta0)) {
#   warning("No parameter values given, using the MPLE for the passed network.\n\t")
# }
# theta0 <- c(theta0[1],rep(0,Clist$nstats-1))
  
  if(!is.null(seed)) set.seed(as.integer(seed))

  if (verb) {
    cat(paste("Starting ",nsim," MCMC iteration", ifelse(nsim>1,"s",""),
        " of ", burnin+interval*(MCMCsamplesize-1), 
        " steps", ifelse(nsim>1, " each", ""), ".\n", sep=""))
  }

  for(i in 1:nsim){
    Clist <- ergm.Cprepare(nw, model)
#   Clist.miss <- ergm.design(nw, model, verbose=verbose)
    maxedges <- max(20000, Clist$nedges)
    if (verb) {
       cat(paste("#", i, " of ", nsim, ": ", sep=""))
     }
    if(i==1 | !sequential){
      use.burnin <- burnin
    }else{
      use.burnin <- interval
    }

    if(is.null(theta0)) {
      fit <- ergm.mple(Clist=Clist, Clist.miss=Clist.miss, 
                          m=model, verbose=verbose, ...)
      theta0 <- fit$coef
      if(is.null(invcov)) { invcov <- fit$covar }
    }else{
      if(is.null(invcov)) { invcov <- diag(length(theta0)) }
    }
    eta0 <- ergm.eta(theta0, model$etamap)
    
    netsumm<-summary(model$formula)
    if(length(netsumm)!=length(meanstats))
      stop("Incorrect length of the meanstats vector: should be ", length(netsumm), " but is ",length(meanstats),".")

    stats <- matrix(netsumm-meanstats,
                    ncol=Clist$nstats,byrow=TRUE,nrow=MCMCsamplesize)
    tau <- rep(tau,length=length(eta0))
#
#   Check for truncation of the returned edge list
#
    z <- list(newnwheads=maxedges+1)
    eta0[is.na(eta0)]<-0
    while(z$newnwheads[1] > maxedges){
     maxedges <- 10*maxedges
     z <- .C("SAN_wrapper",
             as.integer(Clist$heads), as.integer(Clist$tails), 
             as.integer(Clist$nedges), as.integer(Clist$maxpossibleedges),
             as.integer(Clist$n),
             as.integer(Clist$dir), as.integer(Clist$bipartite),
             as.integer(Clist$nterms), 
             as.character(Clist$fnamestring),
             as.character(Clist$snamestring), 
             as.character(MHproposal$name),
             as.character(MHproposal$package),
             as.double(Clist$inputs),
             as.double(eta0),
             as.double(tau),
             as.integer(MCMCsamplesize),
             s = as.double(stats),
             as.integer(use.burnin), as.integer(interval), 
             newnwheads = integer(maxedges),
             newnwtails = integer(maxedges), 
             as.double(invcov),
             as.integer(verb),
             as.integer(MHproposal$bd$attribs), 
             as.integer(MHproposal$bd$maxout), as.integer(MHproposal$bd$maxin),
             as.integer(MHproposal$bd$minout), as.integer(MHproposal$bd$minin),
             as.integer(MHproposal$bd$condAllDegExact),
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
    out.list[[i]] <- newnw.extract(nw, z)
    out.mat <- rbind(out.mat,z$s[(Clist$nstats+1):(2*Clist$nstats)])
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


san.ergm <- function(object, nsim=1, seed=NULL, ..., theta0=NULL,
                       burnin=100000, interval=10000, 
                       meanstats=NULL,
                       sequential=TRUE, 
                       constraints=NULL,
                       control=control.san(),
                       verbose=FALSE) {
  out.list <- vector("list", nsim)
  out.mat <- numeric(0)
  
#  if(is.null(multiplicity) & !is.null(object$multiplicity)){
#    multiplicity <- object$multiplicity
#  }
  if(!is.null(seed)) set.seed(as.integer(seed))
  
  nw <- object$network  
  
  m <- ergm.getmodel(object$formula, nw, drop=FALSE)
  MCMCsamplesize <- 1
  verb <- match(verbose,
                c("FALSE","TRUE", "very"), nomatch=1)-1
  MHproposal<-MHproposal(object,constraints=constraints,arguments=control$prop.args,nw=nw,model=m,weights=control$prop.weights)
  # multiplicity.constrained <- 1  
  if(is.null(meanstats)){
    stop("You need to specify target statistic via",
         " the 'meanstats' argument")
  }
  if(is.null(theta0)) {
    theta0 <- object$coef
  }
  eta0 <- ergm.eta(theta0, m$etamap)
  if (verb & nsim > 1) {
    cat("Starting",nsim,"MCMC iterations of",burnin+interval*(MCMCsamplesize-1),
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
    
    stats <- matrix(summary(m$formula)-meanstats,
                    ncol=Clist$nstats,byrow=TRUE,nrow=MCMCsamplesize)
#
#   Check for truncation of the returned edge list
#
    z <- list(newnwheads=maxedges+1,newnwtails=maxedges+1)
    while(z$newnwheads[1] > maxedges){
     maxedges <- 10*maxedges
     if (verb & nsim > 1) {
       cat(paste("#", i, " of ", nsim, ": ", sep=""))
     }
     z <- .C("SAN_wrapper",
            as.integer(Clist$heads), as.integer(Clist$tails), 
            as.integer(Clist$nedges), as.integer(Clist$maxpossibleedges),
            as.integer(Clist$n),
            as.integer(Clist$dir), as.integer(Clist$bipartite),
            as.integer(Clist$nterms), 
            as.character(Clist$fnamestring),
            as.character(Clist$snamestring), 
            as.character(MHproposal$name),
            as.character(MHproposal$package),
            as.double(Clist$inputs),
            as.double(eta0),
            as.integer(MCMCsamplesize),
            s = as.double(stats),
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
        statsmatrix <- matrix(z$s, MCMCsamplesize, Clist$nstats, byrow = TRUE)
        colnames(statsmatrix) <- m$coef.names
      }else{
        globalstatsmatrix <- rbind(globalstatsmatrix, summary(Clist))
        statsmatrix <- rbind(statsmatrix,
                             matrix(z$s, MCMCsamplesize,
                                    Clist$nstats, byrow = TRUE))
      }
    }
    #
    #   Next update the network to be the final (possibly conditionally)
    #   simulated one

    
    out.list[[i]] <- newnw.extract(nw,z)
    out.mat <- rbind(out.mat,z$s[(Clist$nstats+1):(2*Clist$nstats)])
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
