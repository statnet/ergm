san <- function(object, ...){
 UseMethod("san")
}

san.default <- function(object,...)
{
  stop("Either a ergm object or a formula argument must be given")
}

san.formula <- function(object, nsim=1, seed=NULL, ...,theta0,
                        burnin=10000, interval=10000,
                        meanstats=NULL,
                        sequential=TRUE,
                        constraints=~.,
                        control=san.control(),
                        verbose=FALSE) {
  out.list <- list()
  out.mat <- numeric(0)
  formula <- object

  if(is.null(seed)){seed <- sample(10000000, size=1)}
  set.seed(as.integer(seed))
  nw <- ergm.getnetwork(formula)
  if(class(nw) =="network.series"){
    nw <- nw$networks[[1]]
  }
  nw <- as.network(nw)
  if(missing(meanstats)){
    stop("You need to specify target statistic via",
         " the 'meanstats' argument")
  }
  if(!is.network(nw)){
    stop("A network object on the LHS of the formula ",
         "must be given")
  }

  model <- ergm.getmodel(formula, nw, drop=control$drop)
  Clist <- ergm.Cprepare(nw, model)
  Clist.miss <- ergm.design(nw, model, initialfit=TRUE, verbose=verbose)
  
  MCMCsamplesize <- 1
  verb <- match(verbose,
                c("FALSE","TRUE", "very"), nomatch=1)-1
  MHproposal<-getMHproposal(constraints,control$prop.args,nw,model,weights=control$prop.weights)
  if(missing(theta0)) {
    warning("No parameter values given, using MPLE\n\t")
  }
# theta0 <- c(theta0[1],rep(0,Clist$nparam-1))
  
  if(is.null(seed)){seed <- sample(10000000, size=1)}
  set.seed(as.integer(seed))
    
  for(i in 1:nsim){
    Clist <- ergm.Cprepare(nw, model)
#   Clist.miss <- ergm.design(nw, model, initialfit=TRUE, verbose=verbose)
    maxedges <- max(20000, Clist$nedges)
    if(i==1 | !sequential){
      use.burnin <- burnin
    }else{
      use.burnin <- interval
    }

    if(missing(theta0)) {
      theta0 <- ergm.mple(Clist=Clist, Clist.miss=Clist.miss, 
                          m=model, verbose=verbose, ...)$coef
    }
    stats <- matrix(summary(model$formula)-meanstats,
                    ncol=Clist$nparam,byrow=T,nrow=MCMCsamplesize)
#
#   Check for truncation of the returned edge list
#
    z <- list(newnwheads=maxedges+1)
    while(z$newnwheads[1] > maxedges){
     maxedges <- 10*maxedges
     z <- .C("SAN_wrapper",
             as.integer(Clist$heads), as.integer(Clist$tails), 
             as.integer(Clist$nedges), as.integer(Clist$n),
             as.integer(Clist$dir), as.integer(Clist$bipartite),
             as.integer(Clist$nterms), 
             as.character(Clist$fnamestring),
             as.character(Clist$snamestring), 
             as.character(MHproposal$name),
             as.character(MHproposal$package),
             as.double(Clist$inputs),
             as.double(theta0),
             as.integer(MCMCsamplesize),
             s = as.double(stats),
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
    out.list[[i]] <- newnw.extract(nw, z)
    out.mat <- rbind(out.mat,z$s[(Clist$nparam+1):(2*Clist$nparam)])
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
                       control=san.control(),
                       verbose=FALSE) {
  out.list <- vector("list", nsim)
  out.mat <- numeric(0)
  
#  if(missing(multiplicity) & !is.null(object$multiplicity)){
#    multiplicity <- object$multiplicity
#  }
  if(is.null(seed)){seed <- sample(10000000, size=1)}
  set.seed(as.integer(seed))
  
  nw <- object$network  
  
  m <- ergm.getmodel(object$formula, nw, drop=control$drop)
  MCMCsamplesize <- 1
  verb <- match(verbose,
                c("FALSE","TRUE", "very"), nomatch=1)-1
  MHproposal<-getMHproposal(object,constraints=constraints,arguments=control$prop.args,nw=nw,model=model,weights=control$prop.weights)
  # multiplicity.constrained <- 1  
  if(missing(meanstats)){
    stop("You need to specify target statistic via",
         " the 'meanstats' argument")
  }
  if(missing(theta0)) {
    theta0 <- object$coef
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
                    ncol=Clist$nparam,byrow=T,nrow=MCMCsamplesize)
#
#   Check for truncation of the returned edge list
#
    z <- list(newnwheads=maxedges+1,newnwtails=maxedges+1)
    while(z$newnwheads[1] > maxedges){
     maxedges <- 10*maxedges
     z <- .C("SAN_wrapper",
            as.integer(Clist$heads), as.integer(Clist$tails), 
            as.integer(Clist$nedges), as.integer(Clist$n),
            as.integer(Clist$dir), as.integer(Clist$bipartite),
            as.integer(Clist$nterms), 
            as.character(Clist$fnamestring),
            as.character(Clist$snamestring), 
            as.character(MHproposal$name),
            as.character(MHproposal$package),
            as.double(Clist$inputs),
            as.double(theta0),
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

    
    out.list[[i]] <- newnw.extract(nw,z)
    out.mat <- rbind(out.mat,z$s[(Clist$nparam+1):(2*Clist$nparam)])
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
