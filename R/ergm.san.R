#=========================================================================
# This file contains 4 functions for created "SAN-ed" networks & formulas
#           <san>              <san.formula>
#           <san.default>      <san.ergm>
#=========================================================================




####################################################################
# Each of the <san.X> functions samples one or more networks via
# <SAN_wrapper.C> according to the vector of mean stats given;
# execution will halt
#    - if X is neither an ergm object or a formula
#    - if the formula does not correctly specify a network
#    - no mean stats are given
#
# --PARAMETERS--
#   object     : an ergm object of a formula for such
#   nsim       : the number of sampled networks to return;
#                default=1
#   seed       : the number at which to start the random number
#                generator; default=NULL
#   theta0     : a vector of initial values for the theta coefficients;
#                default=those returned by <ergm.mple>
#   invcov     : the initial inverse covariance matrix used to
#                calculate the Mahalanobis distance; default=that 
#                from the mple fit if 'theta0'=NULL, else default=the
#                identity matrix of size 'theta0'
#   burnin     : the number of proposal to disregard before sampling
#                begins; default=1e4
#   interval   : the number of proposals between sampled statistics;
#                default=1e4
#   meanstats  : a vector of the mean statistics for each model
#                coefficient; default=NULL (which will halt execution)
#   basis      : optionally, a network can be provided in 'basis' and
#                this replaces that given by 'object'; default=NULL
#   sequential : whether subsequent sampling should start with the
#                previously sampled network; the alternative is to
#                always begin sampling from the original network;
#                default=TRUE
#   constraints: a one-sided formula giving one or more constraints on
#                the support of the distribution of the networks being
#                modeled; a list of availabe options is described in the
#                <ergm> R documentation; default=~.
#   control    : a control list for tuning the MHproposals and other
#                elements of the fit; default=<control.san>()
#   verbose    : whether this and the C program should be verbose;
#                default=FALSE
#   ...        : additional parameters that will passed onto <ergm.mple>
#
# --IGNORED--
#   tau:  this is passed along to several C functions; its use in
#         <SANMetropolisHastings.c> is commented out
#
# --RETURNED--
#   outlist: either a single sampled network if 'nsim'=1, else a
#            network.series object as list containing
#              formula :  the formula given by 'object'
#              networks:  the list of sampled networks
#              stats   :  the summary statistics of the sampled networks
#              coef    :  the initial theta coefficients used by
#                         the sampling rountine, i.e. 'theta0'
#            
##############################################################################

san <- function(object, ...){
 UseMethod("san")
}

san.default <- function(object,...)
{
  stop("Either a ergm object or a formula argument must be given")
}

san.formula <- function(object, nsim=1, seed=NULL, theta0=NULL,
                        tau=1, invcov=NULL,
                        burnin=10000, interval=10000,
                        meanstats=NULL,
                        basis=NULL,
                        sequential=TRUE,
                        constraints=~.,
                        control=control.san(),
                        verbose=FALSE, ...) {
  out.list <- list()
  out.mat <- numeric(0)
  formula <- object

  if(!is.null(seed)) set.seed(as.integer(seed))
  if(!is.null(basis)) {
    nw <- basis
    formula <- ergm.update.formula(formula, nw ~ .)
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

  if(MHproposal$name=="CondDegree"){ 
   formula.conddegmple <- ergm.update.formula(formula, ~ conddegmple + .)
   m.conddeg <- ergm.getmodel(formula.conddegmple, nw, drop=FALSE, initialfit=TRUE)
   Clist.conddegmple <- ergm.Cprepare(nw, m.conddeg)
   Clist.conddegmple$meanstats=c(1,meanstats)
   conddeg <- list(m=m.conddeg, Clist=Clist.conddegmple, Clist.miss=ergm.Cprepare(nw, m.conddeg))
  }else{
   conddeg <- NULL
  }

  if (verb) {
    cat(paste("Starting ",nsim," MCMC iteration", ifelse(nsim>1,"s",""),
        " of ", burnin+interval*(MCMCsamplesize-1), 
        " steps", ifelse(nsim>1, " each", ""), ".\n", sep=""))
  }

  MCMCparams=c(control,
     list(samplesize=burnin, burnin=burnin, interval=interval,
          Clist.miss=Clist.miss))

  for(i in 1:nsim){
    Clist <- ergm.Cprepare(nw, model)
#   Clist.miss <- ergm.design(nw, model, verbose=verbose)
    maxedges <- max(control$maxedges, Clist$nedges)
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
                       conddeg=conddeg,
                       MCMCparams=MCMCparams, MHproposal=MHproposal,
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
    z <- list(newnwtails=maxedges+1)
    eta0[is.na(eta0)]<-0
    while(z$newnwtails[1] > maxedges){
     maxedges <- 10*maxedges
     nedges <- c(Clist$nedges,0)
     tails <- Clist$tails
     heads <- Clist$heads
     if(!is.null(MCMCparams$Clist.miss)){
       nedges[2] <- MCMCparams$Clist.miss$nedges
       tails <- c(tails, MCMCparams$Clist.miss$tails)
       heads <- c(heads, MCMCparams$Clist.miss$heads)
     }
     # *** don't forget to pass in tails before heads now.
     z <- .C("SAN_wrapper",
             as.integer(length(nedges)), as.integer(nedges),
             as.integer(tails), as.integer(heads),
             as.integer(Clist$maxpossibleedges),
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
             newnwtails = integer(maxedges),
             newnwheads = integer(maxedges), 
             as.double(invcov),
             as.integer(verb),
             as.integer(MHproposal$bd$attribs), 
             as.integer(MHproposal$bd$maxout), as.integer(MHproposal$bd$maxin),
             as.integer(MHproposal$bd$minout), as.integer(MHproposal$bd$minin),
             as.integer(MHproposal$bd$condAllDegExact),
             as.integer(length(MHproposal$bd$attribs)), 
             as.integer(maxedges), 
             PACKAGE="ergm")
    }
#
#   Next update the network to be the final (possibly conditionally)
#   simulated one
#
    out.list[[i]] <- newnw.extract(nw, z, output=control$network.output)
    out.mat <- rbind(out.mat,z$s[(Clist$nstats+1):(2*Clist$nstats)])
    if(sequential){
      nw <-  as.network.uncompressed(out.list[[i]])
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


san.ergm <- function(object, nsim=1, seed=NULL, theta0=object$coef,
                       burnin=10000, interval=10000, 
                       meanstats=NULL,
                       basis=NULL,
                       sequential=TRUE, 
                       constraints=NULL,
                       control=control.san(),
                       verbose=FALSE, ...) {
  san.formula(object$formula, nsim=nsim, seed=seed, theta0=theta0,
              burnin=burnin, interval=interval,
              meanstats=meanstats,
              basis=basis,
              sequential=sequential,
              constraints=constraints,
              control=control,
              verbose=verbose, ...)
}

