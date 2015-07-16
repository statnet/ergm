#  File R/ergm.san.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2015 Statnet Commons
#######################################################################
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
#   init     : a vector of initial values for the theta coefficients;
#                default=those returned by <ergm.mple>
#   invcov     : the initial inverse covariance matrix used to
#                calculate the Mahalanobis distance; default=that 
#                from the mple fit if 'init'=NULL, else default=the
#                identity matrix of size 'init'
#   burnin     : the number of proposal to disregard before sampling
#                begins; default=1e4
#   interval   : the number of proposals between sampled statistics;
#                default=1e4
#   target.stats  : a vector of the mean statistics for each model
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
#            network.list object as list containing
#              formula :  the formula given by 'object'
#              networks:  the list of sampled networks
#              stats   :  the summary statistics of the sampled networks
#              coef    :  the initial theta coefficients used by
#                         the sampling rountine, i.e. 'init'
#            
##############################################################################

san <- function(object, ...){
 UseMethod("san")
}

san.default <- function(object,...)
{
  stop("Either a ergm object or a formula argument must be given")
}

san.formula <- function(object, response=NULL, reference=~Bernoulli, constraints=~., target.stats=NULL,
                        nsim=1, basis=NULL,
                        sequential=TRUE,
                        control=control.san(),
                        verbose=FALSE, ...) {
  check.control.class("san")
  
  out.list <- list()
  out.mat <- numeric(0)
  formula <- object

  if(!is.null(control$seed)) set.seed(as.integer(control$seed))
  if(!is.null(basis)) {
    nw <- basis
    formula <- ergm.update.formula(formula, nw ~ ., from.new="nw")
    object <- formula
  } else {
    nw <- ergm.getnetwork(formula)
  }
  if(inherits(nw,"network.list")){
    nw <- nw$networks[[1]]
  }
  nw <- as.network(nw)
  if(is.null(target.stats)){
    stop("You need to specify target statistic via",
         " the 'target.stats' argument")
  }
  if(!is.network(nw)){
    stop("A network object on the LHS of the formula ",
         "must be given")
  }

# model <- ergm.getmodel(formula, nw, drop=control$drop)
  model <- ergm.getmodel(formula, nw, response=response)
  Clist <- ergm.Cprepare(nw, model, response=response)
  Clist.miss <- ergm.design(nw, model, verbose=verbose)
  
  verb <- match(verbose,
                c("FALSE","TRUE", "very"), nomatch=1)-1
  MHproposal<-MHproposal(constraints,arguments=control$SAN.prop.args,nw=nw,weights=control$SAN.prop.weights, class="c",reference=reference,response=response)
# if(is.null(control$coef)) {
#   warning("No parameter values given, using the MPLE for the passed network.\n\t")
# }
# control$coef <- c(control$coef[1],rep(0,Clist$nstats-1))
  
  if(MHproposal$name %in% c("CondDegree","CondDegreeMix")){ 
   formula.conddegmple <- ergm.update.formula(formula, . ~ conddegmple + .)
   m.conddeg <- ergm.getmodel(formula.conddegmple, nw)
   Clist.conddegmple <- ergm.Cprepare(nw, m.conddeg)
   Clist.conddegmple$target.stats=c(1,target.stats)
   conddeg <- list(m=m.conddeg, Clist=Clist.conddegmple, Clist.miss=ergm.Cprepare(nw, m.conddeg))
  }else{
   conddeg <- NULL
  }

  if (verb) {
    cat(paste("Starting ",nsim," MCMC iteration", ifelse(nsim>1,"s",""),
        " of ", control$SAN.burnin+control$SAN.interval*(nsim-1), 
        " steps", ifelse(nsim>1, " each", ""), ".\n", sep=""))
  }

  for(i in 1:nsim){
    Clist <- ergm.Cprepare(nw, model,response=response)
#   Clist.miss <- ergm.design(nw, model, verbose=verbose)
    maxedges <- max(control$SAN.init.maxedges, Clist$nedges)
    if (verb) {
       cat(paste("#", i, " of ", nsim, ": ", sep=""))
     }

    if(is.null(control$coef)) {
      if(reference==~Bernoulli){
        fit <- suppressWarnings(try(ergm.mple(Clist=Clist, Clist.miss=Clist.miss, 
                         conddeg=conddeg,
                         control=control, MHproposal=MHproposal,
                         m=model, verbose=verbose, ...)))
        control$coef <- if(inherits(fit, "try-error")) rep(0,length(model$coef.names)) else fit$coef
        if(is.null(control$invcov)) { control$invcov <- fit$covar }
      }else{
        control$coef<-rep(0,length(model$coef.names))
        if(is.null(control$invcov)) control$invcov <- diag(length(control$coef))
      }
    }else{
      if(is.null(control$invcov)) { control$invcov <- diag(length(control$coef)) }
    }
    eta0 <- ifelse(is.na(control$coef), 0, control$coef)
    
    netsumm<-summary(model$formula,response=response)
    target.stats <- vector.namesmatch(target.stats, names(netsumm))
    
    stats <- matrix(netsumm-target.stats,
                    ncol=Clist$nstats,byrow=TRUE,nrow=nsim)
    tau <- rep(control$SAN.tau,length=length(eta0))
#
#   Check for truncation of the returned edge list
#
    repeat{
      nedges <- c(Clist$nedges,0)
      tails <- Clist$tails
      heads <- Clist$heads
      weights <- Clist$weights
      # *** don't forget to pass in tails before heads now.

      if(is.null(Clist$weights)){
        z <- .C("SAN_wrapper",
                as.integer(length(nedges)), as.integer(nedges),
                as.integer(tails), as.integer(heads),
                as.integer(Clist$n),
                as.integer(Clist$dir), as.integer(Clist$bipartite),
                as.integer(Clist$nterms), 
                as.character(Clist$fnamestring),
                as.character(Clist$snamestring), 
                as.character(MHproposal$name),
                as.character(MHproposal$pkgname),
                as.double(c(Clist$inputs,MHproposal$inputs)),
                as.double(.deinf(eta0)),
                as.double(.deinf(tau)),
                as.integer(1), # "samplesize"
                s = as.double(stats),
                as.integer(if(i==1 | !sequential) control$SAN.burnin else control$SAN.interval), as.integer(control$SAN.interval), 
                newnwtails = integer(maxedges),
                newnwheads = integer(maxedges), 
                as.double(control$invcov),
                as.integer(verb),
                as.integer(MHproposal$arguments$constraints$bd$attribs), 
                as.integer(MHproposal$arguments$constraints$bd$maxout), as.integer(MHproposal$arguments$constraints$bd$maxin),
                as.integer(MHproposal$arguments$constraints$bd$minout), as.integer(MHproposal$arguments$constraints$bd$minin),
                as.integer(MHproposal$arguments$constraints$bd$condAllDegExact),
                as.integer(length(MHproposal$arguments$constraints$bd$attribs)), 
                as.integer(maxedges),
                status = integer(1),
                PACKAGE="ergm")
      }else{
        z <- .C("WtSAN_wrapper",
                as.integer(length(nedges)), as.integer(nedges),
                as.integer(tails), as.integer(heads), as.double(weights),
                as.integer(Clist$n),
                as.integer(Clist$dir), as.integer(Clist$bipartite),
                as.integer(Clist$nterms), 
                as.character(Clist$fnamestring),
                as.character(Clist$snamestring), 
                as.character(MHproposal$name),
                as.character(MHproposal$pkgname),
                as.double(c(Clist$inputs,MHproposal$inputs)),
                as.double(.deinf(eta0)),
                as.double(.deinf(tau)),
                as.integer(1), # "samplesize"
                s = as.double(stats),
                as.integer(if(i==1 | !sequential) control$SAN.burnin else control$SAN.interval), as.integer(control$SAN.interval), 
                newnwtails = integer(maxedges),
                newnwheads = integer(maxedges),
                newnwweights = double(maxedges), 
                as.double(control$invcov),
                as.integer(verb),
                as.integer(MHproposal$arguments$constraints$bd$attribs), 
                as.integer(MHproposal$arguments$constraints$bd$maxout), as.integer(MHproposal$arguments$constraints$bd$maxin),
                as.integer(MHproposal$arguments$constraints$bd$minout), as.integer(MHproposal$arguments$constraints$bd$minin),
                as.integer(MHproposal$arguments$constraints$bd$condAllDegExact),
                as.integer(length(MHproposal$arguments$constraints$bd$attribs)), 
                as.integer(maxedges), 
                status = integer(1),
                PACKAGE="ergm")
      }

      if(z$status==1) maxedges <- 5*maxedges
      else if(z$status==0) break
      else stop("Error in SAN C code.")
    }
    #
    #   Next update the network to be the final (possibly conditionally)
    #   simulated one
    #
    out.list[[i]] <- newnw.extract(nw, z, output=control$network.output, response=response)
    out.mat <- rbind(out.mat,z$s[(Clist$nstats+1):(2*Clist$nstats)])
    if(sequential){
      nw <-  as.network.uncompressed(out.list[[i]])
    }
  }
  if(nsim > 1){
    out.list <- list(formula = formula, networks = out.list, 
                     stats = out.mat, coef=control$coef)
    class(out.list) <- "network.list"
  }else{
    out.list <- out.list[[1]]
  }
  return(out.list)
}

san.ergm <- function(object, formula=object$formula, 
                     constraints=object$constraints, 
                     target.stats=object$target.stats,
                     nsim=1, basis=NULL,
                     sequential=TRUE, 
                     control=object$control$SAN.control,
                     verbose=FALSE, ...) {
  if(is.null(control$coef)) control$coef <- coef(object)
  san.formula(formula, nsim=nsim, 
              target.stats=target.stats,
              basis=basis,
              reference = object$reference,
              sequential=sequential,
              constraints=constraints,
              control=control,
              verbose=verbose, ...)
}

