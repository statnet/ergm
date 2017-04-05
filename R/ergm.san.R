#  File ergm/R/ergm.san.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
####################################################################
# Each of the <san.X> functions samples one or more networks via
# <SAN_wrapper.C> according to the vector of mean stats given;
# execution will halt
#    - if X is neither an ergm object or a formula
#    - if the formula does not correctly specify a network
#    - no mean stats are given
#
##############################################################################

san <- function(object, ...){
 UseMethod("san")
}

san.default <- function(object,...)
{
  stop("Either a ergm object or a formula argument must be given")
}

san.formula <- function(object, constraints=~., target.stats=NULL,
                        nsim=1, basis=NULL,
                        sequential=TRUE,
                        control=control.san(),
                        verbose=FALSE, ...) {
  out.list <- list()
  out.mat <- numeric(0)
  formula <- object

  if(!is.null(control$seed)) set.seed(as.integer(control$seed))
  if(!is.null(basis)) {
    nw <- basis
    formula <- ergm.update.formula(formula, nw ~ .)
    object <- formula
  } else {
    nw <- ergm.getnetwork(formula)
  }
  if(class(nw) =="network.list"){
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
  model <- ergm.getmodel(formula, nw)
  Clist <- ergm.Cprepare(nw, model)
  Clist.miss <- ergm.design(nw, model, verbose=verbose)
  
  verb <- match(verbose,
                c("FALSE","TRUE", "very"), nomatch=1)-1
  MHproposal<-MHproposal(constraints,arguments=control$SAN.prop.args,nw=nw,weights=control$SAN.prop.weights, class="c")
# if(is.null(control$coef)) {
#   warning("No parameter values given, using the MPLE for the passed network.\n\t")
# }
# control$coef <- c(control$coef[1],rep(0,Clist$nstats-1))
  
  if(MHproposal$name=="CondDegree"){ 
   formula.conddegmple <- ergm.update.formula(formula, ~ conddegmple + .)
   m.conddeg <- ergm.getmodel(formula.conddegmple, nw, initialfit=TRUE)
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
    Clist <- ergm.Cprepare(nw, model)
#   Clist.miss <- ergm.design(nw, model, verbose=verbose)
    maxedges <- max(control$SAN.init.maxedges, Clist$nedges)
    if (verb) {
       cat(paste("#", i, " of ", nsim, ": ", sep=""))
     }

    if(is.null(control$coef)) {
      fit <- ergm.mple(Clist=Clist, Clist.miss=Clist.miss, 
                       conddeg=conddeg,
                       control=control, MHproposal=MHproposal,
                       m=model, verbose=verbose, ...)
      control$coef <- fit$coef
      if(is.null(control$invcov)) { control$invcov <- fit$covar }
    }else{
      if(is.null(control$invcov)) { control$invcov <- diag(length(control$coef)) }
    }
    eta0 <- ergm.eta(control$coef, model$etamap)
    
    netsumm<-summary(model$formula)
    target.stats <- vector.namesmatch(target.stats, names(netsumm))
    
    stats <- matrix(netsumm-target.stats,
                    ncol=Clist$nstats,byrow=TRUE,nrow=nsim)
    tau <- rep(control$SAN.tau,length=length(eta0))
#
#   Check for truncation of the returned edge list
#
    eta0[is.na(eta0)]<-0
    repeat{
      nedges <- c(Clist$nedges,0)
      tails <- Clist$tails
      heads <- Clist$heads
      # *** don't forget to pass in tails before heads now.

        z <- .C("SAN_wrapper",
                as.integer(length(nedges)), as.integer(nedges),
                as.integer(tails), as.integer(heads),
                as.integer(Clist$n),
                as.integer(Clist$dir), as.integer(Clist$bipartite),
                as.integer(Clist$nterms), 
                as.character(Clist$fnamestring),
                as.character(Clist$snamestring), 
                as.character(MHproposal$name),
                as.character(MHproposal$package),
                as.double(c(Clist$inputs,MHproposal$inputs)),
                as.double(eta0),
                as.double(tau),
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
                PACKAGE="ergm")

      if(z$newnwtails[1] <= maxedges) break
      maxedges <- 5*maxedges
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
                     control=object$SAN.control,
                     verbose=FALSE, ...) {
  if(is.null(control$coef)) control$coef <- coef(object)
  san.formula(formula, nsim=nsim, 
              target.stats=target.stats,
              basis=basis,
              sequential=sequential,
              constraints=constraints,
              control=control,
              verbose=verbose, ...)
}

