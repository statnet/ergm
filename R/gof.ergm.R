#  File R/gof.ergm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
#=============================================================================
# This file contains the following 8 functions for assessing goodness of fit
#         <gof>              <summary.gofobject>
#         <gof.default>      <plot.gofobject>
#         <gof.ergm>         <ergm.get.terms.formula>
#         <gof.formula>      <ergm.rhs.formula>
#=============================================================================



###############################################################################
# Each of the <gof.X> functions assesses the goodness of fit of X by comparison
# with 'control$nsim' ergm simulations of X
#
# --PARAMETERS--
#   object/formula: either an ergm object or a formula
#   ...           : additional parameters passed from within the program;
#                   these are ignored
#   init        : the parameters from which the simulations will be drawn;
#                   default=NULL;
#   control$nsim          : the number of simulated ergms, with which to compare X;
#                   default=100
#   burnin        : the number of proposals to disregard before any MCMC
#                   sampling is done; this is passed along to the simulation
#                   routines; default=10000
#   interval      : the number of proposals between sampled ergm statistics;
#                   this is passed along to the simulation rountines;
#                   default=1000
#   GOF           : a one-sided formula specifying which summary statistics
#                   should be used in the GOF comparison; choices include
#                       distance      espartners    dspartners
#                       odegree       idegree       degree
#                       triadcensus   model
#                   default=NULL; is internally mapped to 
#                   ~degree+espartners+distance if nw is undirected, and
#                   ~idegree+odegree+espartners+distance otherwise
#   constraints   : a one-sided formula of the constraint terms; options are
#                         bd        degrees        nodegrees
#                         edges     degreedist     idegreedist
#                         observed  odegreedist
#                   default="~ ."   
#   control       : a list of parameters for controlling GOF evaluation, as
#                   returned by <control.gof.X>; default=control.gof.X()
#                   (note that <control.gof.X> has different defaults 
#                    depending on the class of X)
#   seed          : an integer value at which to set the random generator;
#                   default=NULL
#   verbose       : whether to print information on the progress of the
#                   simulations; default=FALSE
#
# --RETURNED--
#   returnlist: a list with the following components for each term
#               G given in 'GOF'
#      summary.G: a matrix of summary statistics for the observed and
#                 simulated G's; if G takes on the values {G1, G2,...,Gq},
#                 the entries of 'summary.G' are
#         [i,1]-- the observed frequency of Gi
#         [i,2]-- the minimum value of Gi from the simulations
#         [i,3]-- the mean value of Gi from the simulations
#         [i,4]-- the maximum value of Gi from the simulations
#         [i,5]-- the p-value for the observed Gi estimated from the
#                 distribution of simulations
#      pobs.G   : a vector giving G's observed probability distribution
#      psim.G   : a matrix of G's simulated probability distributions; each
#                 row gives a distribution
#      bds.G    : the estimatd confidence interval, as the .025 and .975
#                 quantile values of the simulations
#      obs.G    : the vector of summary statistics for the observed X
#      sim.G    : the matrix of summary statistics for each simulated
#                 version of X
#
###############################################################################

gof <- function(object, ...){
      UseMethod("gof")
    }


gof.default <- function(object,...) {
  classes <- setdiff(gsub(pattern="^gof.",replacement="",as.vector(methods("gof"))), "default")
  stop("Goodness-of-Fit methods have been implemented only for class(es) ",
       paste.and(paste('"',classes,'"',sep="")), " in the packages loaded.")
}


gof.ergm <- function (object, ..., 
                      coef=NULL,
                      GOF=NULL, 
                      constraints=NULL,
                      control=control.gof.ergm(),
                      verbose=FALSE) {
  check.control.class(c("gof.ergm","gof.formula"), "gof.ergm")
  control.toplevel(...)
  nw <- as.network(object$network)

  if(!is.null(object$response)) stop("GoF for valued ERGMs is not implemented at this time.")
  
  #Set up the defaults, if called with GOF==NULL
  if(is.null(GOF)){
    if(is.directed(nw))
      GOF<- ~idegree + odegree + espartners + distance + model
    else
      GOF<- ~degree + espartners + distance + model
  }
  # Add a model term, unless it is explicitly excluded
  model_trms <- unlist(dimnames(attr(terms(GOF),"factors"))[1])
  if(!("model" %in% model_trms)){
    GOF <- nonsimp.update.formula(GOF, ~ . + model)
  }

  ## FIXME: Need to do this differently. This approach will (probably)
  ## break if any of the terms in the formula have non-constant
  ## arguments.
  ## Also, this is just plain ugly.
  formula <- as.formula(paste("nw ~",paste(ergm.rhs.formula(object$formula),collapse="+")))
# paste("~",paste(unlist(dimnames(attr(terms(formula),"factors"))[-1]),collapse="+"),sep="")
  if(!is.network(nw)){
    stop("A network must be given as part of the network object.")
  }

  if(missing(coef)){coef <- object$coef}

  ## If a different constraint was specified, use it; otherwise, copy
  ## from the ERGM.

  control.transfer <- c("MCMC.burnin", "MCMC.interval", "MCMC.prop.weights", "MCMC.prop.args", "MCMC.packagenames", "MCMC.init.maxedges")
  for(arg in control.transfer)
    if(is.null(control[[arg]]))
      control[arg] <- list(object$control[[arg]])

  control <- set.control.class("control.gof.formula")
  
  if(is.null(constraints)) constraints <- object$constraints
  
  gof.formula(object=formula, coef=coef,
              GOF=GOF,
              constraints=constraints,
              control=control,
              verbose=verbose, ...)
}



gof.formula <- function(object, ..., 
                        coef=NULL,
                        GOF=NULL,
                        constraints=~.,
                        control=control.gof.formula(),
			unconditional=TRUE,
                        verbose=FALSE) {
  check.control.class(c("gof.formula","gof.ergm"), "ERGM gof.formula")
  control.toplevel(...)

  if("response" %in% names(list(...))) stop("GoF for valued ERGMs is not implemented at this time.")

  if(!is.null(control$seed)) {set.seed(as.integer(control$seed))}
  if (verbose) 
    message("Starting GOF for the given ERGM formula.")
  # Unused code
  coefmissing <- NULL
  # get network
  trms <- ergm.getterms(object)
  if(length(trms)>2){
    nw <- eval(trms[[2]], sys.parent())
  }else{
    stop("A network object on the RHS of the formula argument must be given")
  }
  if(is.ergm(nw)){
    all.gof.vars <- ergm.rhs.formula(object)
    object <- nw$formula
    if(missing(coef)){coef <- nw$coef}
    trms <- ergm.getterms(object)
    if(length(trms)>2){
      nw <- eval(trms[[2]], sys.parent())
    }else{
      stop("A network object on the RHS of the formula argument must be given")
    }
  }else{
    nw <- as.network(nw)
    if(!is.network(nw)){
      stop("A network object on the RHS of the formula argument must be given")
    }
    if(is.null(GOF)){
      if(is.directed(nw))
        GOF<- ~idegree + odegree + espartners + distance + model
      else
        GOF<- ~degree + espartners + distance + model
    }

    # Add a model term, unless it is explicitly excluded
    model_trms <- unlist(dimnames(attr(terms(GOF),"factors"))[1])
    if(!("model" %in% model_trms)){
      GOF <- nonsimp.update.formula(GOF, ~ . + model)
    }
  
    all.gof.vars <- ergm.rhs.formula(GOF)
  }

# match variables

  for(i in seq(along=all.gof.vars)){
    all.gof.vars[i] <- match.arg(all.gof.vars[i],
                                 c('distance', 'espartners', 'dspartners', 'odegree', 'idegree', 
                                   'degree', 'triadcensus', 'model'
                                   )
                                 )
  }
  GOF <- as.formula(paste("~",paste(all.gof.vars,collapse="+")))
  
  nw <- as.network(nw)
  
  if(!is.network(nw)){
    stop("A network object on the RHS of the formula argument must be given")
  }

# if(is.bipartite(nw)){
#   object <- ergm.update.formula(object, ~ . + bipartite)
#   trms <- ergm.getterms(object)
#   termnames <- ergm.gettermnames(trms)
# }

  m <- ergm.getmodel(object, nw)
  Clist <- ergm.Cprepare(nw, m)

  if(is.null(coef)){
      coef <- rep(0,Clist$nstats)
      warning("No parameter values given, using 0\n\t")
  }
# if(is.bipartite(nw)){
#     coef <- c(coef,-1)
# }

  # If missing simulate from the conditional model
  if(network.naedgecount(nw) & unconditional){
   if(verbose){message("Conditional simulations for missing fit")}
   if(is.null(coefmissing)){coefmissing <- coef}
   constraints.obs<-ergm.update.formula(constraints,~.+observed)
   SimCond <- gof(object=object, coef=coefmissing,
                  GOF=GOF, 
                  constraints=constraints.obs,
                  control=control,
                  unconditional=FALSE,
                  verbose=verbose)
  }

# test to see which of these is/are necessary
#  pval.model<-pval.triadcensus<-pval.dist<-pval.deg<-pval.espart<-pval.espart<-NULL
##
#  obs.model<-pobs.model<-sim.model<-psim.model<-pval.model<-bds.model<-NULL
#  obs.triadcensus<-pobs.triadcensus<-sim.triadcensus<-psim.triadcensus<-pval.triadcensus<-bds.triadcensus<-NULL
#  obs.dist<-pobs.dist<-sim.dist<-psim.dist<-pval.dist<-bds.dist<-NULL
#  obs.deg<-pobs.deg<-sim.deg<-psim.deg<-pval.deg<-bds.deg<-NULL
#  obs.espart<-pobs.espart<-sim.espart<-psim.espart<-pval.espart<-bds.espart<-NULL
#  obs.dspart<-pobs.dspart<-sim.dspart<-psim.dspart<-pval.dspart<-bds.dspart<-NULL
#
#  obs.ideg<-pobs.ideg<-sim.ideg<-psim.ideg<-pval.ideg<-bds.ideg<-pval.ideg<-NULL
#  obs.odeg<-pobs.odeg<-sim.odeg<-psim.odeg<-pval.odeg<-bds.odeg<-pval.odeg<-NULL

  n <- network.size(nw)

  # Calculate network statistics for the observed graph
  # Set up the output arrays of sim variables
  if(verbose)
    message("Calculating observed network statistics.")
  
  if ('model' %in% all.gof.vars) {
   if(!network.naedgecount(nw) | !unconditional){
    obs.model <- summary(object)
   }else{
    obs.model <- SimCond$obs.model
   }
   sim.model <- array(0,dim=c(control$nsim,length(obs.model)))
   dimnames(sim.model) <- list(paste(c(1:control$nsim)),names(obs.model))
  }

  if ('distance' %in% all.gof.vars) {
   if(!network.naedgecount(nw) | !unconditional){
    obs.dist <- ergm.geodistdist(nw)
    obs.dist[obs.dist==Inf] <- n
   }else{
    obs.dist <- SimCond$summary.dist[,"mean"]
   }
   sim.dist <-array(0,dim=c(control$nsim,n))
   dimnames(sim.dist)  <- list(paste(c(1:control$nsim)),paste(1:n))
  }

  if ('odegree' %in% all.gof.vars) {
   if(!network.naedgecount(nw) | !unconditional){
    mesp <- paste("c(",paste(0:(n-1),collapse=","),")",sep="")
    obs.odeg <- summary(as.formula(paste('nw ~ odegree(',mesp,')',sep="")))
   }else{
    obs.odeg <- SimCond$summary.odeg[,"mean"]
   }
   sim.odeg <- array(0,dim=c(control$nsim,n))
#  obs.odeg <- c(obs.odeg,rep(0,n-length(obs.odeg)))
   dimnames(sim.odeg)   <- list(paste(c(1:control$nsim)),paste(0:(n-1)))
   names(obs.odeg) <- dimnames(sim.odeg)[[2]]
  }

  if ('idegree' %in% all.gof.vars) {
   if(!network.naedgecount(nw) | !unconditional){
    mesp <- paste("c(",paste(0:(n-1),collapse=","),")",sep="")
    obs.ideg <- summary(as.formula(paste('nw ~ idegree(',mesp,')',sep="")))
   }else{
    obs.ideg <- SimCond$summary.ideg[,"mean"]
   }
   sim.ideg <- array(0,dim=c(control$nsim,n))
#  obs.ideg <- c(obs.ideg,rep(0,n-length(obs.ideg)))
   dimnames(sim.ideg)   <- list(paste(c(1:control$nsim)),paste(0:(n-1)))
   names(obs.ideg) <- dimnames(sim.ideg)[[2]]
  }

  if ('degree' %in% all.gof.vars) {
   if(!network.naedgecount(nw) | !unconditional){
    if(is.bipartite(nw)){
     obs.deg <- degreedist(nw, print=FALSE)$b2
     obs.deg <- c(obs.deg,rep(0,n-length(obs.deg)))
    }else{
     mesp <- paste("c(",paste(0:(n-1),collapse=","),")",sep="")
     obs.deg <- summary(as.formula(paste('nw ~ degree(',mesp,')',sep="")))
    }
   }else{
    obs.deg <- SimCond$summary.deg[,"mean"]
   }
   sim.deg <- array(0,dim=c(control$nsim,n))
   dimnames(sim.deg)   <- list(paste(c(1:control$nsim)),paste(0:(n-1)))
   names(obs.deg) <- dimnames(sim.deg)[[2]]
  }
 
  if ('espartners' %in% all.gof.vars) {
#  obs.espart <- espartnerdist(nw, print=verbose)
   if(!network.naedgecount(nw) | !unconditional){
    mesp <- paste("c(",paste(0:(network.size(nw)-2),collapse=","),")",sep="")
    obs.espart <- summary(as.formula(paste('nw ~ esp(',mesp,')',sep="")))
   }else{
    obs.espart <- SimCond$summary.espart[,"mean"]
   }
   sim.espart <- array(0,dim=c(control$nsim,n-1))
   dimnames(sim.espart) <- list(paste(c(1:control$nsim)),paste(0:(n-2)))
  }
 
  if ('dspartners' %in% all.gof.vars) {
   if(!network.naedgecount(nw) | !unconditional){
#   obs.dspart <- dspartnerdist(nw, print=verbose)
    mesp <- paste("c(",paste(0:(network.size(nw)-2),collapse=","),")",sep="")
    obs.dspart <- summary(as.formula(paste('nw ~ dsp(',mesp,')',sep="")))
   }else{
    obs.dspart <- SimCond$summary.dspart[,"mean"]
   }
   sim.dspart <- array(0,dim=c(control$nsim,n-1))
   dimnames(sim.dspart) <- list(paste(c(1:control$nsim)),paste(0:(n-2)))
  }

  if ('triadcensus' %in% all.gof.vars) {
   if(is.directed(nw)){
    triadcensus <- 0:15
    namestriadcensus <- c("003","012", "102", "021D", "021U", "021C",
      "111D", "111U", "030T",
      "030C", "201", "120D", "120U", "120C", "210", "300")
    triadcensus.formula <- "~ triadcensus(0:15)"
   }else{
    triadcensus <- 0:3
    namestriadcensus <- c("0","1","2", "3")
    triadcensus.formula <- "~ triadcensus(0:3)"
   }
   if(!network.naedgecount(nw) | !unconditional){
    obs.triadcensus <- summary(as.formula(paste('nw',triadcensus.formula,sep="")))
   }else{
    obs.triadcensus <- SimCond$summary.triadcensus[,"mean"]
   }
   sim.triadcensus <- array(0,dim=c(control$nsim,length(triadcensus)))
   dimnames(sim.triadcensus) <- list(paste(c(1:control$nsim)), namestriadcensus)
   names(obs.triadcensus) <- namestriadcensus
  }
 
  # Simulate an exponential family random graph model

#  SimNetworkSeriesObj <- simulate(object, control$nsim=control$nsim, seed=seed,
#                                  coef=coef,
#                                  burnin=burnin, interval=interval,
#                                  constraints=constraints,
#                                  control=control.simulate.formula(
#                                   prop.args=control$MCMC.prop.args,
#                                   prop.weights=control$MCMC.prop.weights,
#                                   summarizestats=control$summarizestats,
#                                   drop=control$drop),
#                                  verbose=verbose, basis=nw)
# New approach below avoids having to store gigantic unnecessary
# network.list object

  if(verbose)
    message("Starting simulations.")

  tempnet <- nw
  for (i in 1:control$nsim) {
    if(verbose){
      message("Sim ",i," of ",control$nsim,": ",appendLF=FALSE)
    }
    if(network.naedgecount(nw) & !unconditional){tempnet <- nw}
    tempnet <- simulate(object, nsim=1, coef=coef,
                        constraints=constraints, 
                        control=set.control.class("control.simulate.formula",control),
                        basis=tempnet,
                        verbose=verbose)
    seed <- NULL # Don't re-seed after first iteration   
#    if(verbose){
#     cat(paste("...",i,sep=""))
#     if ((i %% 10 == 0) || (i==control$nsim)) cat("\n")
#    }
    if ('model' %in% all.gof.vars) {
     sim.model[i,] <- summary(ergm.update.formula(object,tempnet ~ ., from.new="tempnet"))
    }
    if ('distance' %in% all.gof.vars) {
     sim.dist[i,] <- ergm.geodistdist(tempnet)
    }
    if ('idegree' %in% all.gof.vars) {
     mesp <- paste("c(",paste(0:(n-1),collapse=","),")",sep="")
     gi <- tempnet
     sim.ideg[i,] <- summary(as.formula(paste('gi ~ idegree(',mesp,')',sep="")))
#    temp <- table(degreedist(tempnet, print=verbose)[1,])
#    sim.ideg[i,] <- c(temp, rep(0, n-length(temp)))
    }
    if ('odegree' %in% all.gof.vars) {
     mesp <- paste("c(",paste(0:(n-1),collapse=","),")",sep="")
     gi <- tempnet
     sim.odeg[i,] <- summary(as.formula(paste('gi ~ odegree(',mesp,')',sep="")))
#    temp <- table(degreedist(tempnet, print=verbose)[2,])
#    sim.odeg[i,] <- c(temp, rep(0, n-length(temp)))
    }
    if ('degree' %in% all.gof.vars) {
     gi <- tempnet
     if(is.bipartite(gi)){
      temp <- degreedist(gi, print=FALSE)$b2
      sim.deg[i,] <- c(temp,rep(0,n-length(temp)))
     }else{                                                
      mesp <- paste("c(",paste(0:(n-1),collapse=","),")",sep="")
      sim.deg[i,] <- summary(as.formula(paste('gi ~ degree(',mesp,')',sep="")))
     }
#    temp <- table(degreedist(tempnet, print=verbose))
#    sim.deg[i,] <- c(temp, rep(0, n-length(temp)))
    }
    if ('espartners' %in% all.gof.vars) {
#    sim.espart[i,] <- espartnerdist(tempnet,
#                                   print=verbose)
     gi <- tempnet
     mesp <- paste("c(",paste(0:(network.size(gi)-2),collapse=","),")",sep="")
     sim.espart[i,] <- summary(as.formula(paste('gi ~ esp(',mesp,')',sep="")))
    }
    if ('dspartners' %in% all.gof.vars) {
#    sim.espart[i,] <- dspartnerdist(tempnet,
#                                   print=verbose)
     gi <- tempnet
     mesp <- paste("c(",paste(0:(network.size(gi)-2),collapse=","),")",sep="")
     sim.dspart[i,] <- summary(as.formula(paste('gi ~ dsp(',mesp,')',sep="")))
    }
    if ('triadcensus' %in% all.gof.vars) {
     gi <- tempnet
     sim.triadcensus[i,] <- summary(as.formula(paste('gi',triadcensus.formula,sep="")))
    }
  }
  if(verbose){
    message("")
  }

  # calculate p-values
  
  returnlist <- list(network.size=n, GOF=GOF)
  
  if ('model' %in% all.gof.vars) {
    pval.model <- apply(sim.model <= obs.model[col(sim.model)],2,mean)
    pval.model.top <- apply(sim.model >= obs.model[col(sim.model)],2,mean)
    pval.model <- cbind(obs.model,apply(sim.model, 2,min), apply(sim.model, 2,mean),
                        apply(sim.model, 2,max), pmin(1,2*pmin(pval.model,pval.model.top)))
    dimnames(pval.model)[[2]] <- c("obs","min","mean","max","MC p-value")
    pobs.model <- pval.model.top
    psim.model <- apply(sim.model,2,rank)/nrow(sim.model)
    psim.model <- matrix(psim.model, ncol=ncol(sim.model)) # Guard against the case of sim.model having only one row.
    bds.model <- apply(psim.model,2,quantile,probs=c(0.025,0.975))

    returnlist$summary.model <- returnlist$pval.model <- pval.model
    returnlist$pobs.model <- pobs.model
    returnlist$psim.model <- psim.model
    returnlist$bds.model <- bds.model
    returnlist$obs.model <- obs.model
    returnlist$sim.model <- sim.model
  }

  if ('distance' %in% all.gof.vars) {
    pval.dist <- apply(sim.dist <= obs.dist[col(sim.dist)],2,mean)
    pval.dist.top <- apply(sim.dist >= obs.dist[col(sim.dist)],2,mean)
    pval.dist <- cbind(obs.dist,apply(sim.dist, 2,min), apply(sim.dist, 2,mean),
                       apply(sim.dist, 2,max), pmin(1,2*pmin(pval.dist,pval.dist.top)))
    dimnames(pval.dist)[[2]] <- c("obs","min","mean","max","MC p-value")
    pobs.dist <- obs.dist/sum(obs.dist)
    psim.dist <- sweep(sim.dist,1,apply(sim.dist,1,sum),"/")
    psim.dist[is.na(psim.dist)] <- 1
    bds.dist <- apply(psim.dist,2,quantile,probs=c(0.025,0.975))
    
    returnlist$summary.dist <- returnlist$pval.dist <- pval.dist
    returnlist$pobs.dist <- pobs.dist
    returnlist$psim.dist <- psim.dist
    returnlist$bds.dist <- bds.dist
    returnlist$obs.dist <- obs.dist
    returnlist$sim.dist <- sim.dist
  }

  if ('idegree' %in% all.gof.vars) {
    pval.ideg <- apply(sim.ideg <= obs.ideg[col(sim.ideg)],2,mean)
    pval.ideg.top <- apply(sim.ideg >= obs.ideg[col(sim.ideg)],2,mean)
    pval.ideg <- cbind(obs.ideg,apply(sim.ideg, 2,min), apply(sim.ideg, 2,mean),
                       apply(sim.ideg, 2,max), pmin(1,2*pmin(pval.ideg,pval.ideg.top)))
    dimnames(pval.ideg)[[2]] <- c("obs","min","mean","max","MC p-value")
    pobs.ideg <- obs.ideg/sum(obs.ideg)
    psim.ideg <- sweep(sim.ideg,1,apply(sim.ideg,1,sum),"/")
    psim.ideg[is.na(psim.ideg)] <- 1
    bds.ideg <- apply(psim.ideg,2,quantile,probs=c(0.025,0.975))
    
    returnlist$summary.ideg <- returnlist$pval.ideg <- pval.ideg
    returnlist$pobs.ideg <- pobs.ideg
    returnlist$psim.ideg <- psim.ideg
    returnlist$bds.ideg <- bds.ideg
    returnlist$obs.ideg <- obs.ideg
    returnlist$sim.ideg <- sim.ideg
  }

  if ('odegree' %in% all.gof.vars) {
    pval.odeg <- apply(sim.odeg <= obs.odeg[col(sim.odeg)],2,mean)
    pval.odeg.top <- apply(sim.odeg >= obs.odeg[col(sim.odeg)],2,mean)
    pval.odeg <- cbind(obs.odeg,apply(sim.odeg, 2,min), apply(sim.odeg, 2,mean),
                       apply(sim.odeg, 2,max), pmin(1,2*pmin(pval.odeg,pval.odeg.top)))
    dimnames(pval.odeg)[[2]] <- c("obs","min","mean","max","MC p-value")
    pobs.odeg <- obs.odeg/sum(obs.odeg)
    psim.odeg <- sweep(sim.odeg,1,apply(sim.odeg,1,sum),"/")
    psim.odeg[is.na(psim.odeg)] <- 1
    bds.odeg <- apply(psim.odeg,2,quantile,probs=c(0.025,0.975))
    
    returnlist$summary.odeg <- returnlist$pval.odeg <- pval.odeg
    returnlist$pobs.odeg <- pobs.odeg
    returnlist$psim.odeg <- psim.odeg
    returnlist$bds.odeg <- bds.odeg
    returnlist$obs.odeg <- obs.odeg
    returnlist$sim.odeg <- sim.odeg
  }

  if ('degree' %in% all.gof.vars) {
    pval.deg <- apply(sim.deg <= obs.deg[col(sim.deg)],2,mean)
    pval.deg.top <- apply(sim.deg >= obs.deg[col(sim.deg)],2,mean)
    pval.deg <- cbind(obs.deg,apply(sim.deg, 2,min), apply(sim.deg, 2,mean),
                      apply(sim.deg, 2,max), pmin(1,2*pmin(pval.deg,pval.deg.top)))
    dimnames(pval.deg)[[2]] <- c("obs","min","mean","max","MC p-value")
    pobs.deg <- obs.deg/sum(obs.deg)
    psim.deg <- sweep(sim.deg,1,apply(sim.deg,1,sum),"/")
    psim.deg[is.na(psim.deg)] <- 1
    bds.deg <- apply(psim.deg,2,quantile,probs=c(0.025,0.975))
    
    returnlist$summary.deg <- returnlist$pval.deg <- pval.deg
    returnlist$pobs.deg <- pobs.deg
    returnlist$psim.deg <- psim.deg
    returnlist$bds.deg <- bds.deg
    returnlist$obs.deg <- obs.deg
    returnlist$sim.deg <- sim.deg
  }

  if ('espartners' %in% all.gof.vars) {
    pval.espart <- apply(sim.espart <= obs.espart[col(sim.espart)],2,mean)
    pval.espart.top <- apply(sim.espart >= obs.espart[col(sim.espart)],2,mean)
    pval.espart <- cbind(obs.espart,apply(sim.espart, 2,min), apply(sim.espart, 2,mean),
                         apply(sim.espart, 2,max), pmin(1,2*pmin(pval.espart,pval.espart.top)))
    dimnames(pval.espart)[[2]] <- c("obs","min","mean","max","MC p-value")
    pobs.espart <- obs.espart/sum(obs.espart)
    psim.espart <- sweep(sim.espart,1,apply(sim.espart,1,sum),"/")
    psim.espart[is.na(psim.espart)] <- 1
    bds.espart <- apply(psim.espart,2,quantile,probs=c(0.025,0.975))
    
    returnlist$summary.espart <- returnlist$pval.espart <- pval.espart
    returnlist$pobs.espart <- pobs.espart
    returnlist$psim.espart <- psim.espart
    returnlist$bds.espart <- bds.espart
    returnlist$obs.espart <- obs.espart
    returnlist$sim.espart <- sim.espart
  }

  if ('dspartners' %in% all.gof.vars) {
    pval.dspart <- apply(sim.dspart <= obs.dspart[col(sim.dspart)],2,mean)
    pval.dspart.top <- apply(sim.dspart >= obs.dspart[col(sim.dspart)],2,mean)
    pval.dspart <- cbind(obs.dspart,apply(sim.dspart, 2,min), apply(sim.dspart, 2,mean),
                         apply(sim.dspart, 2,max), pmin(1,2*pmin(pval.dspart,pval.dspart.top)))
    dimnames(pval.dspart)[[2]] <- c("obs","min","mean","max","MC p-value")
    pobs.dspart <- obs.dspart/sum(obs.dspart)
    psim.dspart <- sweep(sim.dspart,1,apply(sim.dspart,1,sum),"/")
    psim.dspart[is.na(psim.dspart)] <- 1
    bds.dspart <- apply(psim.dspart,2,quantile,probs=c(0.025,0.975))
    
    returnlist$summary.dspart <- returnlist$pval.dspart <- pval.dspart
    returnlist$pobs.dspart <- pobs.dspart
    returnlist$psim.dspart <- psim.dspart
    returnlist$bds.dspart <- bds.dspart
    returnlist$obs.dspart <- obs.dspart
    returnlist$sim.dspart <- sim.dspart
  }

  if ('triadcensus' %in% all.gof.vars) {
    pval.triadcensus <- apply(sim.triadcensus <= obs.triadcensus[col(sim.triadcensus)],2,mean)
    pval.triadcensus.top <- apply(sim.triadcensus >= obs.triadcensus[col(sim.triadcensus)],2,mean)
    pval.triadcensus <- cbind(obs.triadcensus,apply(sim.triadcensus, 2,min), apply(sim.triadcensus, 2,mean),
                              apply(sim.triadcensus, 2,max), pmin(1,2*pmin(pval.triadcensus,pval.triadcensus.top)))
    dimnames(pval.triadcensus)[[2]] <- c("obs","min","mean","max","MC p-value")
    pobs.triadcensus <- obs.triadcensus/sum(obs.triadcensus)
    psim.triadcensus <- sweep(sim.triadcensus,1,apply(sim.triadcensus,1,sum),"/")
    psim.triadcensus[is.na(psim.triadcensus)] <- 1
    bds.triadcensus <- apply(psim.triadcensus,2,quantile,probs=c(0.025,0.975))
    
    returnlist$summary.triadcensus <- returnlist$pval.triadcensus <- pval.triadcensus
    returnlist$pobs.triadcensus <- pobs.triadcensus
    returnlist$psim.triadcensus <- psim.triadcensus
    returnlist$bds.triadcensus <- bds.triadcensus
    returnlist$obs.triadcensus <- obs.triadcensus
    returnlist$sim.triadcensus <- sim.triadcensus
  }
  class(returnlist) <- "gofobject"
  returnlist
}



################################################################
# The <print.gofobject> function prints the summary matrices
# of each GOF term included in the build of the gofobject
#
# --PARAMETERS--
#   x  : a gofobject, as returned by one of the <gof.X> functions
#   ...: additional printing parameters; these are ignored
#
# --RETURNED--
#   NULL
#################################################################

print.gofobject <- function(x, ...){
  all.gof.vars <- ergm.rhs.formula(x$GOF)
  # match variables
  goftypes <- matrix( c(
      "model", "model statistics", "summary.model",
      "distance", "minimum geodesic distance", "summary.dist",
      "idegree", "in-degree", "summary.ideg",
      "odegree", "out-degree", "summary.odeg",
      "degree", "degree", "summary.deg",
      "espartners", "edgewise shared partner", "summary.espart",
      "dspartners", "dyadwise shared partner", "summary.dspart",
      "triadcensus", "triad census", "summary.triadcensus"), 
                      byrow=TRUE, ncol=3)
  for(i in seq(along=all.gof.vars)){
    all.gof.vars[i] <- match.arg(all.gof.vars[i], goftypes[,1])
  }
  for(statname in all.gof.vars){
    r <- match(statname, goftypes[,1])  # find row in goftypes matrix
    cat("\nGoodness-of-fit for", goftypes[r, 2],"\n\n")
    m <- x[[goftypes[r, 3] ]] # get summary statistics
    zerorows <- m[,"obs"]==0 & m[,"min"]==0 & m[,"max"]==0
    print(m[!zerorows,])
  }
  invisible()
}



summary.gofobject <- function(object, ...) {
  print.gofobject(object, ...) # Nothing better for now
}


###################################################################
# The <plot.gofobject> function plots the GOF diagnostics for each
# term included in the build of the gofobject
#
# --PARAMETERS--
#   x          : a gofobject, as returned by one of the <gof.X>
#                functions
#   ...        : additional par arguments to send to the native R
#                plotting functions
#   cex.axis   : the magnification of the text used in axis notation;
#                default=0.7
#   plotlogodds: whether the summary results should be presented
#                as their logodds; default=FALSE
#   main       : the main title; default="Goodness-of-fit diagnostics"
#   verbose    : this parameter is ignored; default=FALSE
#   normalize.reachibility: whether to normalize the distances in
#                the 'distance' GOF summary; default=FALSE
#
# --RETURNED--
#   NULL
#
###################################################################

plot.gofobject <- function(x, ..., 
         cex.axis=0.7, plotlogodds=FALSE,
         main="Goodness-of-fit diagnostics", 
         normalize.reachability=FALSE,
         verbose=FALSE) {

 color <- "gray75"
#par(oma=c(0.5,2,1,0.5))

#statsno <- (sum(stats=='deg')>0) + (sum(stats=='espart')>0) + (sum(stats=='d
 all.gof.vars <- ergm.rhs.formula(x$GOF)
 statsno <- length(all.gof.vars)

# match variables

 for(i in seq(along=all.gof.vars)){
   all.gof.vars[i] <- match.arg(all.gof.vars[i],
    c('distance', 'triadcensus', 'espartners', 'dspartners', 'odegree', 'idegree', 
      'degree', 'model'
     )
                               )
 }
 GOF <- as.formula(paste("~",paste(all.gof.vars,collapse="+")))

 if(statsno==0){
  stop("The gof object does not contain any statistics!\n")
 }
 n <- x$network.size

#attach(x)
  
 ###model####

 for(statname in all.gof.vars){
  if ('model' == statname) {

   nstats <- length(x$obs.model)
   if( min(x$pval.model[,"MC p-value"]) <0) {
    pval.max <- max((1:nstats)[x$pval.model[1:nstats, "MC p-value"] < 1]) + 3
   }
   else {
    pval.max <- max((1:nstats)[x$obs.model[1:nstats] > 0]) + 3
   }

   if (is.finite(pval.max) & pval.max < nstats) {
        model <- c(1:pval.max)
    }
    else {
        model <- c(1:nstats)
    }
    if (plotlogodds) {
        odds <- x$psim.model
#       odds[odds==0] <- NA
        odds[!is.na(odds)] <- log(odds[!is.na(odds)]/(1 - odds[!is.na(odds)]))
        odds.obs <- x$pobs.model
#       odds.obs[odds.obs==0] <- NA
        odds.obs[!is.na(odds.obs)] <- log(odds.obs[!is.na(odds.obs)]/(1 - odds.obs[!is.na(odds.obs)]))
        odds.bds <- x$bds.model
#       odds.bds[odds.bds==0] <- NA
        odds.bds[!is.na(odds.bds)] <- log(odds.bds[!is.na(odds.bds)]/(1 - odds.bds[!is.na(odds.bds)]))
        mininf <- min(min(odds[!is.infinite(odds)]),min(odds.obs[!is.infinite(odds.obs)]),min(odds.bds[!is.infinite(odds.bds)]))
        maxinf <- max(max(odds[!is.infinite(odds)]),max(odds.obs[!is.infinite(odds.obs)]),max(odds.bds[!is.infinite(odds.bds)]))

        odds[is.infinite(odds)&odds>0] <- maxinf
        odds[is.infinite(odds)&odds<0] <- mininf
        odds.obs[is.infinite(odds.obs)&odds.obs>0] <- maxinf
        odds.obs[is.infinite(odds.obs)&odds.obs<0] <- mininf
        odds.bds[is.infinite(odds.bds)&odds.bds>0] <- maxinf
        odds.bds[is.infinite(odds.bds)&odds.bds<0] <- mininf

        out <- odds
        out.obs <- odds.obs
        out.bds <- odds.bds

        ylab <- "log-odds for the statistic"
    }
    else {
        out <- x$psim.model
        out.obs <- x$pobs.model
        out.bds <- x$bds.model
        ylab <- "simulated quantiles"
    }
    pnames <- names(out.obs)
    ymin <- min(min(out,na.rm=TRUE),min(out.obs,na.rm=TRUE))
    ymax <- max(max(out,na.rm=TRUE),max(out.obs,na.rm=TRUE))

    boxplot(data.frame(out[, model]), xlab = "model statistics", 
            ylab = ylab, names = pnames, cex.axis = cex.axis, outline=FALSE,
            ylim=c(ymin,ymax), ...)

    points(seq(along = model), out.bds[1,model], pch = 1,cex=0.75)
    points(seq(along = model), out.bds[2,model], pch = 1,cex=0.75)
    lines(seq(along = model), out.bds[1, model], pch = 18,lty=1,lwd=1,col=color)
    lines(seq(along = model), out.bds[2, model], pch = 18,lty=1,lwd=1,col=color)
    points(seq(along = model), out.obs[model], pch = 16,cex=0.75)
    lines(seq(along = model), out.obs[model], lty = 1,lwd=3)
  }

 ###degree####

  if ('degree' == statname) {

   if( min(x$pval.deg[,"MC p-value"]) <0) {
    pval.max <- max((1:(n - 1))[x$pval.deg[1:(n - 1), "MC p-value"] < 1]) + 3
   }
   else {
    pval.max <- max((1:(n - 1))[x$obs.deg[1:(n - 1)] > 0]) + 3
   }

   if (is.finite(pval.max) & pval.max < n) {
        deg <- c(1:pval.max)
    }
    else {
        deg <- c(1:n)
    }
    if (plotlogodds) {
        odds <- x$psim.deg
#       odds[odds==0] <- NA
        odds[!is.na(odds)] <- log(odds[!is.na(odds)]/(1 - odds[!is.na(odds)]))
        odds.obs <- x$pobs.deg
#       odds.obs[odds.obs==0] <- NA
        odds.obs[!is.na(odds.obs)] <- log(odds.obs[!is.na(odds.obs)]/(1 - odds.obs[!is.na(odds.obs)]))
        odds.bds <- x$bds.deg
#       odds.bds[odds.bds==0] <- NA
        odds.bds[!is.na(odds.bds)] <- log(odds.bds[!is.na(odds.bds)]/(1 - odds.bds[!is.na(odds.bds)]))
        mininf <- min(min(odds[!is.infinite(odds)]),min(odds.obs[!is.infinite(odds.obs)]),min(odds.bds[!is.infinite(odds.bds)]))
        maxinf <- max(max(odds[!is.infinite(odds)]),max(odds.obs[!is.infinite(odds.obs)]),max(odds.bds[!is.infinite(odds.bds)]))

        odds[is.infinite(odds)&odds>0] <- maxinf
        odds[is.infinite(odds)&odds<0] <- mininf
        odds.obs[is.infinite(odds.obs)&odds.obs>0] <- maxinf
        odds.obs[is.infinite(odds.obs)&odds.obs<0] <- mininf
        odds.bds[is.infinite(odds.bds)&odds.bds>0] <- maxinf
        odds.bds[is.infinite(odds.bds)&odds.bds<0] <- mininf

        out <- odds
        out.obs <- odds.obs
        out.bds <- odds.bds

        ylab <- "log-odds for a node"
    }
    else {
        out <- x$psim.deg
        out.obs <- x$pobs.deg
        out.bds <- x$bds.deg
        ylab <- "proportion of nodes"
    }
    pnames <- c(deg)-1
    ymin <- min(min(out,na.rm=TRUE),min(out.obs,na.rm=TRUE))
    ymax <- max(max(out,na.rm=TRUE),max(out.obs,na.rm=TRUE))

    boxplot(data.frame(out[, deg]), xlab = "degree", 
            ylab = ylab, names = pnames, cex.axis = cex.axis, outline=FALSE,
            ylim=c(ymin,ymax), ...)

    points(seq(along = deg), out.bds[1,deg], pch = 1,cex=0.75)
    points(seq(along = deg), out.bds[2,deg], pch = 1,cex=0.75)
    lines(seq(along = deg), out.bds[1, deg], pch = 18,lty=1,lwd=1,col=color)
    lines(seq(along = deg), out.bds[2, deg], pch = 18,lty=1,lwd=1,col=color)
    points(seq(along = deg), out.obs[deg], pch = 16,cex=0.75)
    lines(seq(along = deg), out.obs[deg], lty = 1,lwd=3)
  }

  ###odegree####

  if ('odegree' == statname) {

   if( min(x$pval.odeg[,"MC p-value"]) <0) {
    pval.max <- max((1:(n - 1))[x$pval.odeg[1:(n - 1), "MC p-value"] < 1]) + 3
   }
   else {
    pval.max <- max((1:(n - 1))[x$obs.odeg[1:(n - 1)] > 0]) + 3
   }

   if (is.finite(pval.max) & pval.max < n) {
        odeg <- c(1:pval.max)
    }
    else {
        odeg <- c(1:n)
    }
    if (plotlogodds) {
        odds <- x$psim.odeg
#       odds[odds==0] <- NA
        odds[!is.na(odds)] <- log(odds[!is.na(odds)]/(1 - odds[!is.na(odds)]))
        odds.obs <- x$pobs.odeg
#       odds.obs[odds.obs==0] <- NA
        odds.obs[!is.na(odds.obs)] <- log(odds.obs[!is.na(odds.obs)]/(1 - odds.obs[!is.na(odds.obs)]))
        odds.bds <- x$bds.odeg
#       odds.bds[odds.bds==0] <- NA
        odds.bds[!is.na(odds.bds)] <- log(odds.bds[!is.na(odds.bds)]/(1 - odds.bds[!is.na(odds.bds)]))
        mininf <- min(min(odds[!is.infinite(odds)]),min(odds.obs[!is.infinite(odds.obs)]),min(odds.bds[!is.infinite(odds.bds)]))
        maxinf <- max(max(odds[!is.infinite(odds)]),max(odds.obs[!is.infinite(odds.obs)]),max(odds.bds[!is.infinite(odds.bds)]))

        odds[is.infinite(odds)&odds>0] <- maxinf
        odds[is.infinite(odds)&odds<0] <- mininf
        odds.obs[is.infinite(odds.obs)&odds.obs>0] <- maxinf
        odds.obs[is.infinite(odds.obs)&odds.obs<0] <- mininf
        odds.bds[is.infinite(odds.bds)&odds.bds>0] <- maxinf
        odds.bds[is.infinite(odds.bds)&odds.bds<0] <- mininf

        out <- odds
        out.obs <- odds.obs
        out.bds <- odds.bds

        ylab <- "log-odds for a node"
    }
    else {
        out <- x$psim.odeg
        out.obs <- x$pobs.odeg
        out.bds <- x$bds.odeg
        ylab <- "proportion of nodes"
    }
    pnames <- c(odeg)-1
    ymin <- min(min(out,na.rm=TRUE),min(out.obs,na.rm=TRUE))
    ymax <- max(max(out,na.rm=TRUE),max(out.obs,na.rm=TRUE))

    boxplot(data.frame(out[, odeg]), xlab = "out degree", 
            ylab = ylab, names = pnames, cex.axis = cex.axis, outline=FALSE,
            ylim=c(ymin,ymax), ...)

    points(seq(along = odeg), out.bds[1,odeg], pch = 1,cex=0.75)
    points(seq(along = odeg), out.bds[2,odeg], pch = 1,cex=0.75)
    lines(seq(along = odeg), out.bds[1, odeg], pch = 18,lty=1,lwd=1,col=color)
    lines(seq(along = odeg), out.bds[2, odeg], pch = 18,lty=1,lwd=1,col=color)
    points(seq(along = odeg), out.obs[odeg], pch = 16,cex=0.75)
    lines(seq(along = odeg), out.obs[odeg], lty = 1,lwd=3)
  }

  ###idegree####

  if ('idegree' == statname) {

   if( min(x$pval.ideg[,"MC p-value"]) <0) {
    pval.max <- max((1:(n - 1))[x$pval.ideg[1:(n - 1), "MC p-value"] < 1]) + 3
   }
   else {
    pval.max <- max((1:(n - 1))[x$obs.ideg[1:(n - 1)] > 0]) + 3
   }

   if (is.finite(pval.max) & pval.max < n) {
        ideg <- c(1:pval.max)
    }
    else {
        ideg <- c(1:n)
    }
    if (plotlogodds) {
        odds <- x$psim.ideg
#       odds[odds==0] <- NA
        odds[!is.na(odds)] <- log(odds[!is.na(odds)]/(1 - odds[!is.na(odds)]))
        odds.obs <- x$pobs.ideg
#       odds.obs[odds.obs==0] <- NA
        odds.obs[!is.na(odds.obs)] <- log(odds.obs[!is.na(odds.obs)]/(1 - odds.obs[!is.na(odds.obs)]))
        odds.bds <- x$bds.ideg
#       odds.bds[odds.bds==0] <- NA
        odds.bds[!is.na(odds.bds)] <- log(odds.bds[!is.na(odds.bds)]/(1 - odds.bds[!is.na(odds.bds)]))
        mininf <- min(min(odds[!is.infinite(odds)]),min(odds.obs[!is.infinite(odds.obs)]),min(odds.bds[!is.infinite(odds.bds)]))
        maxinf <- max(max(odds[!is.infinite(odds)]),max(odds.obs[!is.infinite(odds.obs)]),max(odds.bds[!is.infinite(odds.bds)]))

        odds[is.infinite(odds)&odds>0] <- maxinf
        odds[is.infinite(odds)&odds<0] <- mininf
        odds.obs[is.infinite(odds.obs)&odds.obs>0] <- maxinf
        odds.obs[is.infinite(odds.obs)&odds.obs<0] <- mininf
        odds.bds[is.infinite(odds.bds)&odds.bds>0] <- maxinf
        odds.bds[is.infinite(odds.bds)&odds.bds<0] <- mininf

        out <- odds
        out.obs <- odds.obs
        out.bds <- odds.bds

        ylab <- "log-odds for a node"
    }
    else {
        out <- x$psim.ideg
        out.obs <- x$pobs.ideg
        out.bds <- x$bds.ideg
        ylab <- "proportion of nodes"
    }
    pnames <- c(ideg)-1
    ymin <- min(min(out,na.rm=TRUE),min(out.obs,na.rm=TRUE))
    ymax <- max(max(out,na.rm=TRUE),max(out.obs,na.rm=TRUE))

    boxplot(data.frame(out[, ideg]), xlab = "in degree", 
            ylab = ylab, names = pnames, cex.axis = cex.axis, outline=FALSE,
            ylim=c(ymin,ymax), ...)

    points(seq(along = ideg), out.bds[1,ideg], pch = 1,cex=0.75)
    points(seq(along = ideg), out.bds[2,ideg], pch = 1,cex=0.75)
    lines(seq(along = ideg), out.bds[1, ideg], pch = 18,lty=1,lwd=1,col=color)
    lines(seq(along = ideg), out.bds[2, ideg], pch = 18,lty=1,lwd=1,col=color)
    points(seq(along = ideg), out.obs[ideg], pch = 16,cex=0.75)
    lines(seq(along = ideg), out.obs[ideg], lty = 1,lwd=3)
  }

  ###espart####

  if ('espartners' == statname) {

   pval.max <- max((1:(n - 1))[x$pval.espart[1:(n - 1), "MC p-value"] < 
        1]) + 3
    if (is.finite(pval.max) & pval.max < n) {
        espart <- c(1:pval.max)
    }
    else {
        espart <- c(1:(n-1))
    }
    if (plotlogodds) {
        odds <- x$psim.espart
#       odds[odds==0] <- NA
        odds[!is.na(odds)] <- log(odds[!is.na(odds)]/(1 - odds[!is.na(odds)]))
        odds.obs <- x$pobs.espart
#       odds.obs[odds.obs==0] <- NA
        odds.obs[!is.na(odds.obs)] <- log(odds.obs[!is.na(odds.obs)]/(1 - odds.obs[!is.na(odds.obs)]))
        odds.bds <- x$bds.espart
#       odds.bds[odds.bds==0] <- NA
        odds.bds[!is.na(odds.bds)] <- log(odds.bds[!is.na(odds.bds)]/(1 - odds.bds[!is.na(odds.bds)]))
        mininf <- min(min(odds[!is.infinite(odds)]),min(odds.obs[!is.infinite(odds.obs)]),min(odds.bds[!is.infinite(odds.bds)]))
        maxinf <- max(max(odds[!is.infinite(odds)]),max(odds.obs[!is.infinite(odds.obs)]),max(odds.bds[!is.infinite(odds.bds)]))
        odds[is.infinite(odds)&odds>0] <- maxinf
        odds[is.infinite(odds)&odds<0] <- mininf
        odds.obs[is.infinite(odds.obs)&odds.obs>0] <- maxinf
        odds.obs[is.infinite(odds.obs)&odds.obs<0] <- mininf
        odds.bds[is.infinite(odds.bds)&odds.bds>0] <- maxinf
        odds.bds[is.infinite(odds.bds)&odds.bds<0] <- mininf
        out <- odds
        out.obs <- odds.obs
        out.bds <- odds.bds

        ylab <- "log-odds for an edge"
    }
    else {
        out <- x$psim.espart
        out.obs <- x$pobs.espart
        out.bds <- x$bds.espart
        ylab <- "proportion of edges"
        mininf <- min(min(out),min(out.obs),min(out.bds))
        maxinf <- max(max(out),max(out.obs),max(out.bds))
    }
    pnames <- c(espart)-1
    ymin <- min(min(out,na.rm=TRUE),min(out.obs,na.rm=TRUE))
    ymax <- max(max(out,na.rm=TRUE),max(out.obs,na.rm=TRUE))

    boxplot(data.frame(out[, espart]), xlab = "edge-wise shared partners", 
            ylab = ylab, names = pnames, cex.axis = cex.axis, outline=FALSE,
            ylim=c(ymin,ymax), ...)

    points(seq(along = espart), out.bds[1,espart], pch = 1,cex=0.75)
    points(seq(along = espart), out.bds[2,espart], pch = 1,cex=0.75)
    lines(seq(along = espart), out.bds[1, espart], pch = 18,lty=1,lwd=1,col=color)
    lines(seq(along = espart), out.bds[2, espart], pch = 18,lty=1,lwd=1,col=color)
    points(seq(along = espart), out.obs[espart], pch = 16, cex=0.75)
    lines(seq(along = espart), out.obs[espart], lty = 1,lwd=3)

  }

  ###dspart####

  if ('dspartners' == statname) {
   pval.max <- max((1:(n - 1))[x$pval.dspart[1:(n - 1), "MC p-value"] < 
        1]) + 3
    if (is.finite(pval.max) & pval.max < n) {
        dspart <- c(1:pval.max)
    }
    else {
        dspart <- c(1:n)
    }
    if (plotlogodds) {
        odds <- x$psim.dspart
#       odds[odds==0] <- NA
        odds[!is.na(odds)] <- log(odds[!is.na(odds)]/(1 - odds[!is.na(odds)]))
        odds.obs <- x$pobs.dspart
#       odds.obs[odds.obs==0] <- NA
        odds.obs[!is.na(odds.obs)] <- log(odds.obs[!is.na(odds.obs)]/(1 - odds.obs[!is.na(odds.obs)]))
        odds.bds <- x$bds.dspart
#       odds.bds[odds.bds==0] <- NA
        odds.bds[!is.na(odds.bds)] <- log(odds.bds[!is.na(odds.bds)]/(1 - odds.bds[!is.na(odds.bds)]))
        mininf <- min(min(odds[!is.infinite(odds)]),min(odds.obs[!is.infinite(odds.obs)]),min(odds.bds[!is.infinite(odds.bds)]))
        maxinf <- max(max(odds[!is.infinite(odds)]),max(odds.obs[!is.infinite(odds.obs)]),max(odds.bds[!is.infinite(odds.bds)]))
        odds[is.infinite(odds)&odds>0] <- maxinf
        odds[is.infinite(odds)&odds<0] <- mininf
        odds.obs[is.infinite(odds.obs)&odds.obs>0] <- maxinf
        odds.obs[is.infinite(odds.obs)&odds.obs<0] <- mininf
        odds.bds[is.infinite(odds.bds)&odds.bds>0] <- maxinf
        odds.bds[is.infinite(odds.bds)&odds.bds<0] <- mininf
        out <- odds
        out.obs <- odds.obs
        out.bds <- odds.bds

        ylab <- "log-odds for an edge"
    }
    else {
        out <- x$psim.dspart
        out.obs <- x$pobs.dspart
        out.bds <- x$bds.dspart
        ylab <- "proportion of dyads"
        mininf <- min(min(out),min(out.obs),min(out.bds))
        maxinf <- max(max(out),max(out.obs),max(out.bds))
    }
    pnames <- c(dspart)-1
    ymin <- min(min(out,na.rm=TRUE),min(out.obs,na.rm=TRUE))
    ymax <- max(max(out,na.rm=TRUE),max(out.obs,na.rm=TRUE))

    boxplot(data.frame(out[, dspart]), xlab = "dyad-wise shared partners", 
            ylab = ylab, names = pnames, cex.axis = cex.axis, outline=FALSE,
            ylim=c(ymin,ymax), ...)

    points(seq(along = dspart), out.bds[1,dspart], pch = 1,cex=0.75)
    points(seq(along = dspart), out.bds[2,dspart], pch = 1,cex=0.75)
    lines(seq(along = dspart), out.bds[1, dspart], pch = 18,lty=1,lwd=1,col=color)
    lines(seq(along = dspart), out.bds[2, dspart], pch = 18,lty=1,lwd=1,col=color)
    points(seq(along = dspart), out.obs[dspart], pch = 16,cex=0.75)
    lines(seq(along = dspart), out.obs[dspart], lty = 1,lwd=3)
  }

  ###triadcensus####

  if ('triadcensus' == statname) {

    if (plotlogodds) {
        odds <- x$psim.triadcensus
#       odds[odds==0] <- NA
        odds[!is.na(odds)] <- log(odds[!is.na(odds)]/(1 - odds[!is.na(odds)]))
        odds.obs <- x$pobs.triadcensus
#       odds.obs[odds.obs==0] <- NA
        odds.obs[!is.na(odds.obs)] <- log(odds.obs[!is.na(odds.obs)]/(1 - odds.obs[!is.na(odds.obs)]))
        odds.bds <- x$bds.triadcensus
#       odds.bds[odds.bds==0] <- NA
        odds.bds[!is.na(odds.bds)] <- log(odds.bds[!is.na(odds.bds)]/(1 - odds.bds[!is.na(odds.bds)]))
        mininf <- min(min(odds[!is.infinite(odds)]),min(odds.obs[!is.infinite(odds.obs)]),min(odds.bds[!is.infinite(odds.bds)]))
        maxinf <- max(max(odds[!is.infinite(odds)]),max(odds.obs[!is.infinite(odds.obs)]),max(odds.bds[!is.infinite(odds.bds)]))
        odds[is.infinite(odds)&odds>0] <- maxinf
        odds[is.infinite(odds)&odds<0] <- mininf
        odds.obs[is.infinite(odds.obs)&odds.obs>0] <- maxinf
        odds.obs[is.infinite(odds.obs)&odds.obs<0] <- mininf
        odds.bds[is.infinite(odds.bds)&odds.bds>0] <- maxinf
        odds.bds[is.infinite(odds.bds)&odds.bds<0] <- mininf
        out <- odds
        out.obs <- odds.obs
        out.bds <- odds.bds

        ylab <- "log-odds for a triad"
    }
    else {
        out <- x$psim.triadcensus
        out.obs <- x$pobs.triadcensus
        out.bds <- x$bds.triadcensus
        ylab <- "proportion of triads"
        mininf <- min(min(out),min(out.obs),min(out.bds))
        maxinf <- max(max(out),max(out.obs),max(out.bds))
    }
    triadcensus <- dimnames(x$sim.triadcensus)[[2]]
    pnames <- dimnames(x$sim.triadcensus)[[2]]
    ymin <- min(min(out,na.rm=TRUE),min(out.obs,na.rm=TRUE))
    ymax <- max(max(out,na.rm=TRUE),max(out.obs,na.rm=TRUE))

    boxplot(data.frame(out[, triadcensus]), xlab = "triad census", 
            ylab = ylab, names = pnames, cex.axis = cex.axis, outline=FALSE,
            ylim=c(ymin,ymax), ...)

    points(seq(along = triadcensus), out.bds[1,triadcensus], pch = 1,cex=0.75)
    points(seq(along = triadcensus), out.bds[2,triadcensus], pch = 1,cex=0.75)
    lines(seq(along = triadcensus), out.bds[1, triadcensus], pch = 18,lty=1,lwd=1,col=color)
    lines(seq(along = triadcensus), out.bds[2, triadcensus], pch = 18,lty=1,lwd=1,col=color)
    points(seq(along = triadcensus), out.obs[triadcensus], pch = 16, cex=0.75)
    lines(seq(along = triadcensus), out.obs[triadcensus], lty = 1,lwd=3)

  }

  ###distance####

  if ('distance' == statname) {

    pval.max <- max((1:(n - 1))[x$pval.dist[1:(n - 1), "MC p-value"] < 
        1]) + 3
    if (is.finite(pval.max) & pval.max < n) {
        dist <- c(1:pval.max, n)
    }
    else {
        dist <- c(1:n)
    }
    pnames <- paste(dist)
    pnames[length(dist)] <- "NR"
    if (plotlogodds) {
        odds <- x$psim.dist
#       odds[odds==0] <- NA
        odds[!is.na(odds)] <- log(odds[!is.na(odds)]/(1 - odds[!is.na(odds)]))
        odds.obs <- x$pobs.dist
#       odds.obs[odds.obs==0] <- NA
        odds.obs[!is.na(odds.obs)] <- log(odds.obs[!is.na(odds.obs)]/(1 - odds.obs[!is.na(odds.obs)]))
        odds.bds <- x$bds.dist
#       odds.bds[odds.bds==0] <- NA
        odds.bds[!is.na(odds.bds)] <- log(odds.bds[!is.na(odds.bds)]/(1 - odds.bds[!is.na(odds.bds)]))
        oodds <- is.infinite(odds) | is.na(odds)
        oodds.obs <- is.infinite(odds.obs) | is.na(odds.obs)
        oodds.bds <- is.infinite(odds.bds) | is.na(odds.bds)
        mininf <- min(min(odds[!oodds]),min(odds.obs[!oodds.obs]),min(odds.bds[!oodds.bds]))
        maxinf <- max(max(odds[!oodds]),max(odds.obs[!oodds.obs]),max(odds.bds[!oodds.bds]))
        odds[is.infinite(odds)&odds>0] <- maxinf
        odds[is.infinite(odds)&odds<0] <- mininf
        odds.obs[is.infinite(odds.obs)&odds.obs>0] <- maxinf
        odds.obs[is.infinite(odds.obs)&odds.obs<0] <- mininf
        odds.bds[is.infinite(odds.bds)&odds.bds>0] <- maxinf
        odds.bds[is.infinite(odds.bds)&odds.bds<0] <- mininf
        odds.bds[1,][is.na(odds.bds[1,])] <- mininf
        odds.bds[2,][is.na(odds.bds[2,])] <- maxinf
        out <- odds
        out.obs <- odds.obs
        out.bds <- odds.bds

        ylab <- "log-odds for a dyad"
    }
    else {
        out <- x$psim.dist
        out.obs <- x$pobs.dist
        out.bds <- x$bds.dist
        ylab <- "proportion of dyads"
        mininf <- min(min(out),min(out.obs),min(out.bds))
        maxinf <- max(max(out),max(out.obs),max(out.bds))
    }

    if(normalize.reachability){
      mdist <- max(dist,na.rm=TRUE)
      totrange <- range(out.bds[1,][out.bds[1,] > out.bds[1,mdist]],
                        out.bds[2,][out.bds[2,] < out.bds[2,mdist]])
      out[,mdist] <- (out[,mdist]-out.bds[1,mdist]) * 
        diff(totrange) / diff(out.bds[,mdist]) + totrange[1]
      out.obs[mdist] <- (out.obs[mdist]- out.bds[1,mdist]) *
        diff(totrange) / diff(out.bds[,mdist]) + totrange[1]
      out.bds[,mdist] <- totrange
    }

    ymin <- min(min(out,na.rm=TRUE),min(out.obs,na.rm=TRUE))
    ymax <- max(max(out,na.rm=TRUE),max(out.obs,na.rm=TRUE))
    if(!plotlogodds){
     ymin <- max(0,ymin)
     ymax <- min(1,ymax)
    }

    boxplot(data.frame(out[, dist]), xlab = "minimum geodesic distance", 
            ylab = ylab, names = pnames, cex.axis = cex.axis, outline=FALSE,
            ylim=c(ymin,ymax), ...)

    points(seq(along = dist), out.bds[1,dist], pch = 1,cex=0.75)
    points(seq(along = dist), out.bds[2,dist], pch = 1,cex=0.75)
    lines(seq(along = dist)[-length(dist)], out.bds[1, dist][-length(dist)], pch = 18,lty=1,lwd=1,col=color)
    lines(seq(along = dist)[-length(dist)], out.bds[2, dist][-length(dist)], pch = 18,lty=1,lwd=1,col=color)
    points(seq(along = dist), out.obs[dist], pch = 16,cex=0.75)
    lines(seq(along = dist)[-length(dist)], out.obs[dist][-length(dist)],
		 lty = 1,lwd=3)
    }
   }

   mtext(main,side=3,outer=TRUE,cex=1.5,padj=2)
   invisible()
}



#ergm.get.terms.formula <- function(formula){
# trms <- all.names(formula)
# ntrms <- length(trms)
# if(ntrms == 2*trunc(ntrms/2)){
#   ntrms <- ntrms/2
#  }else{
#   ntrms <- (ntrms+1)/2
#  }
# trms[-c(1:ntrms)]
#}



ergm.rhs.formula <- function(formula){
#all.vars(ergm.update.formula(formula, .~0)) 
 unlist(dimnames(attr(terms(formula),"factors"))[-1])
}
