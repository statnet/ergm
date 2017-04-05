#  File ergm/R/ergm.mapl.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
###############################################################################
# The <ergm.mapl> function creates an initial fit for a specified formula
################################################################################

ergm.mapl <- function(formula, init="MPLE", 
                 nsim=25,
                 burnin=10000,
                 maxit=3,
                 constraints=~.,
                 proposaltype="TNT10",
                 target.stats=NULL,
                 control=control.ergm(MPLEtype="penalized"),
                 tau=1, invcov=NULL,
                 verbose=FALSE, ...) {

  current.warn <- options()$warn
  options(warn=0)
  if (verbose) cat("Evaluating network in model\n")

  nw <- ergm.getnetwork(formula)
  
  if (verbose) cat("Fitting initial model.\n")

  proposalclass <- "c"
    
  if(control$drop){
   model.initial <- ergm.getmodel(formula, nw, initialfit=TRUE)
#   obs.stats <- if(!is.null(target.stats)) target.stats else summary(formula)
   obs.stats <- if(!is.null(target.stats)) target.stats else summary(formula)
   extremeval <- +(model.initial$maxval==obs.stats)-(model.initial$minval==obs.stats)
   model.initial$etamap$offsettheta[extremeval!=0] <- TRUE
  }else{
#    model.initial <- ergm.getmodel(formula, nw, initialfit=TRUE)
    model.initial <- ergm.getmodel(formula, nw)
    extremeval <- rep(0, length=length(model.initial$etamap$offsettheta))
  }

  # MPLE & Meanstats -> need fake network
  if(!missing(target.stats)){
    nw<-san(formula, target.stats=target.stats, 
            constraints=~., #constraints=constraints,
            proposaltype=proposaltype,
            tau=tau, invcov=invcov, burnin=10*burnin, verbose=verbose)
    if(verbose){print(summary(formula, basis=nw)-target.stats)}
  }
  
  Clist.initial <- ergm.Cprepare(nw, model.initial)
  Clist.miss.initial <- ergm.design(nw, model.initial, verbose=verbose)
  Clist.initial$target.stats=target.stats
  initcopy <- init
  
  pl <- ergm.pl(Clist=Clist.initial, Clist.miss=Clist.miss.initial,
                m=model.initial,theta.offset=ifelse(extremeval!=0,extremeval*Inf,NA),
                verbose=verbose)
  initialfit <- ergm.maple(pl=pl, model.initial,
                           MPLEtype=control$MPLE.type, 
                           verbose=verbose, ...)
  if("MPLE" %in% init){init <- initialfit$coef}

  if(nsim>0){
   if(missing(target.stats)){
    target.stats <- summary(formula, basis=nw)
    sim<-simulate(formula, constraints=constraints,
                  init=init, burnin=burnin,
                  verbose=verbose)
   }else{
    sim <- nw
   }

   for(i in 1:nsim){
    sim<-simulate(formula, constraints=constraints,
                  init=init, burnin=burnin,
                  verbose=verbose)
    sim <- san(formula, target.stats=target.stats, verbose=verbose,
               proposaltype=proposaltype,
               tau=tau, invcov=invcov, burnin=burnin, 
               constraints=constraints, basis=sim)
    if(verbose){print(summary(formula, basis=sim)-target.stats)}
    if(verbose){print(sum(sim[,] != nw[,]))}
    Clist.initial <- ergm.Cprepare(sim, model.initial)
    Clist.miss.initial <- ergm.design(sim, model.initial, verbose=verbose)
    Clist.initial$target.stats=target.stats
    sim.pl <- ergm.pl(Clist=Clist.initial, Clist.miss=Clist.miss.initial,
                      m=model.initial,theta.offset=ifelse(extremeval!=0,extremeval*Inf,NA),
                      verbose=verbose)
    pl$zy <- c(pl$zy,sim.pl$zy)
    pl$foffset <- c(pl$foffset,sim.pl$foffset)
    pl$xmat <- rbind(pl$xmat,sim.pl$xmat)
    pl$wend <- c(pl$wend,sim.pl$wend)
    pl$zy.full <- c(pl$zy.full,sim.pl$zy.full)
    pl$foffset.full <- c(pl$foffset.full,sim.pl$foffset.full)
    pl$xmat.full <- rbind(pl$xmat.full,sim.pl$xmat.full)
    pl$wend.full <- c(pl$wend.full,sim.pl$wend.full)
   }
   pl$wend <- pl$wend / nsim
   pl$wend.full <- pl$wend.full / nsim
   initialfit <- ergm.maple(pl=pl, model.initial,
                            MPLEtype=control$MPLE.type, 
                            verbose=verbose, ...)
  }

  initialfit$offset <- model.initial$etamap$offsettheta
  initialfit$drop <- extremeval
  initialfit$network <- nw
  initialfit$newnetwork <- nw
  initialfit$formula <- formula
  initialfit$constraints <- constraints
  initialfit$prop.args <- control$MCMC.prop.args
  initialfit$prop.weights <- control$MCMC.prop.weights
  initialfit
}
