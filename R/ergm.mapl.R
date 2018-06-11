#  File R/ergm.mapl.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2018 Statnet Commons
#######################################################################
###############################################################################
# The <ergm.mapl> function creates an initial fit for a specified formula
#
# --PARAMETERS--
#   formula     :  a formula of the form 'nw ~ model term(s)'
#   init      :  a vector of starting values for estimation, or optionally
#                  if these are to be estimated, the string "MPLE";
#                  default="MPLE"
#   nsim        :  the number of simulations to use in forming the initial
#                   fit
#   burnin      :  the number of proposals to ignore before MCMC sampling
#                  begins; default=10,000
#   maxit       :  the number of MCMC parameter updates to the value
#                  maximizing the MCMC likelihood; default=3
#   constraints :  a one-sided formula of the constraint terms; options are
#                      bd        degrees        nodegrees
#                      edges     degreedist     idegreedist
#                      observed  odegreedist
#                  default="~ ."
#   proposaltype:  presumably the proposal type, but this is only used
#                  in calls to <ergm.san>, which doesn't accept a
#                  'proposaltype' argument
#   target.stats   :  a vector of the mean value parameters;
#                  default=the observed statistic from the 'nw' in formula
#   control     :  a list of control parameters returned from <control.ergm>;
#                  default=control.ergm()
#   tau         :  ??, is passed along to <ergm.san> where it is ignored;
#                  see the <ergm.san> header for details; default=1
#   invcov      :  the initial inverse covariance matrix used to calculate
#                  the Mahalanobis distance; default=NULL
#   verbose     :  whether ergm should be verbose (T or F); default=FALSE
#
#
# --RETURNED--
#   v: an ergm object as a list containing several items; for details see
#      the return list in the <ergm> function header (<ergm.mapl>= #);
#
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
  check.control.class("ergm", "ergm.mapl")
  control.toplevel(...,myname="ergm")

  if (verbose) message("Evaluating network in model.")

  nw <- ergm.getnetwork(formula)
  
  if (verbose) message("Fitting initial model.")

  proposalclass <- "c"
    
  if(control$drop){
   model.initial <- ergm_model(formula, nw, initialfit=TRUE, term.options=control$term.options)
   obs.stats <- NVL(target.stats, summary(formula))
   extremeval <- +(model.initial$maxval==obs.stats)-(model.initial$minval==obs.stats)
   model.initial$etamap$offsettheta[extremeval!=0] <- TRUE
  }else{
#    model.initial <- ergm_model(formula, nw, response=response, initialfit=TRUE)
    model.initial <- ergm_model(formula, nw, term.options=control$term.options)
    extremeval <- rep(0, length=length(model.initial$etamap$offsettheta))
  }

  # MPLE & Meanstats -> need fake network
  if(!missing(target.stats)){
    nw<-san(formula, target.stats=target.stats, 
            constraints=~., #constraints=constraints,
            proposaltype=proposaltype,
            tau=tau, invcov=invcov, burnin=10*burnin, verbose=verbose)
    if(verbose){message_print(summary(formula, basis=nw)-target.stats)}
  }
  
  fd.initial <- ergm.design(nw, verbose=verbose)
  initcopy <- init
  
  pl <- ergm.pl(nw=nw, fd=fd.initial,
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
    if(verbose){message_print(summary(formula, basis=sim)-target.stats)}
    if(verbose){message_print(sum(sim[,] != nw[,]))}
    fd.initial <- ergm.design(sim, verbose=verbose)
    sim.pl <- ergm.pl(nw=nw, fd=fd.initial,
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
