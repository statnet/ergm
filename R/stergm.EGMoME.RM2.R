#=============================================================================
# This file contains the 2 following function for fitting a stergm, using the
# Robbins-Monro approach
#       <stergm.RM>
#       <stergm.phase12.C>
#=============================================================================




################################################################################
# The <stergm.RM> function fits a stergm using the Robbins-Monro style of
# methods of moments estimation; this style should only be used if it is known
# a priori that the derivative of each element of the equilibrium expected
# values of the statistics of interest with respect to the corresponding formation
# phase parameter is positive.
#
# --PARAMETERS--
#   theta.form0    : the initial theta formation coefficients
#   nw             : a network object
#   model.form     : a formation model, as returned by <ergm.getmodel>
#   model.diss     : a dissolution model, as returned by <ergm.getmodel>
#   theta.diss     : the initial theta dissolution coefficients
#   control     : the list of parameters which tune the MCMC sampling
#                    processes; the recognized components of 'control'
#                    are those passed to and used by <stergm.phase12.C>
#                    and are described in its function header
#   MHproposal.form: a MHproposal object for the formation process, as
#                    returned by <getMHproposal>
#   MHproposal.diss: a MHproposal object for the dissolution process, as
#                    returned by <getMHproposal>
#   verbose        : whether this and subsequently called R and C code
#                    should be verbose (T or F); default=FALSE
#
# --RETURNED--
#   a stergm object as a list containing:
#    coef.form :   the estimated formation coefficients
#    coef.diss :   the estimated dissolution coefficients
#    newnetwork:   the 'nw' inputted into this function
#    network   :   the 'nw' inputted into this function
#    theta.form.original: the 'theta.form0' inputted into this function
#
################################################################################

stergm.RM2 <- function(theta.form0, nw, model.form, model.diss, 
                            theta.diss,
                            control, MHproposal.form, MHproposal.diss,
                            verbose=FALSE){

  control$nw.diff <- model.form$obs - model.form$target.stats # nw.diff keeps track of the difference between the current network and the target statistics.

  if(verbose) cat("Robbins-Monro algorithm with coef_F_0 = (",theta.form0, ") and coef_D = (",theta.diss,")\n" )
  eta.form <- ergm.eta(theta.form0, model.form$etamap)
  eta.diss <- ergm.eta(theta.diss, model.diss$etamap)

  p <- length(eta.form)
  
  control.phase1<-control
  control.phase1$time.samplesize <- (control$RM.phase1n_base + 3 * p)*control$RM.interval
  control.phase1$time.burnin <- control$RM.burnin
  control.phase1$time.interval <- 1

  # Run Phase 1.
  z <- stergm.getMCMCsample(nw, model.form, model.diss, MHproposal.form, MHproposal.diss, eta.form, eta.diss, control.phase1, verbose)

  nw <- z$newnetwork
  burnin.stats <- sweep(z$statsmatrix.form, 2, control$nw.diff, "+")
  control$nw.diff <- control$nw.diff + z$statsmatrix.form[NROW(z$statsmatrix.form),]
 

  # Not assuming any relationship between parameters and statistics.

  control$invGradient <- matrix(0,p,p)
  
  control$phase2n <- p+7+control$RM.phase2n_base
    
  control$jitter<-rep(min(1/apply(burnin.stats,2,sd))/sqrt(control$phase2n),p)


  oh <- NULL
  
  for(subphase in 1:control$RM.phase2sub){
    for(retry in 1:control$RM.phase2sub_retries){
      
      z <- stergm.EGMoME.RM.Phase2.C(nw, model.form, model.diss, MHproposal.form, MHproposal.diss,
                                     eta.form, eta.diss, control, verbose=verbose)
      nw <- z$newnetwork
      control$nw.diff<-z$nw.diff
      eta.form <- z$coef.form
      oh <- rbind(oh,z$objective.history) # Pool the history.
      
      # Regress statistics on parameters.
      # This uses GLS to account for serial correlation in statistics, and more recent are weighted higher.
      # First row is the intercept.

    oh.wt <- exp(control$RM.grad_decay*seq_len(NROW(oh)))
    x<-oh[,1:p] # #$%^$ gls() doesn't respect I()...
    oh.fit <- apply(oh[,-(1:p)],2,function(y){
      a<-try(coef(gls(y~x,weight=varFixed(~1/oh.wt),correlation=corAR1())))
      if(inherits(a,"try-error")) rep(NA,p+1) else a
    }
                    )
      grad.OK <- !any(is.na(oh.fit))
      if(all(is.na(oh.fit))) error("The search is trapped. Try a different starting value.")
      
      oh.fit[is.na(oh.fit)] <- 0    
      control$invGradient <- robust.inverse(oh.fit[-1,]) * control$RM.init_gain * 2^-subphase

      control$jitter<-rep(0,p)
      if(!grad.OK){
        message("Only partial gradient matrix computed. Redoing the subphase.")
                    
 }else break
    }
    
    control$phase2n <- round(2.52*(control$phase2n-control$phase2n_base)+control$phase2n_base);

  }
  
  #ve<-with(z,list(coef=eta,sample=s$statsmatrix.form,sample.obs=NULL))
  ve<-with(z,list(coef.form=coef.form,coef.diss=theta.diss,objective.history=oh))
  names(ve$coef.form)<-model.form$coef.names
  
  #endrun <- control$MCMC.burnin+control$MCMC.interval*(ve$samplesize-1)
  #attr(ve$sample, "mcpar") <- c(control$MCMC.burnin+1, endrun, control$MCMC.interval)
  #attr(ve$sample, "class") <- "mcmc"
  
  c(ve, list(newnetwork=nw, 
                       init.form=theta.form0,
                       #interval=control$MCMC.interval, burnin=control$MCMC.burnin, 
                       network=nw))
            
}


