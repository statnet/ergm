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

stergm.RM <- function(theta.form0, nw, model.form, model.diss, 
                            theta.diss,
                            control, MHproposal.form, MHproposal.diss,
                            verbose=FALSE){

  control$nw.diff <- model.form$obs - model.form$target.stats # nw.diff keeps track of the difference between the current network and the target statistics.

  if(verbose) cat("Robbins-Monro algorithm with coef_F_0 = (",theta.form0, ") and coef_D = (",theta.diss,")\n" )
  eta.form <- ergm.eta(theta.form0, model.form$etamap)
  eta.diss <- ergm.eta(theta.diss, model.diss$etamap)

  control.phase1<-control
  control.phase1$time.samplesize <- (control$RM.phase1n_base + 3 * length(eta.form))*control$RM.interval
  control.phase1$time.burnin <- control$RM.burnin
  control.phase1$time.interval <- 1

  # Run Phase 1.
  z <- stergm.getMCMCsample(nw, model.form, model.diss, MHproposal.form, MHproposal.diss, eta.form, eta.diss, control.phase1, verbose)

  nw <- z$newnetwork
  burnin.stats <- sweep(z$statsmatrix.form, 2, control$nw.diff, "+")
  control$nw.diff <- control$nw.diff + z$statsmatrix.form[NROW(z$statsmatrix.form),]
 

  # First subphase -- "gradient" matrix is just the covariance matrix.
  control$invGradient <-
    if(length(eta.form)==1) 1/sqrt(var(burnin.stats)) * control$RM.init_gain
    else diag(1/sqrt(diag(cov(burnin.stats)))) * control$RM.init_gain
  control$phase2n <- length(eta.form)+7+control$phase2n_base

  for(subphase in 1:control$RM.phase2sub){
    z <- stergm.EGMoME.RM.Phase2.C(nw, model.form, model.diss, MHproposal.form, MHproposal.diss,
                                   eta.form, eta.diss, control, verbose=verbose)
    nw <- z$newnetwork
    control$nw.diff<-z$nw.diff
    eta.form <- z$coef.form
    oh <- z$objective.history
    control$phase2n <- round(2.52*(control$phase2n-control$phase2n_base)+control$phase2n_base);

    # Regress statistics on parameters.
    # First row is the intercept.
    oh.fit <- coef(lm(oh[,-(1:length(eta.form))]~oh[,1:length(eta.form)]))
    control$invGradient <- robust.inverse(oh.fit[-1,]) * control$RM.init_gain * 2^-subphase
  }
  
  #ve<-with(z,list(coef=eta,sample=s$statsmatrix.form,sample.obs=NULL))
  ve<-with(z,list(coef.form=coef.form,coef.diss=theta.diss,objective.history=objective.history))
  names(ve$coef.form)<-model.form$coef.names
  
  #endrun <- control$MCMC.burnin+control$MCMC.interval*(ve$samplesize-1)
  #attr(ve$sample, "mcpar") <- c(control$MCMC.burnin+1, endrun, control$MCMC.interval)
  #attr(ve$sample, "class") <- "mcmc"
  
  c(ve, list(newnetwork=nw, 
                       init.form=theta.form0,
                       #interval=control$MCMC.interval, burnin=control$MCMC.burnin, 
                       network=nw))
            
}







################################################################################
# The <stergm.phase12.C> function is basically a wrapper for <MCMCDynPhase12.c>,
# which does the Rob-Mon sampling and estimation 
#
# --PARAMETERS--
#   g              : a network object
#   target.stats      : the mean statistics to be subtracted from the observed
#                    statistics
#   model.form     : a formation model, as returned by <ergm.getmodel>
#   model.diss     : a dissolution model, as returned by <ergm.getmodel>
#   MHproposal.form: a MHproposal object for the formation process, as
#                    returned by <getMHproposal>
#   MHproposal.diss: a MHproposal object for the dissolution process, as
#                    returned by <getMHproposal>
#   eta.form0      : the initial and canonical eta formation parameters
#   eta.diss       : the initial and canonical eta dissolution parameters
#   control     : the list of parameters which tune the MCMC sampling
#                    processes; recognized components include:
#       target.stats      : presumably the mean statistics, but this isn't used
#                        other than to return it
#       maxchanges     : 5 times the maximum number of changes to allocate 
#                        space for; this value is ignored if the number of edges 
#                        in 'g' is greater than 'maxchanges'; this value is 
#                        divided by 5 before being used to deteremine the max
#                        number of changes
#       RM.init_gain   : this is only used to adjust 'aDdiaginv'in phase1,
#                        in particular:
#                             aDdiaginv = gain/sqrt(aDdiaginv)
#       RM.phase1n_base: this helps define the 'phase1n' param, which in turn
#                        multiplies 'RM.interval' to control the number of
#                        phase1 iterations; this is the base portion of 'phase1n',
#                        which is added to 3*(the number of formation coefficients)
#                        to form 'phase1n'
#       RM.phase2sub   : phase2 is a 3-deep nested for-loop and 'RM.phase2sub' limits
#                        the outer loop counter
#       RM.phase2n_base: this helps define the 'phase2n' param, which in turn
#                        limits the phase2 middle loop counter; this is the
#                        base portion of 'phase2n', which is added to 7+(the number
#                        of formation coefficients) to form 'phase2n'
#       RM.burnin      : the number of MCMC steps to disregard for the burn-in
#                        period
#       RM.interval    : like the SPSA.interval, this seems a little more like
#                        a sample size, than an interval, it helps control the 
#                        number of MCMCsteps used in phase1 and phase2; in
#                        phase2, this limits the innermost loop counter
#       MH.burnin      : this is received as MH_interval and is used to
#                        control the number of proposals in each MCMC step
#   verbose        : whether this and subsequently called R and C code
#                    should be verbose (T or F); default=FALSE
#   
# --RETURNED--
#   a list with the 2 following components:
#      target.stats: the 'target.stats' from the 'control'; note that this is NOT
#                 the 'target.stats' inputted directly to this function
#      eta.form : the estimated? eta formation coefficients
#
################################################################################

stergm.EGMoME.RM.Phase2.C <- function(nw, model.form, model.diss, 
                             MHproposal.form, MHproposal.diss, eta.form0, eta.diss,
                             control, verbose) {
  Clist.form <- ergm.Cprepare(nw, model.form)
  Clist.diss <- ergm.Cprepare(nw, model.diss)
  maxedges <- max(control$MCMC.init.maxedges, Clist.form$nedges)
  if(verbose){cat(paste("MCMCDyn workspace is",maxedges,"\n"))}
  
  z <- .C("MCMCDynRMPhase2_wrapper",
          # Observed/starting network. 1
          as.integer(Clist.form$tails), as.integer(Clist.form$heads), 
          as.integer(Clist.form$nedges),
          as.integer(Clist.form$n),
          as.integer(Clist.form$dir), as.integer(Clist.form$bipartite),
          # Formation terms and proposals. 8
          as.integer(Clist.form$nterms), as.character(Clist.form$fnamestring), as.character(Clist.form$snamestring), as.integer(model.form$offset),
          as.character(MHproposal.form$name), as.character(MHproposal.form$package),
          as.double(Clist.form$inputs), eta.form=as.double(eta.form0),
          # Formation parameter fitting. 16
          nw.diff=as.double(control$nw.diff),
          as.integer(control$RM.phase2n),
          # Dissolution terms and proposals. 21
          as.integer(Clist.diss$nterms), as.character(Clist.diss$fnamestring), as.character(Clist.diss$snamestring),
          as.character(MHproposal.diss$name), as.character(MHproposal.diss$package),
          as.double(Clist.diss$inputs), as.double(eta.diss),
          # Degree bounds.
          as.integer(MHproposal.form$arguments$constraints$bd$attribs), 
          as.integer(MHproposal.form$arguments$constraints$bd$maxout), as.integer(MHproposal.form$arguments$constraints$bd$maxin),
          as.integer(MHproposal.form$arguments$constraints$bd$minout), as.integer(MHproposal.form$arguments$constraints$bd$minin),
          as.integer(MHproposal.form$arguments$constraints$bd$condAllDegExact), as.integer(length(MHproposal.form$arguments$constraints$bd$attribs)), 
          # MCMC settings.              
          as.integer(control$RM.burnin),
          as.integer(control$RM.interval),
          as.integer(control$MCMC.burnin),
          as.double(control$invGradient),
          # Space for output.
          as.integer(maxedges),
          newnwtails = integer(maxedges), newnwheads = integer(maxedges), 
          objective.history=double(length(eta.form0)*2*100000), #FIXME: Figure out how much space is needed.
          # Verbosity.
          as.integer(verbose), 
          PACKAGE="ergm") 

  eta.form <- z$eta.form
  names(eta.form) <- names(eta.form0)

  newnetwork<-newnw.extract(nw,z)
  
  oh <-  matrix(z$objective.history,ncol=length(eta.form0)*2,byrow=TRUE)
  
  list(nw.diff=z$nw.diff,
       newnetwork=newnetwork,
       coef.form=eta.form,
       objective.history=oh[apply(oh,1,any),])
}
