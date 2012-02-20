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

stergm.RM <- function(theta.form0, theta.diss0, nw, model.form, model.diss, model.mon,
                            control, MHproposal.form, MHproposal.diss,
                            verbose=FALSE){

  control$nw.diff <- model.mon$nw.stats - model.mon$target.stats # nw.diff keeps track of the difference between the current network and the target statistics.

  if(verbose) cat("Robbins-Monro algorithm with coef_F_0 = (",theta.form0, ") and coef_D = (",theta.diss0,")\n" )
  eta.form <- ergm.eta(theta.form0, model.form$etamap)
  eta.diss <- ergm.eta(theta.diss0, model.diss$etamap)

  offsets <- c(model.form$offset, model.diss$offset) # which parameters are offsets?
  p.form.free <- sum(!model.form$offset) # number of free formation parameters
  p.form <- length(model.form$offset) # total number of formation parameters
  p.diss.free <- sum(!model.diss$offset) # number of free dissolution parameters
  p.diss <- length(model.diss$offset) # total number of dissolution parameters
  p.free <- p.form.free+p.diss.free  # number of free parameters (formation and dissolution)
  p <- p.form+p.diss # total number of parameters (free and offset)
  
  q <- length(model.mon$offset) # number of target statistics/estimating equations --- may change when/if overspecified/weighted fitting is introduced

  control$collect.form <- control$collect.diss <- FALSE
  
  control.phase1<-control
  control.phase1$time.samplesize <- 1
  control.phase1$time.burnin <- control$RM.burnin
  control.phase1$time.interval <- 1

  # Run burn-in.
  z <- stergm.getMCMCsample(nw, model.form, model.diss, model.mon, MHproposal.form, MHproposal.diss, eta.form, eta.diss, control.phase1, verbose)

  nw <- z$newnetwork
  control$nw.diff <- control$nw.diff + z$statsmatrix.mon[NROW(z$statsmatrix.mon),]

  # Not assuming any relationship between parameters and statistics -> pure jittering at first.
  control$invGradient <- matrix(0, nrow=q, ncol=p)
  
  control$phase2n <- p+7+control$RM.phase2n_base

  control$jitter<-rep(0,p)

  oh <- NULL
  
  for(subphase in 1:control$RM.phase2sub){
    control$phase2n <- round(2.52*(control$phase2n-control$phase2n_base)+control$phase2n_base)
    control$jitter[!offsets]<- control$RM.init.jitter / subphase
    control$gain <- control$RM.init_gain / subphase
    
    
    for(retry in 1:control$RM.phase2sub_retries){
      
      z <- stergm.EGMoME.RM.Phase2.C(nw, model.form, model.diss, model.mon, MHproposal.form, MHproposal.diss,
                                     eta.form, eta.diss, control, verbose=verbose)
      nw <- z$newnetwork
      control$nw.diff<-z$nw.diff
      eta.form <- z$eta.form
      eta.diss <- z$eta.diss
      oh <- rbind(oh,z$opt.history) # Pool the history.

      # Figure out if any statistics are trapped.
      stats.min <- apply(oh[,-(1:p),drop=FALSE],2,min)
      stats.max <- apply(oh[,-(1:p),drop=FALSE],2,max)
      extreme.oh <- (sweep(oh[,-(1:p),drop=FALSE],2,stats.min,"-")==0) | (sweep(oh[,-(1:p),drop=FALSE],2,stats.max,"-")==0)
      oh.trustworthy <- apply(!extreme.oh,2,mean) > control$RM.prop.var.grad.OK

      # Regress statistics on parameters.
      # This uses GLS to account for serial correlation in statistics, and more recent are weighted higher.
      # First row is the intercept.

      oh.wt <- exp(control$RM.grad_decay*seq_len(NROW(oh)))
      x<-oh[,1:p,drop=FALSE] # #$%^$ gls() doesn't respect I()...
      x<-x[,!offsets,drop=FALSE]
      oh.fit <- sapply(1:q,
                       function(i){
                         y<-oh[,-(1:p),drop=FALSE][,i]
                         a<-try(coef(gls(y~x,subset=!extreme.oh[,i],weight=varFixed(~1/oh.wt),correlation=corAR1())))
                         if(inherits(a,"try-error")) rep(NA,p.free+1) else a
                       }
                       )
      
      if(all(is.na(oh.fit))) stop("The search is trapped. Try a different starting value, longer RM.interval, etc..")
      if(verbose){
        cat("Estimated gradient:\n")
        gr <- oh.fit[-1,]
        rownames(gr)<-c(paste("f.(",model.form$coef.names,")",sep=""),paste("d.(",model.diss$coef.names,")",sep=""))[!offsets]
        colnames(gr)<-model.mon$coef.names
        print(t(gr))
      }
      oh.fit[is.na(oh.fit)] <- 0
      oh.fit[,!oh.trustworthy] <- 0
      control$invGradient <- matrix(0, nrow=q, ncol=p)
      control$invGradient[,!offsets] <- robust.inverse(oh.fit[-1,]) * control$gain

      #control$jitter<-rep(0,p)
      if(!all(oh.trustworthy)){
        message("Only partial gradient matrix computed. Redoing the subphase.")
                    
      }else break
    }
  }

  if(control$RM.refine){
    eta.free <- -solve(t(oh.fit[-1,]),oh.fit[1,])
    if(p.form.free) eta.form[!model.form$offset] <- eta.free[seq_len(p.form.free)]
    if(p.diss.free) eta.diss[!model.diss$offset] <- eta.free[p.form.free+seq_len(p.diss.free)]
  }

  if(control$RM.se){
    control.phase3<-control
    control.phase3$time.burnin <- control$RM.burnin
    control.phase3$time.samplesize <- control$RM.phase3n*control$RM.interval
    control.phase3$time.interval <- 1
    
    # Run Phase 3.
    z <- stergm.getMCMCsample(nw, model.form, model.diss, MHproposal.form, MHproposal.diss, eta.form, eta.diss, control.phase3, verbose)
    G <- t(oh.fit[-1,])
    V.stat<-cov(z$statsmatrix.mon)
    V.par<-matrix(NA,p,p)
    V.par[!offsets,!offsets]<-solve(t(G)%*%G)%*%t(G)%*%V.stat%*%G%*%solve(t(G)%*%G)
  }else V.par <- NULL
  
  #ve<-with(z,list(coef=eta,sample=s$statsmatrix.form,sample.obs=NULL))
  names(eta.form)<-model.form$coef.names
  names(eta.diss)<-model.diss$coef.names
  
  #endrun <- control$MCMC.burnin+control$MCMC.interval*(ve$samplesize-1)
  #attr(ve$sample, "mcpar") <- c(control$MCMC.burnin+1, endrun, control$MCMC.interval)
  #attr(ve$sample, "class") <- "mcmc"
  
  list(newnetwork=nw, 
       init.form=theta.form0,
       init.diss=theta.diss0,
       covar=V.par,
       covar.form=V.par[seq_len(p.form),seq_len(p.form)],
       covar.diss=V.par[p.form+seq_len(p.diss),p.form+seq_len(p.diss)],
       eta.form=eta.form,
       eta.diss=eta.diss,
       opt.history=oh,
       sample=z$statsmatrix.mon,
       network=nw)            
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

stergm.EGMoME.RM.Phase2.C <- function(nw, model.form, model.diss, model.mon,
                             MHproposal.form, MHproposal.diss, eta.form0, eta.diss0,
                             control, verbose) {
  Clist.form <- ergm.Cprepare(nw, model.form)
  Clist.diss <- ergm.Cprepare(nw, model.diss)
  Clist.mon <- ergm.Cprepare(nw, model.mon)
  maxedges <- max(control$MCMC.init.maxedges, Clist.form$nedges)
  if(verbose){cat(paste("MCMCDyn workspace is",maxedges,"\n"))}
  
  z <- .C("MCMCDynRMPhase2_wrapper",
          # Observed/starting network. 
          as.integer(Clist.form$tails), as.integer(Clist.form$heads), 
          as.integer(Clist.form$nedges),
          as.integer(Clist.form$n),
          as.integer(Clist.form$dir), as.integer(Clist.form$bipartite),
          # Formation terms and proposals. 
          as.integer(Clist.form$nterms), as.character(Clist.form$fnamestring), as.character(Clist.form$snamestring),
          as.character(MHproposal.form$name), as.character(MHproposal.form$package),
          as.double(Clist.form$inputs), eta.form=as.double(eta.form0),
          # Dissolution terms and proposals. 
          as.integer(Clist.diss$nterms), as.character(Clist.diss$fnamestring), as.character(Clist.diss$snamestring),
          as.character(MHproposal.diss$name), as.character(MHproposal.diss$package),
          as.double(Clist.diss$inputs), eta.diss=as.double(eta.diss0),
          # Parameter fitting.
          as.integer(Clist.mon$nterms), as.character(Clist.mon$fnamestring), as.character(Clist.mon$snamestring),
          as.double(Clist.mon$inputs), 
          nw.diff=as.double(control$nw.diff),
          as.integer(control$RM.phase2n),
          as.double(control$invGradient),
          as.double(control$jitter), # Add a little bit of noise to parameter guesses.
          # Degree bounds.
          as.integer(MHproposal.form$arguments$constraints$bd$attribs), 
          as.integer(MHproposal.form$arguments$constraints$bd$maxout), as.integer(MHproposal.form$arguments$constraints$bd$maxin),
          as.integer(MHproposal.form$arguments$constraints$bd$minout), as.integer(MHproposal.form$arguments$constraints$bd$minin),
          as.integer(MHproposal.form$arguments$constraints$bd$condAllDegExact), as.integer(length(MHproposal.form$arguments$constraints$bd$attribs)), 
          # MCMC settings.              
          as.integer(control$RM.burnin),
          as.integer(control$RM.interval),
          as.integer(control$MCMC.burnin),
          # Space for output.
          as.integer(maxedges),
          newnwtails = integer(maxedges), newnwheads = integer(maxedges), 
          opt.history=double((Clist.form$nstats+Clist.diss$nstats+Clist.mon$nstats)*control$RM.phase2n),
          # Verbosity.
          as.integer(verbose), 
          PACKAGE="ergm") 

  eta.form <- z$eta.form
  names(eta.form) <- names(eta.form0)
  eta.diss <- z$eta.diss
  names(eta.diss) <- names(eta.diss0)

  newnetwork<-newnw.extract(nw,z)
  
  list(nw.diff=z$nw.diff,
       newnetwork=newnetwork,
       eta.form=eta.form,
       eta.diss=eta.diss,
       opt.history=matrix(z$opt.history,ncol=Clist.form$nstats+Clist.diss$nstats+Clist.mon$nstats,byrow=TRUE))
}
