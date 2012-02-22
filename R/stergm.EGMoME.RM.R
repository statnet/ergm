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

  if(verbose) cat("Robbins-Monro algorithm with coef_F_0 = (",theta.form0, ") and coef_D_0 = (",theta.diss0,")\n" )
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
  
  nw %n% "lasttoggle" <- rep(0, network.dyadcount(nw))
  nw %n% "time" <- 0

  # Run burn-in.
  z <- stergm.getMCMCsample(nw, model.form, model.diss, model.mon, MHproposal.form, MHproposal.diss, eta.form, eta.diss, control.phase1, verbose)

  nw <- z$newnetwork
  control$nw.diff <- control$nw.diff + z$statsmatrix.mon[NROW(z$statsmatrix.mon),]

  # Not assuming any relationship between parameters and statistics -> pure jittering at first.
  control$WinvGradient <- matrix(0, nrow=p, ncol=q)
  control$dejitter <- matrix(0, nrow=p, ncol=p) # Dejitter tries to cancel the effect of jitter on the optimizer.
  
  control$jitter<-rep(0,p)
  control$jitter[!offsets]<- control$RM.init.jitter
  oh.all <- NULL
  jitters.all <- NULL
  
  for(subphase in 1:control$RM.phase2sub){
    if(verbose) cat('======== Subphase ',subphase,' ========\n',sep="")
    control$gain <- control$RM.init.gain / subphase

    for(regain in 1:control$RM.phase2regain){
    
    for(retry in 1:control$RM.phase2sub_retries){
      if(verbose) cat("Running stochastic optimization... ")
      z <- stergm.EGMoME.RM.Phase2.C(nw, model.form, model.diss, model.mon, MHproposal.form, MHproposal.diss,
                                     eta.form, eta.diss, control, verbose=verbose)
      if(verbose) cat("Finished.\n")

      # Extract the history of jitters.
      jitters.all <- rbind(jitters.all,z$opt.history[,p+1:p,drop=FALSE])
      z$opt.history <- z$opt.history[,-(p+1:p),drop=FALSE]
      
      oh.all <- rbind(oh.all,z$opt.history) # Pool the history.
      colnames(oh.all) <- c(paste("f.(",model.form$coef.names,")",sep=""),paste("d.(",model.diss$coef.names,")",sep=""),model.mon$coef.names)

      oh <- oh.all[nrow(oh.all)+1 - max(nrow(oh.all)*control$RM.keep.oh,min(control$RM.phase2n*2,nrow(oh.all))):1,,drop=FALSE]
      jitters <- jitters.all[nrow(jitters.all)+1 - max(nrow(jitters.all)*control$RM.keep.oh,min(control$RM.phase2n*2,nrow(jitters.all))):1,,drop=FALSE]
      
      if(control$RM.plot.progress){library(lattice); print(xyplot(mcmc(oh,start=nrow(oh.all)-nrow(oh)+1)))}

      # Figure out if any statistics are trapped.
      stats.min <- apply(oh.all[,-(1:p),drop=FALSE],2,min)
      stats.max <- apply(oh.all[,-(1:p),drop=FALSE],2,max)
      extreme.oh <- (sweep(oh[,-(1:p),drop=FALSE],2,stats.min,"-")==0) | (sweep(oh[,-(1:p),drop=FALSE],2,stats.max,"-")==0)
      oh.trustworthy <- apply(!extreme.oh,2,mean) > control$RM.prop.var.grad.OK

      if(mean(!apply(extreme.oh,1,any))>control$RM.prop.var.grad.OK.ignore){
        # If vast majority are OK, include everything.
        extreme.oh[]<-FALSE
        oh.trustworthy[]<-TRUE
      }
        

      # Regress statistics on parameters.
      # This uses GLS to account for serial correlation in statistics, and more recent are weighted higher.
      # First row is the intercept.

      oh.wt <- rep(1,NROW(oh))#seq_len(NROW(oh)) #exp(control$RM.grad_decay*seq_len(NROW(oh)))
      x<-oh[,1:p,drop=FALSE][,!offsets,drop=FALSE] # #$%^$ gls() doesn't respect I()...

      ys <- oh[,-(1:p),drop=FALSE]
      oh.fits <- sapply(1:q,
                       function(i){
                         y<-ys[,i]
                         #try(gls(y~x,subset=!extreme.oh[,i],weight=varFixed(~1/oh.wt),correlation=corAR1()))
                         try(lm(y~x,subset=!extreme.oh[,i],weight=oh.wt))
                       },simplify=FALSE)

      bad.fits <- sapply(oh.fits, inherits, "try-error")

      oh.fit <- t.vals <- matrix(NA, nrow=p.free+1,ncol=q)
      t.vals[,!bad.fits] <- sapply(oh.fits[!bad.fits], function(a) summary(a)$coefficients[,3]) # sapply(oh.fits[!bad.fits], function(a) summary(a)$tTable[,3])
      oh.fit[,!bad.fits] <- sapply(oh.fits[!bad.fits], coef)
      
      if(all(is.na(oh.fit))) stop("The search is trapped. Try a different starting value, longer RM.interval, etc..")
      oh.resid <- matrix(NA, nrow=NROW(oh),ncol=q)
      for(i in seq_along(oh.fits))
        if(!bad.fits[i])
          oh.resid[!extreme.oh[,i],i] <- resid(oh.fits[[i]])
      
      oh.r.na <- apply(oh.resid[,!bad.fits,drop=FALSE],1,function(x) any(is.na(x)))
      v <- crossprod(sweep(oh.resid[!oh.r.na,],1,sqrt(oh.wt[!oh.r.na]),"*"))/sum(oh.wt[!oh.r.na])
      v[is.na(v)] <- 0
      
      w <- if(q==1) 1/v else robust.inverse(v) #robust.inverse(diag(diag(v)))
      rownames(v)<-colnames(v)<-rownames(w)<-colnames(w)<-model.mon$coef.names
      
      G <- t(oh.fit[-1,,drop=FALSE])
      colnames(G)<-c(paste("f.(",model.form$coef.names,")",sep=""),paste("d.(",model.diss$coef.names,")",sep=""))[!offsets]
      rownames(G)<-model.mon$coef.names
      if(verbose){
        cat("Most recent parameters:\n")
        cat("Formation:\n")
        print(z$eta.form)
        cat("Dissolution:\n")
        print(z$eta.diss)
        cat("Target differences:\n")
        print(z$nw.diff)
        cat("Approximate objective function:\n")
        print(mahalanobis(oh[nrow(oh),-(1:p),drop=FALSE],0,cov=w,inverted=TRUE))
        cat("Estimated gradient:\n")
        print(G)
        cat("Estimated covariance of statistics:\n")
        print(v)
      }

      redo.subphase <- FALSE

      oh.obj <- mahalanobis(oh[,-(1:p),drop=FALSE], 0, w, inverted=TRUE)
      oh.obj.sub <- oh.obj[length(oh.obj) - control$RM.phase2n:1 + 1]
      oh.obj.prev <- oh.obj[-(length(oh.obj) - control$RM.phase2n:1 + 1)]
      if(length(oh.obj.prev) && t.test(oh.obj.sub,oh.obj.prev,alternative="greater")$p.value<control$RM.revert.threshold.p){
        message("Objective function did not improve.")
        redo.subphase <- TRUE
      }
      
      if(!all(oh.trustworthy)){
        message("Only partial gradient matrix computed reliably.")
        redo.subphase <- TRUE
      }
      
      G[is.na(G)] <- 0
      G[!oh.trustworthy,] <- 0
      # Better fail to update at all than to update in the wrong direction.
      G[abs(t(t.vals[-1,]))<control$RM.gls.min.t] <- 0

      ineffectual <- apply(G==0, 2, all)
      if(any(ineffectual)){
        message("Unable to reliably estimate the effect some parameters.")
        control$jitter[!offsets] <- control$jitter[!offsets]*(1+ineffectual)/mean(1+ineffectual)
        redo.subphase<-TRUE
      }

      control$WinvGradient <- matrix(0, nrow=p, ncol=q)
      control$WinvGradient[!offsets,] <- robust.inverse(G) %*% chol(w,pivot=TRUE) #%*%cw0
      control$WinvGradient[!offsets,] <- control$WinvGradient[!offsets,]/sum(abs(control$WinvGradient[!offsets,])) * control$gain
      control$dejitter <- matrix(0, nrow=p, ncol=p)
      control$dejitter[!offsets,!offsets] <- control$WinvGradient[!offsets,]%*%G
      
      rownames(control$WinvGradient) <- rownames(control$dejitter) <- colnames(control$dejitter) <- c(paste("f.(",model.form$coef.names,")",sep=""),paste("d.(",model.diss$coef.names,")",sep=""))
      colnames(control$WinvGradient) <- model.mon$coef.names

      if(verbose){
        cat("New deviation -> coefficient map:\n")
        print(control$WinvGradient)
        cat("New jitter cancelation matrix:\n")
        print(control$dejitter)
      }
     
      if(!redo.subphase){
        nw <- z$newnetwork
        control$nw.diff<-z$nw.diff
        eta.form <- z$eta.form
        eta.diss <- z$eta.diss
        names(eta.form)<-model.form$coef.names
        names(eta.diss)<-model.diss$coef.names
        break
      }else{message("Redoing subphase.")}
    }

    control$jitter[!offsets] <- apply(oh[,1:p,drop=FALSE][,!offsets,drop=FALSE]-jitters[,!offsets,drop=FALSE],2,sd)*control$RM.jitter.mul
    
    names(control$jitter) <- c(paste("f.(",model.form$coef.names,")",sep=""),paste("d.(",model.diss$coef.names,")",sep=""))
    
    if(verbose){
      cat("New jitter values:\n")
      print(control$jitter)
    }


    ys <- oh[,(1:p),drop=FALSE][,!offsets,drop=FALSE]
    x <- 1:NROW(oh)
    
    eta.trend.ps <- apply(ys,2,
                          function(y)
                          summary(gls(y~x,correlation=corAR1(),
                                      subset=max(NROW(oh)-control$RM.stepdown.phases*control$RM.phase2n+1,1):NROW(oh)))$tTable[2,4])
    cat("Parameter trend p-values:\n")
    names(eta.trend.ps)<-c(paste("f.(",model.form$coef.names,")",sep=""),paste("d.(",model.diss$coef.names,")",sep=""))[!offsets]
    print(eta.trend.ps)
    eta.trend.p <- pchisq(-2 * sum(log(eta.trend.ps)),length(eta.trend.ps)*2,lower.tail=FALSE)
    if(eta.trend.p>control$RM.stepdown.p){
      cat("No trend in parameters detected. Reducing gain.\n")
      break
    }else{
      cat("Trend in parameters detected. Continuing with current gain.\n")
    }
    
    
    
  }
    
  }

  if(control$RM.refine){
    last.inds <- nrow(oh) + 1 - control$RM.refine:1
    eta.free <- colSums(sweep(oh[last.inds,1:p,drop=FALSE][,!offsets,drop=FALSE],1,oh.wt[last.inds],"*"))/sum(oh.wt[last.inds])
    if(p.form.free) eta.form[!model.form$offset] <- eta.free[seq_len(p.form.free)]
    if(p.diss.free) eta.diss[!model.diss$offset] <- eta.free[p.form.free+seq_len(p.diss.free)]
    if(verbose){
      cat("Pooling optimization history to refine the estimate. New estimate:\n")
      cat("Formation:\n")
      print(eta.form)
      cat("Dissolution:\n")
      print(eta.diss)
    }
  }

  if(control$RM.se){
    control.phase3<-control
    control.phase3$time.burnin <- control$RM.burnin
    control.phase3$time.samplesize <- control$RM.phase3n*control$RM.interval
    control.phase3$time.interval <- 1
    
    # Run Phase 3.
    if(verbose)cat("Simulating to estimate standard errors... ")
    z <- stergm.getMCMCsample(nw, model.form, model.diss, model.mon, MHproposal.form, MHproposal.diss, eta.form, eta.diss, control.phase3, verbose)
    if(verbose)cat("Finished.\n")
    V.stat<-cov(z$statsmatrix.mon)
    V.par<-matrix(NA,p,p)
    V.par[!offsets,!offsets]<-solve(t(G)%*%w%*%G)%*%t(G)%*%w%*%V.stat%*%w%*%G%*%solve(t(G)%*%w%*%G)
  }else V.par <- NULL
  
  #ve<-with(z,list(coef=eta,sample=s$statsmatrix.form,sample.obs=NULL))
  
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
       opt.history=oh.all,
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
  maxedges <- max(control$MCMC.init.maxedges, Clist.mon$nedges)
  maxchanges <- max(control$MCMC.init.maxchanges, Clist.mon$nedges)

  repeat{
    z <- .C("MCMCDynRMPhase2_wrapper",
            # Observed/starting network. 
            as.integer(Clist.form$tails), as.integer(Clist.form$heads),
            time = if(is.null(Clist.form$time)) as.integer(0) else as.integer(Clist.form$time),
            lasttoggle = if(is.null(Clist.form$time)) integer(network.dyadcount(nw)) else as.integer(Clist.form$lasttoggle),
            as.integer(Clist.form$nedges),
            as.integer(Clist.form$n),
            as.integer(Clist.form$dir), as.integer(Clist.form$bipartite),
            # Formation terms and proposals. 
            as.integer(Clist.form$nterms), as.character(Clist.form$fnamestring), as.character(Clist.form$snamestring),
            as.character(MHproposal.form$name), as.character(MHproposal.form$package),
            as.double(Clist.form$inputs),
            # Dissolution terms and proposals. 
            as.integer(Clist.diss$nterms), as.character(Clist.diss$fnamestring), as.character(Clist.diss$snamestring),
            as.character(MHproposal.diss$name), as.character(MHproposal.diss$package),
            as.double(Clist.diss$inputs),
            # Parameter fitting.
            eta=as.double(c(eta.form0,eta.diss0)),
            as.integer(Clist.mon$nterms), as.character(Clist.mon$fnamestring), as.character(Clist.mon$snamestring),
            as.double(Clist.mon$inputs), 
            nw.diff=as.double(control$nw.diff),
            as.integer(control$RM.phase2n),
            as.double(control$WinvGradient),
            as.double(control$jitter), as.double(control$dejitter), # Add a little bit of noise to parameter guesses.
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
            as.integer(maxchanges),
            newnwtails = integer(maxedges), newnwheads = integer(maxedges), 
            opt.history=double(((Clist.form$nstats+Clist.diss$nstats)*2+Clist.mon$nstats)*control$RM.phase2n),
            # Verbosity.
            as.integer(verbose),
            status = integer(1), # 0 = OK, MCMCDyn_TOO_MANY_EDGES = 1, MCMCDyn_MH_FAILED = 2, MCMCDyn_TOO_MANY_CHANGES = 3
            PACKAGE="ergm")
    if(z$status==0) break;
    if(z$status==1){
      maxedges <- 5*maxedges
      message("Too many edges encountered in the simulation. Increasing capacity to ", maxedges)
    }
    if(z$status==3){
      maxchanges <- 5*maxchanges
      message("Too many changes elapsed in the simulation. Increasing capacity to ", maxchanges)
    }
  }

  eta.form <- z$eta[seq_len(Clist.form$nstats)]
  names(eta.form) <- names(eta.form0)
  eta.diss <- z$eta[-seq_len(Clist.form$nstats)]
  names(eta.diss) <- names(eta.diss0)

  newnetwork<-newnw.extract(nw,z)
  newnetwork %n% "time" <- z$time
  newnetwork %n% "lasttoggle" <- z$lasttoggle
  
  list(nw.diff=z$nw.diff,
       newnetwork=newnetwork,
       eta.form=eta.form,
       eta.diss=eta.diss,
       opt.history=matrix(z$opt.history,ncol=(Clist.form$nstats+Clist.diss$nstats)*2+Clist.mon$nstats,byrow=TRUE))
}
