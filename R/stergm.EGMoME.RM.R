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

  if(verbose) cat("Starting optimization with with coef_F_0 = (",theta.form0, ") and coef_D_0 = (",theta.diss0,")\n" )
  
  ###### Set the constants and convenience variables. ######
  offsets <- c(model.form$offset, model.diss$offset) # which parameters are offsets?
  p.form.free <- sum(!model.form$offset) # number of free formation parameters
  p.form <- length(model.form$offset) # total number of formation parameters
  p.diss.free <- sum(!model.diss$offset) # number of free dissolution parameters
  p.diss <- length(model.diss$offset) # total number of dissolution parameters
  p.free <- p.form.free+p.diss.free  # number of free parameters (formation and dissolution)
  p <- p.form+p.diss # total number of parameters (free and offset)
  p.names<-c(paste("f.(",model.form$coef.names,")",sep=""),paste("d.(",model.diss$coef.names,")",sep=""))
  
  q <- length(model.mon$offset) # number of target statistics
  q.names<-model.mon$coef.names

  ###### Define the optimization run function. ######
  
  do.optimization<-function(state, control){
    oh.all <- get("oh.all",parent.frame())
    jitters.all <- get("jitters.all",parent.frame())

    if(verbose) cat("Running stochastic optimization... ")
    z<-stergm.EGMoME.RM.Phase2.C(state, model.form, model.diss, model.mon, MHproposal.form, MHproposal.diss,
                                 control, verbose=verbose)
    if(verbose) cat("Finished. Extracting.\n")
    
    # Extract and store the history of jitters.
    jitters.all <- rbind(jitters.all,z$opt.history[,p+1:p,drop=FALSE])
    colnames(jitters.all) <- p.names
    assign("jitters.all",jitters.all,envir=parent.frame())
    z$opt.history <- z$opt.history[,-(p+1:p),drop=FALSE]

    # Extract and store histroy of trials.
    oh.all <- rbind(oh.all,z$opt.history)
    colnames(oh.all) <- c(p.names,q.names)
    assign("oh.all",oh.all,envir=parent.frame())

    
    inds.last <- nrow(oh.all) + 1 - control$RM.runlength:1
    inds.keep <- nrow(oh.all) + 1 - max(nrow(oh.all)*control$RM.keep.oh,min(control$RM.runlength*2,nrow(oh.all))):1

    # Extract and store subhistories of interest.
    assign("oh",oh.all[inds.keep,,drop=FALSE],envir=parent.frame())
    assign("oh.last",oh.all[inds.last,,drop=FALSE],envir=parent.frame())
    assign("jitters",jitters.all[inds.keep,,drop=FALSE],envir=parent.frame())
    assign("jitters.last",jitters.all[inds.last,,drop=FALSE],envir=parent.frame())

    # Plot if requested.
    if(control$RM.plot.progress){library(lattice); print(xyplot(mcmc(oh),as.table=TRUE))}
    
    # Extract and return the "state".
    list(nw = z$newnetwork,
         nw.diff = z$nw.diff,
         eta.form = z$eta.form,
         eta.diss = z$eta.diss)
  }

  eval.optpars <- function(test.G,window,update.jitter){

     ## Regress statistics on parameters.
    # This uses GLS to account for serial correlation in statistics,
    # since we want p-values. First row is the intercept.

    h <- get(if(window) "oh" else "oh.all",parent.frame())
    control <- get("control",parent.frame())
    state <- get("state",parent.frame())
    
    x<-h[,1:p,drop=FALSE][,!offsets,drop=FALSE] # #$%^$ gls() doesn't respect I()...
    ys <- h[,-(1:p),drop=FALSE]

    h.fits <-
      if(test.G)
        sapply(1:q,
               function(i){
                 y<-ys[,i]
                 try(gls(y~x, correlation=corAR1()),silent=TRUE)
               },simplify=FALSE)
      else
        sapply(1:q,
             function(i){
               y<-ys[,i]
               try(lm(y~x))
             },simplify=FALSE)
    
    bad.fits <- sapply(h.fits, inherits, "try-error")

    ## Grab the coefficients, t-values, and residuals.
    
    h.fit <- h.pvals <- matrix(NA, nrow=p.free+1,ncol=q)
    
    h.pvals[,!bad.fits] <- if(test.G) sapply(h.fits[!bad.fits],function(fit) summary(fit)$tTable[,4]) else 0
    h.fit[,!bad.fits] <- sapply(h.fits[!bad.fits], coef)      

    h.resid <- matrix(NA, nrow=nrow(h), ncol=q)
    h.resid[,!bad.fits] <- sapply(h.fits[!bad.fits], resid)
    
    G.signif <- t(h.pvals[-1,,drop=FALSE] < 1-(1-control$RM.phase1.max.p)^(p*q))
    G.signif[is.na(G.signif)] <- FALSE

    ## Compute the variances and the statistic weights.
    
    v <- matrix(NA, q,q)
    v[!bad.fits,!bad.fits] <- crossprod(h.resid[,!bad.fits])/nrow(h.resid)
    v[is.na(v)] <- 0

    w <- robust.inverse(v)

    ## Adjust the number of time steps between jumps.
    
    stat.ac <- rep(NA,q)
    inds <- max(NROW(h)-control$RM.runlength+1,1):NROW(h)
    stat.ac[!bad.fits] <- diag(as.matrix(autocorr(mcmc(h.resid[inds,!bad.fits,drop=FALSE]),lags=1)[1,,]))
    if(verbose>1){
      cat("Lag-1 autocorrelation:\n")
      print(stat.ac)
    }
    control$RM.interval<-max(round(control$RM.interval*(1+sqrt(max(stat.ac,0,na.rm=TRUE))-sqrt(control$RM.target.ac))),control$RM.min.interval)
    if(verbose>1){
      cat("New interval:",control$RM.interval ,"\n")
    }

    
    ## Detect parameters whose effect we aren't able to reliably detect.
    ineffectual.pars <- !apply(G.signif,2,any)

    if(any(ineffectual.pars)){
      cat("Parameters ",paste(p.names[!offsets][ineffectual.pars],sep=", ")," do not have a detectable effect. Shifting jitter to them.\n" ,sep="")
      control$jitter[!offsets] <- control$jitter[!offsets] * (ineffectual.pars+1/2) / mean(control$jitter[!offsets] * (ineffectual.pars+1/2))
    }

    ## Evaluate the dstat/dpar gradient matrix.
    G <- t(h.fit[-1,,drop=FALSE])
    G[!G.signif] <- 0
    G[is.na(G)] <- 0

    rownames(v)<-colnames(v)<-rownames(w)<-colnames(w)<-q.names
    
    colnames(G)<-p.names[!offsets]
    rownames(G)<-q.names
    if(verbose>1){
      cat("Most recent parameters:\n")
      cat("Formation:\n")
      print(state$eta.form)
      cat("Dissolution:\n")
      print(state$eta.diss)
      cat("Target differences:\n")
      print(state$nw.diff)
      cat("Approximate objective function:\n")
      print(mahalanobis(oh[nrow(oh),-(1:p),drop=FALSE],0,cov=w,inverted=TRUE))
      cat("Estimated gradient:\n")
      print(G)
      cat("Estimated covariance of statistics:\n")
      print(v)
    }
          
    control$WinvGradient <- matrix(0, nrow=p, ncol=q)
    control$WinvGradient[!offsets,] <- robust.inverse(G) %*% chol(w,pivot=TRUE) * control$gain
    control$dejitter <- matrix(0, nrow=p, ncol=p)
    control$dejitter[!offsets,!offsets] <- control$WinvGradient[!offsets,]%*%G
    
    rownames(control$WinvGradient) <- rownames(control$dejitter) <- colnames(control$dejitter) <- p.names
    colnames(control$WinvGradient) <- q.names
    
    if(verbose>1){
      cat("New deviation -> coefficient map:\n")
      print(control$WinvGradient)
      cat("New jitter cancelation matrix:\n")
      print(control$dejitter)
    }

    if(update.jitter){
      control$jitter[!offsets] <- apply(oh[,1:p,drop=FALSE][,!offsets,drop=FALSE]-jitters[,!offsets,drop=FALSE],2,sd)*control$RM.jitter.mul
      names(control$jitter) <- p.names
    }
    
    if(verbose>1){
      cat("New jitter values:\n")
      print(control$jitter)
    }

    list(control=control,
         G=G, w=w, v=v, oh.fit=h.fit, ineffectual.pars=ineffectual.pars, bad.fits=bad.fits)
  }

  interpolate.par <- function(h.fit, w=diag(1,nrow=ncol(h.fit))){
    x <- t(h.fit[-1,,drop=FALSE])
    y <- -cbind(h.fit[1,])

    c(solve(t(x)%*%w%*%x)%*%t(x)%*%w%*%y)
  }


  
  ##### Construct the initial state. ######

  nw %n% "lasttoggle" <- rep(0, network.dyadcount(nw))
  nw %n% "time" <- 0
  
  state <- list(nw=nw,
                eta.form = ergm.eta(theta.form0, model.form$etamap),
                eta.diss = ergm.eta(theta.diss0, model.diss$etamap),
                nw.diff  = model.mon$nw.stats - model.mon$target.stats) # nw.diff keeps track of the difference between the current network and the target statistics.

  oh.all <- NULL
  jitters.all <- NULL

  ###### Set up and run the burn-in. ######
  
  control$collect.form <- control$collect.diss <- FALSE
  control.phase1<-control
  control.phase1$time.samplesize <- 1
  control.phase1$time.burnin <- control$RM.burnin
  control.phase1$time.interval <- 1
  
  z <- stergm.getMCMCsample(state$nw, model.form, model.diss, model.mon, MHproposal.form, MHproposal.diss, state$eta.form, state$eta.diss, control.phase1, verbose)

  # Update the state with burn-in results.
  state$nw <- z$newnetwork
  state$nw.diff <- state$nw.diff + z$statsmatrix.mon[NROW(z$statsmatrix.mon),]

  ###### Gradient estimation and getting to a decent configuration. ######

  ## Here, "decent configuration" is defined as one where all target
  ## statistics have some variability --- none are stuck.

  # Start with jitter-only.

  control$WinvGradient <- matrix(0, nrow=p, ncol=q)
  control$dejitter <- matrix(0, nrow=p, ncol=p) # Dejitter tries to cancel the effect of jitter on the optimizer.
  
  control$jitter<-rep(0,p)
  control$jitter[!offsets]<- control$RM.phase1.jitter

  control$RM.interval <- control$RM.init.interval
  
  if(verbose) cat('========  Phase 1: Get initial gradient values and find a configuration under which all targets vary. ========\n',sep="")
  
  for(try in 1:control$RM.phase1.tries){
    if(verbose) cat('======== Attempt ',try,' ========\n',sep="")
    state <- do.optimization(state, control)
    control$gain <- control$RM.init.gain
    out <- eval.optpars(TRUE,FALSE,FALSE)
    control <- out$control
    if(all(!out$ineffectual.pars) && all(!out$bad.fits)){
      if(verbose) cat("Gradients estimated and all statistics are moving. Moving on to Phase 2.\n")
      break
    }
  }

  ###### Main optimization run. ######
  
  for(subphase in 1:control$RM.phase2sub){
    if(verbose) cat('======== Subphase ',subphase,' ========\n',sep="")
    control$gain <- control$RM.init.gain*control$RM.gain.decay^(subphase-1)
    for(regain in 1:control$RM.phase2regain){
      state <-do.optimization(state, control)
      out <- eval.optpars(FALSE,TRUE,TRUE)
      control <- out$control

      if(control$RM.phase2.refine){
        ## Interpolate the estimate.
        eta.int <- try(interpolate.par(out$oh.fit,out$w),silent=TRUE)
        if(!inherits(eta.int,"try-error")){
          eta.last <- c(state$eta.form[!model.form$offset],state$eta.diss[!model.diss$offset])
          eta.int <- eta.last*(1-control$gain) + eta.int*control$gain

          eta.sds <- apply(oh[1:p,],2,sd)[!offsets]
          eta.zs <- (eta.int - eta.last)/eta.sds
          eta.zs <- pmax(pmin(eta.zs, control$RM.phase2.refine.maxjump), -control$RM.phase2.refine.maxjump)
          eta.int <- eta.zs*eta.sds + eta.last

          if(p.form.free) state$eta.form[!model.form$offset] <- eta.int[seq_len(p.form.free)]
          if(p.diss.free) state$eta.diss[!model.diss$offset] <- eta.int[p.form.free+seq_len(p.diss.free)]
          if(verbose){
            cat("Interpolated parameters:\n")
            cat("Formation:\n")
            print(state$eta.form)
            cat("Dissolution:\n")
            print(state$eta.diss)
          }
        }
      }

      ## If the optimization appears to be "going somewhere", keep
      ## going, without reducing the gain.

      ys <- oh[,(1:p),drop=FALSE][,!offsets,drop=FALSE]
      x <- 1:NROW(oh)
      
      eta.trend.ps <- apply(ys,2,
                            function(y)
                            summary(gls(y~x,correlation=corAR1(),
                                        subset=max(NROW(oh)-control$RM.stepdown.subphases*control$RM.runlength+1,1):NROW(oh)))$tTable[2,4])
      names(eta.trend.ps)<-c(paste("f.(",model.form$coef.names,")",sep=""),paste("d.(",model.diss$coef.names,")",sep=""))[!offsets]
      if(verbose>1){
        cat("Parameter trend p-values:\n")
        print(eta.trend.ps)
      }
      eta.trend.p <- pchisq(-2 * sum(log(eta.trend.ps)),length(eta.trend.ps)*2,lower.tail=FALSE)
      if(eta.trend.p>control$RM.stepdown.p){
        cat("No trend in parameters detected. Reducing gain.\n")
        break
      }else{
        cat("Trend in parameters detected. Continuing with current gain.\n")
      }
    }
  }

  ## Refine the estimate.
  
  eta.form <- state$eta.form
  eta.diss <- state$eta.diss

  eta.free <- switch(control$RM.refine,
                     mean = colMeans(oh[,1:p,drop=FALSE][,!offsets,drop=FALSE]),
                     linear =
                       if(p.free==q) solve(t(out$oh.fit[-1,,drop=FALSE]),-out$oh.fit[1,])
                       else interpolate.par(out$oh.fit,out$w),
                     none = c(eta.form,eta.diss)[!offsets])
  if(p.form.free) eta.form[!model.form$offset] <- eta.free[seq_len(p.form.free)]
  if(p.diss.free) eta.diss[!model.diss$offset] <- eta.free[p.form.free+seq_len(p.diss.free)]

  if(verbose){
    cat("Refining the estimate using the", control$RM.refine,"method. New estimate:\n")
    cat("Formation:\n")
    print(eta.form)
    cat("Dissolution:\n")
    print(eta.diss)
  }


  if(control$RM.se){
    G <- out$G
    w <- out$w
    
    control.phase3<-control
    control.phase3$time.burnin <- control$RM.burnin
    control.phase3$time.samplesize <- control$RM.phase3n
    control.phase3$time.interval <- control$RM.interval
    
    ## Estimate standard errors.
    if(verbose)cat("Simulating to estimate standard errors... ")
    z <- stergm.getMCMCsample(state$nw, model.form, model.diss, model.mon, MHproposal.form, MHproposal.diss, eta.form, eta.diss, control.phase3, verbose)
    sm.mon <- sweep(z$statsmatrix.mon,2,state$nw.diff,"+")
    if(verbose)cat("Finished.\n")
    V.stat<-cov(sm.mon)
    V.par<-matrix(NA,p,p)
    V.par[!offsets,!offsets]<-solve(t(G)%*%w%*%G)%*%t(G)%*%w%*%V.stat%*%w%*%G%*%solve(t(G)%*%w%*%G)
  }else{
    V.par <- NULL
    sm.mon <- NULL
  }
  
  #ve<-with(z,list(coef=eta,sample=s$statsmatrix.form,sample.obs=NULL))
  
  #endrun <- control$MCMC.burnin+control$MCMC.interval*(ve$samplesize-1)
  #attr(ve$sample, "mcpar") <- c(control$MCMC.burnin+1, endrun, control$MCMC.interval)
  #attr(ve$sample, "class") <- "mcmc"
  
  list(newnetwork=if(control$RM.se) z$newnetwork else state$nw,
       init.form=theta.form0,
       init.diss=theta.diss0,
       covar=V.par,
       covar.form=V.par[seq_len(p.form),seq_len(p.form),drop=FALSE],
       covar.diss=V.par[p.form+seq_len(p.diss),p.form+seq_len(p.diss),drop=FALSE],
       eta.form=eta.form,
       eta.diss=eta.diss,
       opt.history=oh.all,
       sample=sm.mon,
       network=nw)            
}

stergm.EGMoME.RM.Phase2.C <- function(state, model.form, model.diss, model.mon,
                             MHproposal.form, MHproposal.diss, control, verbose) {
  Clist.form <- ergm.Cprepare(state$nw, model.form)
  Clist.diss <- ergm.Cprepare(state$nw, model.diss)
  Clist.mon <- ergm.Cprepare(state$nw, model.mon)
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
            eta=as.double(c(state$eta.form,state$eta.diss)),
            as.integer(Clist.mon$nterms), as.character(Clist.mon$fnamestring), as.character(Clist.mon$snamestring),
            as.double(Clist.mon$inputs), 
            nw.diff=as.double(state$nw.diff),
            as.integer(control$RM.runlength),
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
            opt.history=double(((Clist.form$nstats+Clist.diss$nstats)*2+Clist.mon$nstats)*control$RM.runlength),
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
  names(eta.form) <- model.form$coef.names
  eta.diss <- z$eta[-seq_len(Clist.form$nstats)]
  names(eta.diss) <- model.diss$coef.names

  newnetwork<-newnw.extract(state$nw,z)
  newnetwork %n% "time" <- z$time
  newnetwork %n% "lasttoggle" <- z$lasttoggle
  
  list(nw.diff=z$nw.diff,
       newnetwork=newnetwork,
       eta.form=eta.form,
       eta.diss=eta.diss,
       opt.history=matrix(z$opt.history,ncol=(Clist.form$nstats+Clist.diss$nstats)*2+Clist.mon$nstats,byrow=TRUE))
}
