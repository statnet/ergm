stergm.EGMME.SA <- function(theta.form0, theta.diss0, nw, model.form, model.diss, model.mon,
                            control, MHproposal.form, MHproposal.diss, eval.optpars,
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
    z<-stergm.EGMME.SA.Phase2.C(state, model.form, model.diss, model.mon, MHproposal.form, MHproposal.diss,
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

    inds.last <- nrow(oh.all) + 1 - control$SA.runlength:1
    inds.keep <- nrow(oh.all) + 1 - max(nrow(oh.all)*control$SA.keep.oh,min(control$SA.runlength*4,nrow(oh.all))):1

    # Extract and store subhistories of interest.
    assign("oh",oh.all[inds.keep,,drop=FALSE],envir=parent.frame())
    assign("oh.last",oh.all[inds.last,,drop=FALSE],envir=parent.frame())
    assign("jitters",jitters.all[inds.keep,,drop=FALSE],envir=parent.frame())
    assign("jitters.last",jitters.all[inds.last,,drop=FALSE],envir=parent.frame())

    # Plot if requested.
    if(control$SA.plot.progress){
      if(!dev.interactive(TRUE)){
        warning("Progress plot requested on a non-interactive graphics device. Ignoring.")
        control$SA.plot.progress <- FALSE # So that we don't print a warning every step.
      }else{
        library(lattice)
        thin <- (nrow(oh)-1)%/%control$SA.max.plot.points + 1
        print(xyplot(mcmc(oh),panel=function(...){panel.xyplot(...);panel.abline(0,0)},thin=thin,as.table=TRUE))
      }
    }
    
    # Extract and return the "state".
    list(nw = z$newnetwork,
         nw.diff = z$nw.diff,
         eta.form = z$eta.form,
         eta.diss = z$eta.diss,
         oh.last = oh.last)
  }

  backoff.test <- function(state, state.bak, threshold, print.full=FALSE){
    if(is.null(state$oh.last) || is.null(state.bak$oh.last)) return(state)
    
    stats.last <- state$oh.last
    stats.last3 <- state.bak$oh.last

    w <- robust.inverse((cov(stats.last)/mean(diag(cov(stats.last)))+cov(stats.last3)/mean(diag(cov(stats.last3))))/2)
    
    m.last <- mahalanobis(stats.last, 0, w, inverted=TRUE)
    m.last3 <- mahalanobis(stats.last3, 0, w, inverted=TRUE)

    if(verbose) cat("Scaled objective function",mean(m.last3),"->",mean(m.last),"\n")
    if(mean(m.last)/mean(m.last3)>threshold^(if(is.null(state.bak$backoffs)) 1 else state.bak$backoffs)){
      if(print.full || verbose) cat("Objective function got significantly worse. Backing off.\n") else cat("!")
      state.bak$backoffs <- if(is.null(state.bak$backoffs)) 1 else state.bak$backoffs+1
      state.bak
    }else state
  }

  accel.brake <- function(state,state.bak){
    
    
    ph <- oh[-nrow(oh),1:p,drop=FALSE][,!offsets,drop=FALSE]
    ph.jumps <- ph[seq(from=1+control$SA.runlength,to=nrow(ph),by=control$SA.runlength),,drop=FALSE]-ph[seq(from=1,to=nrow(oh)-control$SA.runlength,by=control$SA.runlength),,drop=FALSE]

    last.jump <- c(state$eta.form,state$eta.diss)[!offsets] - c(state.bak$eta.form,state.bak$eta.diss)[!offsets]

    if(nrow(ph.jumps)<3) return(state)

    jump.w<-robust.inverse(cov(ph.jumps)/sum(diag(cov(ph.jumps))))
    
    jump.sizes <- mahalanobis(ph.jumps, 0, jump.w, inverted=TRUE)
    last.jump.size <- mahalanobis(rbind(last.jump), 0, jump.w, inverted=TRUE)

    if(verbose) cat("Most recent jump size relative to average is",last.jump.size/mean(jump.sizes),", with average being",mean(jump.sizes),"\n")
    
    if(last.jump.size>mean(jump.sizes)*control$SA.phase2.maxreljump){
      if(verbose) cat("Jumps accelerating too fast. Truncating.\n") else cat("?")
      par <- c(state.bak$eta.form,state.bak$eta.diss)[!offsets] + last.jump/(last.jump.size/(mean(jump.sizes)*control$SA.phase2.maxreljump))
      if(p.form.free) state$eta.form[!model.form$offset] <- par[seq_len(p.form.free)]
      if(p.diss.free) state$eta.diss[!model.diss$offset] <- par[p.form.free+seq_len(p.diss.free)]
    }
    state
  }
  
  interpolate.par <- function(h.fit, w=diag(1,nrow=ncol(h.fit))){
    x <- t(h.fit[-1,,drop=FALSE])
    y <- -cbind(h.fit[1,])

    c(solve(t(x)%*%w%*%x)%*%t(x)%*%w%*%y)
  }

  best.state <- function(){
    w <- robust.inverse(cov(oh.all[,-(1:p),drop=FALSE]))
    best.i <- which.min(mahalanobis((oh.all[-1,-(1:p),drop=FALSE]+oh.all[-nrow(oh.all),-(1:p),drop=FALSE])/2,0,w,inverted=TRUE))
    best.par <- oh.all[best.i,1:p][!offsets]
    if(p.form.free) state$eta.form[!model.form$offset] <- best.par[seq_len(p.form.free)]
    if(p.diss.free) state$eta.diss[!model.diss$offset] <- best.par[p.form.free+seq_len(p.diss.free)]
    state
  }


  
  ##### Construct the initial state. ######

  oh.all <- NULL
  jitters.all <- NULL
  state <- NULL

  for(restart in 1:control$SA.restarts){
  
    if(is.null(nw %n% "lasttoggle")) nw %n% "lasttoggle" <- rep(0, network.dyadcount(nw))
    if(is.null(nw %n% "time")) nw %n% "time" <- 0
    
    nw.diff <- model.mon$nw.stats - model.mon$target.stats  # nw.diff keeps track of the difference between the current network and the target statistics.
    
    state <- list(nw=nw,
                  eta.form = ergm.eta(theta.form0, model.form$etamap),
                  eta.diss = ergm.eta(theta.diss0, model.diss$etamap),
                  nw.diff  = nw.diff)
    
    cat('========  Phase 1: Burn in, get initial gradient values, and find a configuration under which all targets vary. ========\n',sep="")
    
    ###### Set up and run the burn-in. ######
    
    control$collect.form <- control$collect.diss <- FALSE
    control.phase1<-control
    control.phase1$time.samplesize <- 1
    control.phase1$time.burnin <- control$SA.burnin
    control.phase1$time.interval <- 1
    
    cat("Burning in... ")
    
    z <- stergm.getMCMCsample(state$nw, model.form, model.diss, model.mon, MHproposal.form, MHproposal.diss, state$eta.form, state$eta.diss, control.phase1, verbose)
    
    cat("Done.\n")
    
    # Update the state with burn-in results.
    state$nw <- z$newnetwork
    state$nw.diff <- state$nw.diff + z$statsmatrix.mon[NROW(z$statsmatrix.mon),]
    
    # Set backup states
    state.bak <- state
    
    ###### Gradient estimation and getting to a decent configuration. ######
    
    ## Here, "decent configuration" is defined as one where all target
    ## statistics have some variability --- none are stuck.
    
    # Start with jitter-only.
    
    control$GainM <- matrix(0, nrow=p, ncol=q)
    control$dejitter <- matrix(0, nrow=p, ncol=p) # Dejitter tries to cancel the effect of jitter on the optimizer.
    
    control$jitter<-rep(0,p)
    control$jitter[!offsets]<- control$SA.phase1.jitter
    
    ## Adjust the number of time steps between jumps using burn-in.
    edge.ages <- state$nw%n%"time"-ergm.el.lasttoggle(state$nw)[,3]+1
    
    control$SA.interval<- max(control$SA.min.interval, if(length(edge.ages)>0) control$SA.interval.mul*mean(edge.ages))
    if(verbose>1){
      cat("New interval:",control$SA.interval ,"\n")
    }  
    
    for(try in 1:control$SA.phase1.tries){
      if(verbose) cat('======== Attempt ',try,' ========\n',sep="") else cat('Attempt',try,':\n')
      state <- try(do.optimization(state, control), silent=!verbose)
      if(inherits(state, "try-error")){
        cat("Something went very wrong. Restarting with smaller gain.\n")
        control$SA.init.gain <- control$SA.init.gain / control$SA.gain.decay
        do.restart <- TRUE
        break
      }else do.restart <- FALSE

      #state <- backoff.test(state, state.bak, control$SA.phase1.backoff.rat, TRUE)
      state <- best.state()
      
      control$gain <- control$SA.init.gain
      out <- try(eval.optpars(TRUE,restart>1,FALSE), silent=!verbose)
      if(inherits(out, "try-error")){
        cat("Something went very wrong. Restarting with smaller gain.\n")
        control$SA.init.gain <- control$SA.init.gain / control$SA.gain.decay
        do.restart <- TRUE
        break
      }else do.restart <- FALSE
      control <- out$control

      if(all(!out$ineffectual.pars) && all(!out$bad.fits)){
        cat("Gradients estimated and all statistics are moving. Proceeding to Phase 2.\n")
        break
      }
      if(try==control$SA.phase1.tries) stop("The optimizer was unable to find a reasonable configuration: one or more statistics are still stuck after multiple tries, and one or more parameters do not appear to have any robust effect.")
    }
    if(do.restart) next
    
    ###### Main optimization run. ######
    
    cat('========  Phase 2: Find and refine the estimate. ========\n',sep="")
    
    for(subphase in 1:control$SA.phase2.levels){
      if(verbose) cat('======== Subphase ',subphase,' ========\n',sep="") else cat('Subphase 2.',subphase,' ',sep="")
      
      control$gain <- control$SA.init.gain*control$SA.gain.decay^(subphase-1)
      stepdown.count <- control$SA.stepdown.ct.base + round(control$SA.stepdown.ct.subphase*subphase)
      
      for(regain in 1:control$SA.phase2.repeats){
        if(verbose==0) cat(".")
        state <- try(do.optimization(state, control), silent=!verbose)
        if(inherits(state, "try-error")){
          cat("Something went very wrong. Restarting with smaller gain.\n")
          control$SA.init.gain <- control$SA.init.gain / control$SA.gain.decay
          do.restart <- TRUE
          break
        }else do.restart <- FALSE
        
        if(verbose){
          cat("New parameters:\n")
          cat("Formation:\n")
          print(state$eta.form)
          cat("Dissolution:\n")
          print(state$eta.diss)
        }
        
        # Perform backoff test and update backup
        state <- backoff.test(state, state.bak, control$SA.phase2.backoff.rat, FALSE)
        
        out <- try(eval.optpars(FALSE,TRUE,TRUE), silent=!verbose)
        if(inherits(out, "try-error")){
          cat("Something went very wrong. Restarting with smaller gain.\n")
          control$SA.init.gain <- control$SA.init.gain / control$SA.gain.decay
          do.restart <- TRUE
          break
        }else do.restart <- FALSE
        control <- out$control
        
        if(control$SA.phase2.refine){
          ## Interpolate the estimate.
          eta.int <- try(interpolate.par(out$oh.fit,out$w),silent=TRUE)
          if(!inherits(eta.int,"try-error")){
            if(p.form.free) state$eta.form[!model.form$offset] <- state$eta.form[!model.form$offset]*(1-min(control$gain*control$SA.runlength,1)) + eta.int[seq_len(p.form.free)]*min(control$gain*control$SA.runlength,1)
            if(p.diss.free) state$eta.diss[!model.diss$offset] <- state$eta.diss[!model.diss$offset]*(1-min(control$gain*control$SA.runlength,1)) + eta.int[p.form.free+seq_len(p.diss.free)]*min(control$gain*control$SA.runlength,1)
            if(verbose){
              cat("Interpolated parameters:\n")
              cat("Formation:\n")
              print(state$eta.form)
              cat("Dissolution:\n")
              print(state$eta.diss)
            }
          }
        }
        
        state <- accel.brake(state,state.bak)
        state.bak <- state
        
        ## If the optimization appears to be actively reducing the objective function, keep
        ## going, without reducing the gain.
        
        x <- unique(round(seq(from=1,to=NROW(oh),length.out=control$SA.stepdown.maxn)))
        ys <- oh[x,-(1:p),drop=FALSE]
        ys <- mahalanobis(ys,0,robust.inverse(cov(ys)),inverted=TRUE)
        
        fit <- try(summary(gls(ys~x,correlation=corAR1()))$tTable[2,c(1,4)])
        if(!inherits(fit, "try-error")){
          p.val <- fit[2]/2 # We are interested in one-sided decline here.
          est <- fit[1]
          if(est>0) p.val <- 1-p.val # If it's actually getting worse one-sided p-value is thus.
          if(verbose){
            cat("Trend in the objective function p-value:",p.val,". ")
          }
          if(p.val>control$SA.stepdown.p){
            stepdown.count <- stepdown.count - 1
            if(stepdown.count<=0){
              if(verbose) cat("No trend in objective function detected. Reducing gain.\n")
              stepdown.count <- control$SA.stepdown.ct.base + round((subphase+1)*control$SA.stepdown.ct.subphase)
              if(!verbose) cat("\n")
              break
            }else if(verbose) cat("No trend in objective function detected.",stepdown.count,"to go.\n")
          }else{
            stepdown.count <- control$SA.stepdown.ct.base + round(subphase*control$SA.stepdown.ct.subphase)
            if(verbose) cat("Trend in objective function detected. Resetting counter.\n")
          }
        }else{
          if(verbose) cat("Problem testing trend in objective function detected. Continuing with current gain.\n")
        }
        if(do.restart) break
      }

      mc.se <- {
        h <- oh[,1:p,drop=FALSE][,!offsets,drop=FALSE]
        apply(h,2,sd)/sqrt(effectiveSize(h))
      }

      if(verbose){
        cat("Approximate precision of the estimate:\n")
        print(mc.se)
      }

      if(all(mc.se < control$SA.phase2.max.mc.se)){
        if(verbose) cat("EGMME appears to be estimated to the desired precision level. Stopping.\n")
        break
      }
      
      if(do.restart) break      
    }
    if(!do.restart) break # If We've gotten this far, no restart conditions have been triggered, so we are good.
  }
  
  ## Refine the estimate.

  if(inherits(state, "try-error")) stop("Something went wrong too many times. Try better starting values or reducing control$SA.init.gain.") 
  
  eta.form <- state$eta.form
  eta.diss <- state$eta.diss

  eta.free <- switch(control$SA.refine,
                     mean = colMeans(oh[,1:p,drop=FALSE][,!offsets,drop=FALSE]),
                     linear =
                       if(p.free==q) solve(t(out$oh.fit[-1,,drop=FALSE]),-out$oh.fit[1,])
                       else interpolate.par(out$oh.fit,out$w),
                     none = c(eta.form,eta.diss)[!offsets])
  if(p.form.free) eta.form[!model.form$offset] <- eta.free[seq_len(p.form.free)]
  if(p.diss.free) eta.diss[!model.diss$offset] <- eta.free[p.form.free+seq_len(p.diss.free)]

  if(verbose){
    cat("Refining the estimate using the", control$SA.refine,"method. New estimate:\n")
    cat("Formation:\n")
    print(eta.form)
    cat("Dissolution:\n")
    print(eta.diss)
  }


  if(control$SA.se){
    G <- out$G
    w <- out$w
    
    control.phase3<-control
    control.phase3$time.burnin <- control$SA.burnin
    control.phase3$time.samplesize <- control$SA.phase3.samplesize
    control.phase3$time.interval <- control$SA.interval
    
    ## Estimate standard errors.
    cat('========  Phase 3: Simulate from the fit and estimate standard errors. ========\n',sep="")
    z <- stergm.getMCMCsample(state$nw, model.form, model.diss, model.mon, MHproposal.form, MHproposal.diss, eta.form, eta.diss, control.phase3, verbose)
    sm.mon <- sweep(z$statsmatrix.mon,2,state$nw.diff,"+")
    if(verbose)cat("Finished.\n")
    V.stat<-cov(sm.mon)
    V.par<-matrix(NA,p,p)
    V.par[!offsets,!offsets]<-solve(t(G)%*%w%*%G)%*%t(G)%*%w%*%V.stat%*%w%*%G%*%solve(t(G)%*%w%*%G)
    sm.mon <- mcmc(sm.mon, start=control$SA.burnin, thin=control$SA.interval)
  }else{
    V.par <- NULL
    sm.mon <- NULL
  }
  
  #ve<-with(z,list(coef=eta,sample=s$statsmatrix.form,sample.obs=NULL))
  
  #endrun <- control$MCMC.burnin+control$MCMC.interval*(ve$samplesize-1)
  #attr(ve$sample, "mcpar") <- c(control$MCMC.burnin+1, endrun, control$MCMC.interval)
  #attr(ve$sample, "class") <- "mcmc"
  
  list(newnetwork=if(control$SA.se) z$newnetwork else state$nw,
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

stergm.EGMME.SA.Phase2.C <- function(state, model.form, model.diss, model.mon,
                             MHproposal.form, MHproposal.diss, control, verbose) {
  Clist.form <- ergm.Cprepare(state$nw, model.form)
  Clist.diss <- ergm.Cprepare(state$nw, model.diss)
  Clist.mon <- ergm.Cprepare(state$nw, model.mon)
  maxedges <- max(control$MCMC.init.maxedges, Clist.mon$nedges)
  maxchanges <- max(control$MCMC.init.maxchanges, Clist.mon$nedges)

  repeat{
    z <- .C("MCMCDynSArun_wrapper",
            # Observed/starting network. 
            as.integer(Clist.form$tails), as.integer(Clist.form$heads),
            time = if(is.null(Clist.form$time)) as.integer(0) else as.integer(Clist.form$time),
            lasttoggle = if(is.null(Clist.form$time)) integer(network.dyadcount(state$nw)) else as.integer(Clist.form$lasttoggle),
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
            as.integer(control$SA.runlength),
            as.double(control$GainM),
            as.double(control$jitter), as.double(control$dejitter), # Add a little bit of noise to parameter guesses.
            # Degree bounds.
            as.integer(MHproposal.form$arguments$constraints$bd$attribs), 
            as.integer(MHproposal.form$arguments$constraints$bd$maxout), as.integer(MHproposal.form$arguments$constraints$bd$maxin),
            as.integer(MHproposal.form$arguments$constraints$bd$minout), as.integer(MHproposal.form$arguments$constraints$bd$minin),
            as.integer(MHproposal.form$arguments$constraints$bd$condAllDegExact), as.integer(length(MHproposal.form$arguments$constraints$bd$attribs)), 
            # MCMC settings.              
            as.integer(control$SA.burnin),
            as.integer(control$SA.interval),
            as.integer(control$MCMC.burnin),
            # Space for output.
            as.integer(maxedges),
            as.integer(maxchanges),
            newnwtails = integer(maxedges), newnwheads = integer(maxedges), 
            opt.history=double(((Clist.form$nstats+Clist.diss$nstats)*2+Clist.mon$nstats)*control$SA.runlength),
            # Verbosity.
            as.integer(max(verbose-1,0)),
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
