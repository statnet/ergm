stergm.EGMME.GD <- function(theta.form0, theta.diss0, nw, model.form, model.diss, model.mon,
                            control, MHproposal.form, MHproposal.diss,
                            verbose=FALSE){

  eval.optpars <- function(test.G,window,update.jitter){
    for(name in ls(pos=parent.frame())) assign(name, get(name, parent.frame()))
    
    ## Regress statistics on parameters.
    # This uses GLS to account for serial correlation in statistics,
    # since we want p-values. First row is the intercept.

    h <- get(if(window) "oh" else "oh.all")
    control <- get("control")
    state <- get("state")
    
    x<-h[,1:p,drop=FALSE][,!offsets,drop=FALSE] # #$%^$ gls() doesn't respect I()...
    ys <- h[,-(1:p),drop=FALSE]

    h.fits <-
      if(test.G)
        sapply(1:q,
               function(i){
                 y<-ys[,i]
                 try(gls(y~x, correlation=corARMA(p=2)),silent=TRUE)
               },simplify=FALSE)
      else
        sapply(1:q,
             function(i){
               y<-ys[,i]
               suppressWarnings(try(lmrob(y~x), silent=TRUE))
             },simplify=FALSE)
    
    bad.fits <- sapply(h.fits, inherits, "try-error")
    bad.fits <-     # Also, ignore fits where the statistics are too concentrated.    
      bad.fits | (apply(ys,2,function(y){
        freqs <- table(y)
        sum(freqs[-which.max(freqs)])
      })<nrow(h)/2)

    if(all(bad.fits)) stop("The optimization appears to be stuck. Try better starting parameters, lower SA.init.gain, etc.")
    
    ## Grab the coefficients, t-values, and residuals.
    
    h.fit <- h.pvals <- matrix(NA, nrow=p.free+1,ncol=q)
    
    h.pvals[,!bad.fits] <- if(test.G) sapply(h.fits[!bad.fits],function(fit) summary(fit)$tTable[,4]) else 0
    h.fit[,!bad.fits] <- sapply(h.fits[!bad.fits], coef)      

    h.resid <- matrix(NA, nrow=nrow(h), ncol=q)
    h.resid[,!bad.fits] <- sapply(h.fits[!bad.fits], resid)
    
    G.signif <- t(h.pvals[-1,,drop=FALSE] < 1-(1-control$SA.phase1.max.p)^(p*q))
    G.signif[is.na(G.signif)] <- FALSE

    ## Compute the variances and the statistic weights.
    
    v <- matrix(NA, q,q)
    v[!bad.fits,!bad.fits] <- crossprod(h.resid[,!bad.fits,drop=FALSE])/nrow(h.resid)
    v[is.na(v)] <- 0

    w <- robust.inverse(v)
    
    ## Adjust the number of time steps between jumps.
    edge.ages <- state$nw%n%"time"-ergm.el.lasttoggle(state$nw)[,3]+1
    
    control$SA.interval<- max(control$SA.min.interval, if(length(edge.ages)>0) control$SA.interval.mul*mean(edge.ages))
    if(verbose>1){
      cat("New interval:",control$SA.interval ,"\n")
    }

    
    ## Detect parameters whose effect we aren't able to reliably detect.
    ineffectual.pars <- !apply(G.signif,2,any)

    if(all(ineffectual.pars)){
      cat("None of the parameters have a detectable effect. Increasing jitter.\n" )
      control$jitter[!offsets] <- control$jitter[!offsets]*2
    }
    else if(any(ineffectual.pars)){
      cat("Parameters",paste.and(p.names[!offsets][ineffectual.pars]),"do not have a detectable effect. Shifting jitter to them.\n" )
      control$jitter[!offsets] <- control$jitter[!offsets] * (ineffectual.pars+1/2) / mean(control$jitter[!offsets] * (ineffectual.pars+1/2))
    }

    ## Evaluate the dstat/dpar gradient matrix.
    G <- t(h.fit[-1,,drop=FALSE])
    G[!G.signif] <- 0
    G[is.na(G)] <- 0

    rownames(w)<-colnames(w)<-rownames(v)<-colnames(v)<-q.names
    
    colnames(G)<-p.names[!offsets]
    rownames(G)<-q.names
    if(verbose>1){
      cat("Most recent parameters:\n")
      cat("Formation:\n")
      print(state$eta.form)
      cat("Dissolution:\n")
      print(state$eta.diss)
      cat("Target differences (most recent):\n")
      print(state$nw.diff)
      cat("Target differences (last run):\n")
      print(colMeans(oh.last[,-(1:p),drop=FALSE]))
      cat("Approximate objective function (most recent):\n")
      print(mahalanobis(oh[nrow(oh),-(1:p),drop=FALSE],0,cov=w,inverted=TRUE))
      cat("Approximate objective function (last run):\n")
      print(mahalanobis(colMeans(oh.last[,-(1:p),drop=FALSE]),0,cov=w,inverted=TRUE))
      cat("Estimated gradient:\n")
      print(G)
      cat("Estimated covariance of statistics:\n")
      print(v)
    }

    
    control$GainM <- matrix(0, nrow=p, ncol=q)
    control$GainM[!offsets,] <- t(G) %*% w  
    wscale <- 1 / sqrt(mean(diag(control$GainM[!offsets,,drop=FALSE]%*%v%*%t(control$GainM[!offsets,,drop=FALSE])))) * min(1,sqrt(mahalanobis(colMeans(oh.last[,-(1:p),drop=FALSE]),0,cov=w,inverted=TRUE)))
    #wscale <-  #/ sqrt(mahalanobis(colMeans(oh.last[,-(1:p),drop=FALSE]),0,cov=w,inverted=TRUE))
    control$GainM[!offsets,] <- control$GainM[!offsets,] * wscale * control$gain
    control$GainM[!is.finite(control$GainM)] <- 0

    control$dejitter <- matrix(0, nrow=p, ncol=p)
    control$dejitter[!offsets,!offsets] <- control$GainM[!offsets,,drop=FALSE]%*%G
    
    rownames(control$GainM) <- rownames(control$dejitter) <- colnames(control$dejitter) <- p.names
    colnames(control$GainM) <- q.names
    
    if(verbose>1){
      cat("New deviation -> coefficient map:\n")
      print(control$GainM)
      cat("New jitter cancelation matrix:\n")
      print(control$dejitter)
    }

    if(update.jitter){
      control$jitter[!offsets] <- apply(oh[,1:p,drop=FALSE][,!offsets,drop=FALSE]-jitters[,!offsets,drop=FALSE],2,sd)*control$SA.phase2.jitter.mul
      names(control$jitter) <- p.names
    }
    
    if(verbose>1){
      cat("New jitter values:\n")
      print(control$jitter)
    }

    list(control=control,
         G=G, w=w, v=v, oh.fit=h.fit, ineffectual.pars=ineffectual.pars, bad.fits=bad.fits)
  }

    stergm.EGMME.SA(theta.form0, theta.diss0, nw, model.form, model.diss, model.mon,
                            control, MHproposal.form, MHproposal.diss, eval.optpars,
                            verbose)
}
