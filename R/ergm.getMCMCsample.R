#  File ergm/R/ergm.getMCMCsample.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
#########################################################################################
# The <ergm.getMCMCsample> function samples networks using an MCMC algorithm via
# <MCMC_wrapper.C>. Unlike its <ergm.getMCMCsample> counterpart, this function is
# caple of running in multiple threads.  Note that the returned stats will be relative to
# the original network, i.e., the calling function must shift the statistics if required. 
# The calling function must also attach column names to the statistics matrix if required.
#########################################################################################

ergm.getMCMCsample <- function(nw, model, MHproposal, eta0, control, 
                                        verbose) {
  
  Clist <- ergm.Cprepare(nw, model)

  z <- NULL

  if(control$parallel==0){
    flush.console()
    z <- ergm.mcmcslave(Clist,MHproposal,eta0,control,verbose)
    
    # Note that status code 1 (MCMC_TOO_MANY_EDGES) is handled in ergm.mcmcslave.
    
    if(z$status == 2){ # MCMC_MH_FAILED
      # MH proposal failed somewhere. Throw an error.
      stop("Sampling failed due to a Metropolis-Hastings proposal failing.")
    }
    
    if(!is.null(z$burnin.failed) && z$burnin.failed) warning("Burn-in failed to converge after retries.")
    
    statsmatrix <- matrix(z$s, nrow=control$MCMC.samplesize,
                          ncol=Clist$nstats,
                          byrow = TRUE)
    newnetwork <- newnw.extract(nw,z)
    newnetworks <- list(newnetwork)
  }else{
    control.parallel <- control
    control.parallel$samplesize <- round(control$MCMC.samplesize / control$parallel)
    
    cl <- ergm.getCluster(control, verbose)
    #
    #   Run the jobs with rpvm or Rmpi
    #
    flush.console()
    outlist <- clusterCall(cl,ergm.mcmcslave,
                           Clist,MHproposal,eta0,control.parallel,verbose)
    #
    #   Process the results
    #
    statsmatrix <- NULL
    newnetworks <- list()
    for(i in (1:control$parallel)){
      z <- outlist[[i]]
      # Note that status code 1 (MCMC_TOO_MANY_EDGES) is handled in ergm.mcmcslave.

      if(z$status == 2){ # MCMC_MH_FAILED
        # MH proposal failed somewhere. Throw an error.
        stop("Sampling failed due to a Metropolis-Hastings proposal failing.")
      }
      
      if(!is.null(z$burnin.failed) && z$burnin.failed) warning("Burn-in failed to converge after retries.")
      
      statsmatrix <- rbind(statsmatrix,
                           matrix(z$s, nrow=control.parallel$samplesize,
                                  ncol=Clist$nstats,
                                  byrow = TRUE))
    }

    newnetworks[[i]]<-newnetwork<-newnw.extract(nw,z)
    if(verbose){cat("parallel samplesize=",nrow(statsmatrix),"by",
                    control.parallel$samplesize,"\n")}
    
    ergm.stopCluster(cl)
  }

  colnames(statsmatrix) <- model$coef.names

  statsmatrix[is.na(statsmatrix)] <- 0
  list(statsmatrix=statsmatrix, newnetwork=newnetwork, newnetworks=newnetworks)
}





###############################################################################
# The <ergm.mcmcslave> function is that which the slaves will call to perform
# a validation on the mcmc equal to their slave number. It also returns an
# MCMC sample.
###############################################################################

ergm.mcmcslave <- function(Clist,MHproposal,eta0,control,verbose) {
  # A subroutine to allow caller to override some settings or resume
  # from a pervious run.
  dorun <- function(prev.run=NULL, burnin=NULL, samplesize=NULL, interval=NULL, maxedges=NULL){
    
    numnetworks <- 0

    if(is.null(prev.run)){ # Start from Clist
      nedges <- c(Clist$nedges,0,0)
      tails <- Clist$tails
      heads <- Clist$heads
      stats <- rep(0, Clist$nstats)
    }else{ # Pick up where we left off
      nedges <- prev.run$newnwtails[1]
      tails <- prev.run$newnwtails[2:(nedges+1)]
      heads <- prev.run$newnwheads[2:(nedges+1)]
      nedges <- c(nedges,0,0)
      stats <- matrix(prev.run$s,
                      ncol=Clist$nstats,
                      byrow = TRUE)
      stats <- stats[nrow(stats),]
    }

    if(is.null(burnin)) burnin <- control$MCMC.burnin
    if(is.null(samplesize)) samplesize <- control$MCMC.samplesize
    if(is.null(interval)) interval <- control$MCMC.interval
    if(is.null(maxedges)) maxedges <- control$MCMC.init.maxedges

    repeat{
        z <- .C("MCMC_wrapper",
                as.integer(numnetworks), as.integer(nedges),
                as.integer(tails), as.integer(heads),
                as.integer(Clist$n),
                as.integer(Clist$dir), as.integer(Clist$bipartite),
                as.integer(Clist$nterms),
                as.character(Clist$fnamestring),
                as.character(Clist$snamestring),
                as.character(MHproposal$name), as.character(MHproposal$package),
                as.double(c(Clist$inputs,MHproposal$inputs)), as.double(eta0),
                as.integer(samplesize),
                s = as.double(rep(stats, samplesize)),
                as.integer(burnin), 
                as.integer(interval),
                newnwtails = integer(maxedges),
                newnwheads = integer(maxedges),
                as.integer(verbose), as.integer(MHproposal$arguments$constraints$bd$attribs),
                as.integer(MHproposal$arguments$constraints$bd$maxout), as.integer(MHproposal$arguments$constraints$bd$maxin),
                as.integer(MHproposal$arguments$constraints$bd$minout), as.integer(MHproposal$arguments$constraints$bd$minin),
                as.integer(MHproposal$arguments$constraints$bd$condAllDegExact), as.integer(length(MHproposal$arguments$constraints$bd$attribs)),
                as.integer(maxedges),
                status = integer(1),
                PACKAGE="ergm")
        
        # save the results
        z<-list(s=z$s, newnwtails=z$newnwtails, newnwheads=z$newnwheads, status=z$status, maxedges=maxedges)

      if(z$status!=1) return(z) # Handle everything except for MCMC_TOO_MANY_EDGES elsewhere.

      # The following is only executed (and the loop continued) if too many edges.
      maxedges <- maxedges * 10

    }
  }
  
  if(!is.null(control$MCMC.burnin.retries) && control$MCMC.burnin.retries>0){
    library(coda)

    out <- NULL

    for(try in seq_len(control$MCMC.burnin.retries+1)){
      samplesize <- min(control$MCMC.samplesize,control$MCMC.burnin*control$MCMC.burnin.check.last)
      burnin <- ceiling(control$MCMC.burnin*(1-control$MCMC.burnin.check.last))
      interval <- ceiling(control$MCMC.burnin*control$MCMC.burnin.check.last/samplesize)
      out<-dorun(prev.run=out,
                 burnin = burnin,
                 samplesize = samplesize,
                 interval = interval,
                 maxedges = out$maxedges # note that the first time through, maxedges=NULL$maxedges, which is NULL.
                 )
      # Stop if something went wrong.
      if(out$status!=0) return(out)

      # Get the vector of the last burnin draws.
      burnin.stats <- matrix(out$s, nrow=samplesize,
                             ncol=Clist$nstats,
                             byrow = TRUE)[,Clist$diagnosable,drop=FALSE]
      colnames(burnin.stats) <- names(Clist$diagnosable)[Clist$diagnosable==TRUE]

      if(control$MCMC.runtime.traceplot) plot(mcmc(burnin.stats,start=burnin+1,burnin+samplesize*interval,thin=interval),ask=FALSE,smooth=TRUE,density=FALSE)
      
      ## The following function is a modified version of geweke.diag
      ## from the coda R package. The original code is Copyright (C)
      ## 2005-2011 Martyn Plummer, Nicky Best, Kate Cowles, Karen
      ## Vines
      ##
      ## It is incorporated into the ergm package under the terms of
      ## the GPL v3 license.
      ##
      ## coda's implementation uses spectrum0, which is not robust
      ## enough.
      my.geweke.diag<-function (x, frac1 = 0.1, frac2 = 0.5){
        x <- as.mcmc(x)
        xstart <- c(start(x), end(x) - frac2 * (end(x) - start(x)))
        xend <- c(start(x) + frac1 * (end(x) - start(x)), end(x))
        y.variance <- y.mean <- vector("list", 2)
        for (i in 1:2) {
          y <- window(x, start = xstart[i], end = xend[i])
          y.mean[[i]] <- apply(as.matrix(y), 2, mean)
          y.variance[[i]] <- spectrum0.ar(y)$spec/niter(y)
        }
        z <- (y.mean[[1]] - y.mean[[2]])/sqrt(y.variance[[1]] + y.variance[[2]])

        # Output 2-sided P-value, rather than Z score.
        pnorm(abs(z),0,1,lower.tail=FALSE)*2
               
      }

      # Bonferroni adjustment
      failed <- my.geweke.diag(burnin.stats) < control$MCMC.burnin.check.alpha/ncol(burnin.stats)
      # If a statistic didn't mix at all, fail it as well.
      failed[is.na(failed)] <- TRUE
      if(any(failed)){
        cat("Burn-in failed to converge or mixed very poorly for statistics", paste.and(names(Clist$diagnosable[Clist$diagnosable])[failed]), ". Rerunning.\n")
        if(try == control$MCMC.burnin.retries+1) burnin.failed <- TRUE
      }
      else{
        burnin.failed <- FALSE
        break
      }
    }

    # Do the actual sampling run. Note that we've already done the burnin.
    out <- dorun(prev.run=out, burnin=0, maxedges=out$maxedges)
    if(control$MCMC.runtime.traceplot) {
      stats <- matrix(out$s, nrow=control$MCMC.samplesize,
                      ncol=Clist$nstats,
                      byrow = TRUE)[,Clist$diagnosable,drop=FALSE]
      colnames(stats) <- names(Clist$diagnosable)[Clist$diagnosable==TRUE]
      
      plot(mcmc(stats,start=1,end=control$MCMC.samplesize*control$MCMC.interval,thin=control$MCMC.interval),ask=FALSE,smooth=TRUE,density=FALSE)
    }
    out$burnin.failed <- burnin.failed
    out
  } else {
    out <- dorun()
    if(control$MCMC.runtime.traceplot) {
      stats <- matrix(out$s, nrow=control$MCMC.samplesize,
                      ncol=Clist$nstats,
                      byrow = TRUE)[,Clist$diagnosable,drop=FALSE]
      colnames(stats) <- names(Clist$diagnosable)[Clist$diagnosable==TRUE]
      
      plot(mcmc(stats,start=control$MCMC.burnin+1,control$MCMC.burnin+control$MCMC.samplesize*control$MCMC.interval,thin=control$MCMC.interval),ask=FALSE,smooth=TRUE,density=FALSE)
    }
    out
  }
}
