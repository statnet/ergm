#  File ergm/R/stergm.getMCMCsample.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
############################################################################
# The <stergm.getMCMCsample> function collects a sample of networks and
# returns the formation and dissolution statistics of each sample, along with
# a toggle matrix of the changes needed from the original network to each
# in the sample
############################################################################
  
stergm.getMCMCsample <- function(nw, model.form, model.diss, model.mon,
                                  MHproposal.form, MHproposal.diss, eta.form, eta.diss, control, 
                                  verbose){


  #
  #   Check for truncation of the returned edge list
  #
  Clist.form <- ergm.Cprepare(nw, model.form)
  Clist.diss <- ergm.Cprepare(nw, model.diss)
  if(!is.null(model.mon)) Clist.mon <- ergm.Cprepare(nw, model.mon)

  collect.form<-if(!is.null(control$collect.form)) control$collect.form else TRUE
  collect.diss<-if(!is.null(control$collect.diss)) control$collect.diss else TRUE

  maxedges <- control$MCMC.init.maxedges
  maxchanges <- control$MCMC.init.maxchanges
  
  repeat{
    #FIXME: Separate MCMC control parameters and properly attach them.
    z <- .C("MCMCDyn_wrapper",
            # Observed network.
            as.integer(Clist.form$tails), as.integer(Clist.form$heads),
            time = if(is.null(Clist.form$time)) as.integer(0) else as.integer(Clist.form$time),
            lasttoggle = if(is.null(Clist.form$lasttoggle)) as.integer(NULL) else as.integer(Clist.form$lasttoggle), 
            as.integer(Clist.form$nedges),
            as.integer(Clist.form$n),
            as.integer(Clist.form$dir), as.integer(Clist.form$bipartite),
            # Formation terms and proposals.
            as.integer(Clist.form$nterms), 
            as.character(Clist.form$fnamestring),
            as.character(Clist.form$snamestring),
            as.character(MHproposal.form$name), as.character(MHproposal.form$package),
            as.double(Clist.form$inputs), as.double(eta.form),
            # Dissolution terms and proposals.
            as.integer(Clist.diss$nterms), 
            as.character(Clist.diss$fnamestring),
            as.character(Clist.diss$snamestring),
            as.character(MHproposal.diss$name), as.character(MHproposal.diss$package),
            as.double(Clist.diss$inputs), as.double(eta.diss),
            # Monitored terms.
            if(!is.null(model.mon)) as.integer(Clist.mon$nterms) else as.integer(0), 
            if(!is.null(model.mon)) as.character(Clist.mon$fnamestring),
            if(!is.null(model.mon)) as.character(Clist.mon$snamestring),
            if(!is.null(model.mon)) as.double(Clist.mon$inputs),
            # Degree bounds.
            as.integer(MHproposal.form$arguments$constraints$bd$attribs), 
            as.integer(MHproposal.form$arguments$constraints$bd$maxout), as.integer(MHproposal.form$arguments$constraints$bd$maxin),
            as.integer(MHproposal.form$arguments$constraints$bd$minout), as.integer(MHproposal.form$arguments$constraints$bd$minin),
            as.integer(MHproposal.form$arguments$constraints$bd$condAllDegExact), as.integer(length(MHproposal.form$arguments$constraints$bd$attribs)),
            # MCMC settings.
            as.double(control$time.samplesize), as.integer(control$MCMC.burnin),
            as.double(control$time.burnin), as.double(control$time.interval),
            # Space for output.
            s.form = if(collect.form) double(Clist.form$nstats*(control$time.samplesize+1)) else double(0),
            s.diss = if(collect.diss) double(Clist.diss$nstats*(control$time.samplesize+1)) else double(0),
            s.mon = if(!is.null(model.mon)) double(Clist.mon$nstats*(control$time.samplesize+1)) else double(0),
            as.integer(maxedges),
            newnwtails = integer(maxchanges), newnwheads = integer(maxchanges), 
            as.integer(maxchanges),
            diffnwtime = integer(maxchanges),
            diffnwtails = integer(maxchanges),
            diffnwheads = integer(maxchanges),
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
  
  statsmatrix.form <-
    if(collect.form) matrix(z$s.form, nrow=control$time.samplesize+1,
                            ncol=Clist.form$nstats,
                            byrow = TRUE)[-1,,drop=FALSE]
    else
      NULL
  
  statsmatrix.diss <-
    if(collect.diss) matrix(z$s.diss, nrow=control$time.samplesize+1,
                            ncol=Clist.diss$nstats,
                            byrow = TRUE)[-1,,drop=FALSE]
  else
    NULL

  statsmatrix.mon <-
    if(!is.null(model.mon))
      matrix(z$s.mon, nrow=control$time.samplesize+1,
             ncol=Clist.mon$nstats,
             byrow = TRUE)[-1,,drop=FALSE]
    else
      NULL
  

  newnetwork<-newnw.extract(nw,z)
  newnetwork %n% "time" <- z$time
  newnetwork %n% "lasttoggle" <- z$lasttoggle

  diffedgelist<-if(!is.null(control$toggles)) {
    if(z$diffnwtails[1]>0){
      tmp <- cbind(z$diffnwtime[2:(z$diffnwtime[1]+1)],z$diffnwtails[2:(z$diffnwheads[1]+1)],z$diffnwheads[2:(z$diffnwtails[1]+1)])
      colnames(tmp) <- c("time","tail","head")
      tmp
    }else{
      tmp <- matrix(0, ncol=3, nrow=0)
      colnames(tmp) <- c("time","tail","head")
      tmp
    }
  }else{
    NULL
  }

  if(!is.null(statsmatrix.form)) colnames(statsmatrix.form) <- model.form$coef.names
  if(!is.null(statsmatrix.diss)) colnames(statsmatrix.diss) <- model.diss$coef.names
  if(!is.null(model.mon)) colnames(statsmatrix.mon) <- model.mon$coef.names

  list(statsmatrix.form=statsmatrix.form, statsmatrix.diss=statsmatrix.diss, statsmatrix.mon=statsmatrix.mon,
       newnetwork=newnetwork,
       changed=diffedgelist,
       maxchanges=control$MCMC.maxchanges)
}
