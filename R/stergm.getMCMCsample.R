############################################################################
# The <stergm.getMCMCsample> function collects a sample of networks and
# returns the formation and dissolution statistics of each sample, along with
# a toggle matrix of the changes needed from the original network to each
# in the sample
#
# --PARAMETERS--
#   nw             : a network object
#   model.form     : a formation model, as returned by <ergm.getmodel>
#   model.diss     : a dissolution model, as returned by <ergm.getmodel>
#   MHproposal.form: a list of parameters needed for MHproposals of the
#                    formations
#   MHproposal.diss: a list of parameters needed for MHproposals of the
#                    dissolutions
#   init.form    : the initial model coefficients for 'model.form'
#   init.diss    : the initial model coefficients for 'model.diss'
#   control     : the list of parameters controlling the MCMC algorithm;
#                    recognized components include:
#      maxchanges    :  this ends up defining the "MCMCDyn workspace", but
#                       in a peculiar way: the 'maxchanges' you give, let's
#                       call this number 'mc' is fed into this equation:
#                           mc < 5^x * mc - 55
#                       the 'x' that solves the equation is then used and
#                       and 'maxchanges' and the MCMCDyn workspace is defined as
#                           5^(x-1) * (the original mc)
#                       if the computed Clist.form$nedges is > than the original
#                       'maxchanges', 'maxchanges' is ignored and this odd little
#                       process occurs to the 'nedeges'
#      parallel      :  the number of threads in which to run sampling
#      samplesize    :  the desired MCMC sample size
#      toggles       :  whether to return the toggle matrix; non-NULL values
#                       will include 'changed' in the return list, NULL
#                       will not
#      MH.burnin     :  this is accepted as 'MH_interval' which determines
#                       the number of proposals in each MCMC step
#      time.burnin   :  the number of MCMC steps to disregard for the burnin
#                       period 
#      time.interval :  the number of MCMC steps to disregard in between
#                       sampled networks
#   verbose        : whether this and other functions should be verbose
#
# Note:  In reality, there should be many fewer arguments to this function,
# since most info should be passed via Clist (this is, after all, what Clist
# is for:  Holding all arguments required for the .C call).  In particular,
# the elements of MHproposal.form, control, verbose, should certainly
# be part of Clist.  But this is a project for another day!
#
# --RETURNED--
#   the MCMC sample as a list containing:
#     statsmatrix.form: the matrix of sampled statistics for 'model.form' RELATIVE TO INITIAL NETWORK
#     statsmatrix.diss: the matrix of sampled statistics for 'model.form' RELATIVE TO INITIAL NETWORK
#     statsmatrix.mon: the matrix of sampled statistics for 'model.mon' RELATIVE TO INITIAL NETWORK
#     newnetwork      : the final network from the sampling process
#     changed         : a toggle matrix, where the first column is
#                       the timestamp of the toggle and the 2nd and 3rd
#                       columns are the head & tail of the toggle; this
#                       is only returned if 'control'$toggles is not NULL
#     maxchanges      : the "MCMC Dyn workspace"; see 'maxchanges' in the
#                       input param list
#
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
            lasttoggle = if(is.null(Clist.form$time)) as.integer(NULL) else as.integer(Clist.form$lasttoggle), 
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
      cbind(z$diffnwtime[2:(z$diffnwtime[1]+1)],z$diffnwtails[2:(z$diffnwheads[1]+1)],z$diffnwheads[2:(z$diffnwtails[1]+1)])
    }else{
      matrix(0, ncol=3, nrow=0)
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
