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


  if(verbose) cat("Robbins-Monro algorithm with coef_F_0 = (",theta.form0, ") and coef_D = (",theta.diss,")\n" )
  eta.form0 <- ergm.eta(theta.form0, model.form$etamap)
  eta.diss <- ergm.eta(theta.diss, model.diss$etamap)


  z <- stergm.phase12.C(nw, model.form$target.stats, model.form, model.diss, MHproposal.form, MHproposal.diss,
                        eta.form0, eta.diss, control, verbose=verbose)
  
  #ve<-with(z,list(coef=eta,sample=s$statsmatrix.form,sample.obs=NULL))
  ve<-with(z,list(coef.form=coef.form,coef.diss=theta.diss))
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

stergm.phase12.C <- function(g, target.stats, model.form, model.diss, 
                             MHproposal.form, MHproposal.diss, eta.form0, eta.diss,
                             control, verbose) {
  # ms <- model$target.stats
  # if(!is.null(ms)) {
  #   if (is.null(names(ms)) && length(ms) == length(model.form$coef.names))
  #     names(ms) <- model.form$coef.names
  #   obs <- control$orig.obs
  #   obs <- obs[match(names(ms), names(obs))]
  #   ms  <-  ms[match(names(obs), names(ms))]
  #   matchcols <- match(names(ms), names(obs))
  #   if (any(!is.na(matchcols))) {
  #     ms[!is.na(matchcols)] <- ms[!is.na(matchcols)] - obs[matchcols[!is.na(matchcols)]]
  #   }
  # }
  Clist.form <- ergm.Cprepare(g, model.form)
  Clist.diss <- ergm.Cprepare(g, model.diss)
  maxedges <- max(control$MCMC.init.maxedges, Clist.form$nedges)
  if(verbose){cat(paste("MCMCDyn workspace is",maxedges,"\n"))}
  
  z <- .C("MCMCDynPhase12",
          # Observed/starting network. 1
          as.integer(Clist.form$tails), as.integer(Clist.form$heads), 
          as.integer(Clist.form$nedges), as.integer(Clist.form$maxpossibleedges),
          as.integer(Clist.form$n),
          as.integer(Clist.form$dir), as.integer(Clist.form$bipartite),
          # Formation terms and proposals. 8
          as.integer(Clist.form$nterms), as.character(Clist.form$fnamestring), as.character(Clist.form$snamestring), as.integer(model.form$offset),
          as.character(MHproposal.form$name), as.character(MHproposal.form$package),
          as.double(Clist.form$inputs), eta.form=as.double(eta.form0),
          # Formation parameter fitting. 16
          as.double(summary(model.form$formula)-target.stats),
          as.double(control$RM.init_gain),
          as.integer(control$RM.phase1n_base),
          as.integer(control$RM.phase2n_base),
          as.integer(control$RM.phase2sub),              
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
          as.integer(control$EGMoME.MCMC.burnin),
          # Space for output.
          as.integer(maxedges),
          # Verbosity.
          as.integer(verbose), 
          PACKAGE="ergm") 

  eta.form <- z$eta
  names(eta.form) <- names(eta.form0)

  list(target.stats=model.form$target.stats,
       coef.form=eta.form)
}
