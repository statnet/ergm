#============================================================================
# This files contains the 2 following functions for fitting stergms using
# simultaneous perturbation stochastic approximation (SPSA)
#      <stergm.SPSA>
#      <stergm.SPSA.C>
#============================================================================




#############################################################################
# The <stergm.SPSA> function fits a stergm using the SPSA method of
# moments estimation.  This is similar to the Robbins-Monro style, but is
# less precise and makes no assumption about the sign of the derivatives of
# each element of the equilibrium expected values of the statistics of
# interest wrt the corresponding formation phase parameter.
#
# --PARAMETERS--
#   theta.form0    : the initial theta formation coefficients
#   nw             : a network object
#   model.form     : a formation model, as returned by <ergm.getmodel>
#   model.diss     : a dissolution model, as returned by <ergm.getmodel>
#   Clist          : the list of inputs that are used by the C code and
#                    returned by <ergm.Cprepare>
#   gamma0         : the intial theta dissolution coefficients
#   control     : the list of parameters which tune the MCMC sampling
#                    processes; the recognized components of 'control'
#                    are those passed to and used by the <stergm.phase12.C>
#                    function below and are described in its function header
#   MHproposal.form: a MHproposal object for the formation process, as
#                    returned by <getMHproposal>
#   MHproposal.diss: a MHproposal object for the dissolution process, as
#                    returned by <getMHproposal>
#   MT             : whether to use a multithreaded or single threaded
#                    SPSA implementation (T or F); TRUE uses multiple
#                    threads, FALSE uses one; default=FALSE
#   verbose        : whether this and subsequently called R and C code
#                    should be verbose (T or F); default=FALSE
#
# --RETURNED--
#   a stergm object as a list containing:
#    eta        : the estimated? eta formation?? coefficients
#    sample     : NULL
#    sample.obs: NULL
#    newnetwork : the 'nw' inputted into this function
#    network    : the 'nw' inputted into this function
#    theta.original   : the 'init' inputted into this function
#    objective.history: the number of SPSA iterations used
#
################################################################################

stergm.SPSA <- function(init, nw, model.form, model.diss,
                      gamma0,
                      control, MHproposal.form, MHproposal.diss, MT=FALSE,
                      verbose=FALSE){
  eta0 <- ergm.eta(init, model.form$etamap)

  if(verbose){
    cat("SPSA algorithm with theta_0 equal to\n")
    print(init)
  }
  eta0 <- ergm.eta(init, model.form$etamap)

  z <- stergm.SPSA.C(nw, model.form$target.stats, model.form, model.diss, MHproposal.form, MHproposal.diss,
                        eta0, gamma0, control, MT=MT, verbose=verbose)

  ve<-with(z,list(coef=eta,sample=NULL,sample.obs=NULL,objective.history=objective.history, newnetwork=nw, 
                       init=init.form,
                       network=nw))
}






#############################################################################
# The <stergm.SPSA.C> function serves as a wrapper to the
# <MCMCDynSPSA_MT_wrapper.C> or <MCMCDynSPSA_wrapper.C> C functions, which
# 
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
#   eta0           : the initial and canonical eta formation parameters
#   gamma0         : the intial theta dissolution coefficients
#   control     : the list of parameters which tune the MCMC sampling
#                    processes; recognized components include:
#       SPSA.iterations: the number of iterations to use in the SPSA sampling
#       SPSA.a         : see the next 2 params
#       SPSA.alpha     : see the next 2 params
#       SPSA.A         : this and the 2 params above help to define the
#                        'gain' paramater as
#                            SPSA.a/(SPSA.A +i +1)^(SPSA.alpha)
#                        where i is indexed from 0 to SPSA.iterations
#       SPSA.c         : see the next param
#       SPSA.gamma     : this and the param above help to define the 'diff'
#                        parameter as
#                             SPSA.c/(i+1)^(SPSA.gamma)
#                        where i is indexed from 0 to SPSA.iterations
#       SPSA.burnin    : the number of MCMC steps to disregard for the burnin
#                        period
#       SPSA.interval  : this is eventually received as 'S' and looks like a
#                        a sample size, rather than an interval, since 'S' controls
#                        the number of MCMC steps that contribute to the stats vector
#       burnin         : this is received as 'dyninterval' by the C code and
#                        eventually used as 'MH_interval' to control the number
#                        of proposals between sampled networks
#       target.stats      : the mean statistics presumably; these are ignored
#                        except to use them in the returned list
#   MT             : whether to use a multithreaded or single threaded
#                    SPSA implementation (T or F); TRUE uses multiple
#                    threads, FALSE uses one
#   verbose        : whether this and subsequently called R and C code
#                    should be verbose (T or F); default=FALSE
#   
# --RETURNED--
#   a list with the 3 following components:
#      target.stats: the 'target.stats' from the 'control'; note that this is NOT
#                 the 'target.stats' inputted directly to this function
#      eta      : the estimated? eta formation?? coefficients
#      objective.history: the number of SPSA iterations used
#
################################################################################

stergm.SPSA.C <- function(g, target.stats, model.form, model.diss, 
                          MHproposal.form, MHproposal.diss, eta0, gamma0,
                          control, MT, verbose) {

  Clist.form <- ergm.Cprepare(g, model.form)
  Clist.diss <- ergm.Cprepare(g, model.diss)
  maxchanges <- max(control$maxchanges, Clist.form$nedges)/5
  control$maxchanges <- control$maxchanges/5
  if(verbose){cat(paste("MCMCDyn workspace is",maxchanges,"\n"))}

  if(MT){
    cat("Using multithreaded SPSA.\n")
    f<-"MCMCDynSPSA_MT_wrapper"
  }
  else{
    cat("Using single-threaded SPSA.\n")
    f<-"MCMCDynSPSA_wrapper"
  }
  
  z <- .C(f,
          # Observed/starting network. 1
          as.integer(Clist.form$tails), as.integer(Clist.form$heads), 
          as.integer(Clist.form$nedges), 
          as.integer(Clist.form$n),
          as.integer(Clist.form$dir), as.integer(Clist.form$bipartite),
          # Formation terms and proposals. 8
          as.integer(Clist.form$nterms), as.character(Clist.form$fnamestring), as.character(Clist.form$snamestring), as.integer(model.form$offset),
          as.character(MHproposal.form$name), as.character(MHproposal.form$package),
          as.double(Clist.form$inputs), eta=as.double(eta0),
          # Formation parameter fitting. 16
          as.double(summary(model.form$formula)-target.stats),
          # Initial deviation. 17
          as.double(control$SPSA.a),
          as.double(control$SPSA.alpha),
          as.double(control$SPSA.A),
          as.double(control$SPSA.c),
          as.double(control$SPSA.gamma),
          as.integer(control$SPSA.iterations),              
          # Dissolution terms and proposals. 21
          as.integer(Clist.diss$nterms), as.character(Clist.diss$fnamestring), as.character(Clist.diss$snamestring),
          as.character(MHproposal.diss$name), as.character(MHproposal.diss$package),
          as.double(Clist.diss$inputs), as.double(gamma0),
          # Degree bounds.
          as.integer(MHproposal.form$arguments$constraints$bd$attribs), 
          as.integer(MHproposal.form$arguments$constraints$bd$maxout), as.integer(MHproposal.form$arguments$constraints$bd$maxin),
          as.integer(MHproposal.form$arguments$constraints$bd$minout), as.integer(MHproposal.form$arguments$constraints$bd$minin),
          as.integer(MHproposal.form$arguments$constraints$bd$condAllDegExact), as.integer(length(MHproposal.form$arguments$constraints$bd$attribs)), 
          # MCMC settings.              
          as.integer(control$SPSA.burnin),
          as.integer(control$SPSA.interval),
          as.integer(control$MCMC.burnin),
          # Space for output.
          as.integer(maxchanges),
          objective.history=double(control$SPSA.iterations),
          # Verbosity.
          as.integer(verbose), 
          PACKAGE="ergm") 

  eta <- z$eta
  names(eta) <- names(eta0)

  list(target.stats=model.form$target.stats,
       coef.form=eta,
       objective.history=z$objective.history)
}
