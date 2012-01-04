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
#   MCMCparams     : the list of parameters which tune the MCMC sampling
#                    processes; the recognized components of 'MCMCparams'
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
#    theta.original   : the 'theta0' inputted into this function
#    objective.history: the number of SPSA iterations used
#
################################################################################

stergm.SPSA <- function(theta0, nw, model.form, model.diss, Clist, 
                      gamma0,
                      MCMCparams, MHproposal.form, MHproposal.diss, MT=FALSE,
                      verbose=FALSE){
  eta0 <- ergm.eta(theta0, model.form$etamap)

  if(verbose){
    cat("SPSA algorithm with theta_0 equal to\n")
    print(theta0)
  }
  eta0 <- ergm.eta(theta0, model.form$etamap)

  z <- stergm.SPSA.C(nw, Clist$meanstats, model.form, model.diss, MHproposal.form, MHproposal.diss,
                        eta0, gamma0, MCMCparams, MT=MT, verbose=verbose)

  ve<-with(z,list(coef=eta,sample=NULL,sample.obs=NULL,objective.history=objective.history))
    
  structure(c(ve, list(newnetwork=nw, 
                       theta.original=theta0,
                       network=nw)),
            class="stergm")
}






#############################################################################
# The <stergm.SPSA.C> function serves as a wrapper to the
# <MCMCDynSPSA_MT_wrapper.C> or <MCMCDynSPSA_wrapper.C> C functions, which
# 
#
# --PARAMETERS--
#   g              : a network object
#   meanstats      : the mean statistics to be subtracted from the observed
#                    statistics
#   model.form     : a formation model, as returned by <ergm.getmodel>
#   model.diss     : a dissolution model, as returned by <ergm.getmodel>
#   MHproposal.form: a MHproposal object for the formation process, as
#                    returned by <getMHproposal>
#   MHproposal.diss: a MHproposal object for the dissolution process, as
#                    returned by <getMHproposal>
#   eta0           : the initial and canonical eta formation parameters
#   gamma0         : the intial theta dissolution coefficients
#   MCMCparams     : the list of parameters which tune the MCMC sampling
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
#       meanstats      : the mean statistics presumably; these are ignored
#                        except to use them in the returned list
#   MT             : whether to use a multithreaded or single threaded
#                    SPSA implementation (T or F); TRUE uses multiple
#                    threads, FALSE uses one
#   verbose        : whether this and subsequently called R and C code
#                    should be verbose (T or F); default=FALSE
#   
# --RETURNED--
#   a list with the 3 following components:
#      meanstats: the 'meanstats' from the 'MCMCparams'; note that this is NOT
#                 the 'meanstats' inputted directly to this function
#      eta      : the estimated? eta formation?? coefficients
#      objective.history: the number of SPSA iterations used
#
################################################################################

stergm.SPSA.C <- function(g, meanstats, model.form, model.diss, 
                          MHproposal.form, MHproposal.diss, eta0, gamma0,
                          MCMCparams, MT, verbose) {

  Clist.form <- ergm.Cprepare(g, model.form)
  Clist.diss <- ergm.Cprepare(g, model.diss)
  maxchanges <- max(MCMCparams$maxchanges, Clist.form$nedges)/5
  MCMCparams$maxchanges <- MCMCparams$maxchanges/5
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
          as.integer(Clist.form$nedges), as.integer(Clist.form$maxpossibleedges),
          as.integer(Clist.form$n),
          as.integer(Clist.form$dir), as.integer(Clist.form$bipartite),
          # Order code. 8
          as.integer(Clist.diss$order.code),
          # Formation terms and proposals. 9
          as.integer(Clist.form$nterms), as.character(Clist.form$fnamestring), as.character(Clist.form$snamestring),
          as.character(MHproposal.form$name), as.character(MHproposal.form$package),
          as.double(Clist.form$inputs), eta=as.double(eta0),
          # Formation parameter fitting. 16
          as.double(summary(model.form$formula)-meanstats),
          # Initial deviation. 17
          as.double(MCMCparams$SPSA.a),
          as.double(MCMCparams$SPSA.alpha),
          as.double(MCMCparams$SPSA.A),
          as.double(MCMCparams$SPSA.c),
          as.double(MCMCparams$SPSA.gamma),
          as.integer(MCMCparams$SPSA.iterations),              
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
          as.integer(MCMCparams$SPSA.burnin),
          as.integer(MCMCparams$SPSA.interval),
          as.integer(MCMCparams$burnin),
          # Space for output.
          as.integer(maxchanges),
          objective.history=double(MCMCparams$SPSA.iterations),
          # Verbosity.
          as.integer(verbose), 
          PACKAGE="ergm") 

  eta <- z$eta
  names(eta) <- names(eta0)

  list(meanstats=MCMCparams$meanstats,
       eta=eta,
       objective.history=z$objective.history)
}
