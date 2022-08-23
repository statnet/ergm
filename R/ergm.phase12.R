#  File R/ergm.phase12.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
###############################################################################
# The <ergm.phase12> function is a wrapper for the <MCMC.phase12.C> method,
# which collects a sample of networks and returns the matrix of summary
# statistics
#
# --PARAMETERS--
#   g         : a network object
#   model     : a model for 'g', as returned by <ergm_model>
#   proposal: an proposal object, as returned by <proposal>
#   eta0      : the vector of initial eta coefficients
#   control: a list of control parameters for the MCMC algorithm;
#               recognized components include:
#                     'maxedges'     'samplesize'     'gain'
#                     'stats'        'phase1'         'nsub'
#                     'burnin'       'interval'       'target.stats'
#               the purpose of most of these variables is given in the
#               <control.ergm> function header; 'stats' seems to be
#                used as the mean statistics; 'target.stats' is merely
#                returned.
#   verbose   : whether the C functions should be verbose (T or F)
#
# --RETURNED--
#   a list containing
#     statsmatrix: the matrix of summary statistics
#     newnetwork : the final network sampled
#     target.stats  : the 'target.stats' from 'control'
#     maxedges   : the 'maxedges' from 'control'
#     eta        : the parameters used to produce the sample given
#                  by 'statsmatrix'
#
###############################################################################

ergm.phase12 <- function(s, theta0,
                        control, verbose) {
  on.exit(ergm_Cstate_clear())

  z <-
    if(!is.valued(s))
      .Call("MCMCPhase12",
            s,
            # Phase12 settings
            as.double(deInf(theta0)),
            as.integer(control$MCMC.samplesize), as.integer(control$MCMC.burnin), as.integer(control$MCMC.interval),
            as.double(control$gain), as.integer(control$phase1), as.integer(control$nsub),
            as.integer(deInf(NVL(control$MCMC.maxedges,Inf),"maxint")),
            as.integer(verbose),
            PACKAGE="ergm")
    else
      .Call("WtMCMCPhase12",
            s,
            # Phase12 settings
            as.double(deInf(theta0)),
            as.integer(control$MCMC.samplesize), as.integer(control$MCMC.burnin), as.integer(control$MCMC.interval),
            as.double(control$gain), as.integer(control$phase1), as.integer(control$nsub),
            as.integer(deInf(NVL(control$MCMC.maxedges,Inf),"maxint")),
            as.integer(verbose),
            PACKAGE="ergm")

  statsmatrix <- matrix(z$s, nrow=control$MCMC.samplesize,
                        ncol=nparam(s,canonical=TRUE),
                        byrow = TRUE)
  theta <- z$theta
  names(theta) <- names(theta0)

  z$state <- update(z$state)
  newnetwork<-as.network(z$state)
  
  colnames(statsmatrix) <- param_names(s,canonical=TRUE)
  list(statsmatrix=statsmatrix, newnetwork=newnetwork, target.stats=as.ergm_model(s)$target.stats, nw.stats=as.ergm_model(s)$nw.stats,
       theta=theta, state=z$state)
}
