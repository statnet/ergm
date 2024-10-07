#  File R/ergm.phase12.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################
###############################################################################
# The <ergm.phase12> function is a wrapper for the <MCMC.phase12.C> method,
# which collects a sample of networks and returns the matrix of summary
# statistics
#
# --PARAMETERS--
#   s         : ergm state object
#   theta0    : the vector of initial theta coefficients
#   control   : a list of control parameters for the MCMC algorithm;
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
            as.integer(control$MCMC.burnin), as.integer(control$MCMC.interval),
            as.double(control$SA.initial_gain), as.integer(control$SA.phase1_n), as.integer(control$SA.nsubphases),
            as.integer(control$SA.min_iterations), as.integer(control$SA.max_iterations),
            as.integer(deInf(NVL(control$MCMC.maxedges,Inf),"maxint")),
            as.integer(verbose),
            PACKAGE="ergm")
    else
      .Call("WtMCMCPhase12",
            s,
            # Phase12 settings
            as.double(deInf(theta0)),
            as.integer(control$MCMC.burnin), as.integer(control$MCMC.interval),
            as.double(control$SA.initial_gain), as.integer(control$SA.phase1_n), as.integer(control$SA.nsubphases),
            as.integer(control$SA.min_iterations), as.integer(control$SA.max_iterations),
            as.integer(deInf(NVL(control$MCMC.maxedges,Inf),"maxint")),
            as.integer(verbose),
            PACKAGE="ergm")

  if(z$status) return(z) # If there is an error.

  theta <- z$theta
  names(theta) <- names(theta0)
  list(status = z$status, theta=theta, state=update(z$state))
}
