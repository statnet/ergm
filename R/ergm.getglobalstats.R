#  File R/ergm.getglobalstats.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
###########################################################################
# The <ergm.getglobalstats> function calculates and returns the global
# statistics for a given network via <network_stats_wrapper.C>
#
# --PARAMETERS--
#   nw:  a network object
#   m :  the model in use with network nw, as returned by <ergm.getmodel>
#
#
# --RETURNED--
#   gs:  a vector of the global statistics
#
#############################################################################

ergm.getglobalstats <- function(nw, m, response=NULL) {
  Clist <- ergm.Cprepare(nw, m, response=response)

  # Adjust to global values. This needs to happen before the C call,
  # so that an s_function, if exists could override.
                                                                
  # New method:  Use $emptynwstats added to m$terms by the InitErgmTerm function
  # Read the comments at the top of InitErgm.R or InitErgmTerm.R for 
  # an explanation of the $emptynwstats mechanism
  gs <- rep(0, Clist$nstats)
  i <- 1
  for (j in 1:length(m$terms)) {
    tmp <- m$terms[[j]]
    k <- tmp$inputs[2] # Number of statistics for this model term
    if (!is.null(tmp$emptynwstats)) {
      gs[i:(i+k-1)] <- gs[i:(i+k-1)] + tmp$emptynwstats
    }
    i <- i + k
  }

  # Note that the empty network statistics are passed to the C
  # code. The reason is that if an s_??? function exists, it can
  # overwrite them, since it can compute the whole thing, while if
  # only the d_??? function exists, it needs to add on to empty
  # network statistics.
  
  # *** don't forget, tails are passes in first now, notheads  
  gs <- if(is.null(response))
         .C("network_stats_wrapper",
            as.integer(Clist$tails), as.integer(Clist$heads), as.integer(!is.null(Clist$time)), as.integer(Clist$time), as.integer(NVL(Clist$lasttoggle,0)),
            as.integer(Clist$nedges),
            as.integer(Clist$n),
            as.integer(Clist$dir), as.integer(Clist$bipartite), 
            as.integer(Clist$nterms), 
            as.character(Clist$fnamestring), as.character(Clist$snamestring), 
            as.double(Clist$inputs),
            gs = as.double(gs),
            PACKAGE="ergm"
            )$gs
         else
         .C("wt_network_stats_wrapper",
            as.integer(Clist$tails), as.integer(Clist$heads), as.double(Clist$weights), as.integer(!is.null(Clist$time)), as.integer(Clist$time), as.integer(NVL(Clist$lasttoggle,0)),
            as.integer(Clist$nedges),
            as.integer(Clist$n),
            as.integer(Clist$dir), as.integer(Clist$bipartite), 
            as.integer(Clist$nterms), 
            as.character(Clist$fnamestring), as.character(Clist$snamestring), 
            as.double(Clist$inputs),
            gs = as.double(gs),
            PACKAGE="ergm"
            )$gs
  names(gs) <- m$coef.names

  gs
}


