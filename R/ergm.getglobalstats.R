#  File ergm/R/ergm.getglobalstats.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
###########################################################################
# The <ergm.getglobalstats> function calculates and returns the global
# statistics for a given network via <network_stats_wrapper.C>
#############################################################################

ergm.getglobalstats <- function(nw, m) {
  Clist <- ergm.Cprepare(nw, m)

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
  
  # *** don't forget, tails are passes in first now, not heads  
  gs <- .C("network_stats_wrapper",
           as.integer(Clist$tails), as.integer(Clist$heads),
           as.integer(!is.null(Clist$time) && !is.null(Clist$lasttoggle)), as.integer(Clist$time), as.integer(Clist$lasttoggle),
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


