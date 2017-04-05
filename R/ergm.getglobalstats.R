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

ergm.getglobalstats <- function(nw, m) {
  Clist <- ergm.Cprepare(nw, m)
  # *** don't forget, tails are passes in first now, notheads
  gs <- .C("network_stats_wrapper",
           as.integer(Clist$tails), as.integer(Clist$heads), 
           as.integer(Clist$nedges),
           as.integer(Clist$n),
           as.integer(Clist$dir), as.integer(Clist$bipartite), 
           as.integer(Clist$nterms), 
           as.character(Clist$fnamestring), as.character(Clist$snamestring), 
           as.double(Clist$inputs),
           gs = double(Clist$nstats),
           PACKAGE="ergm"
           )$gs
  names(gs) <- m$coef.names

  # Adjust to global values
                                                                
  # New method:  Use $emptynwstats added to m$terms by the InitErgm function
  # Read the comments at the top of InitErgm.R or InitErgmTerm.R for 
  # an explanation of the $emptynwstats mechanism
  i <- 1
  for (j in 1:length(m$terms)) {
    tmp <- m$term[[j]]
    k <- tmp$inputs[2] # Number of statistics for this model term
    if (!is.null(tmp$emptynwstats)) {
      gs[i:(i+k-1)] <- gs[i:(i+k-1)] + tmp$emptynwstats
    }
    i <- i + k
  }

  tase <- grep("duration",names(gs)) # not currently part of CRAN version
  if(length(tase) >0){
    gs[tase] <- -gs[tase]
  }

  gs
}


