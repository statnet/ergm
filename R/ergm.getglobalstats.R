ergm.getglobalstats <- function(nw, m) {
  Clist <- ergm.Cprepare(nw, m)
  #
  #    Calculate the global statistics
  #
  gs <- .C("network_stats_wrapper",
           as.integer(Clist$heads), as.integer(Clist$tails), 
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
  
  #
  # Adjust to global values
  #
                                                                
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

  gs
}


