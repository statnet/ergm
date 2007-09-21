ergm.godfather <- function(formula, timestamps=NULL, toggles=NULL, sim=NULL,
                           accumulate=FALSE,
                           verbose=FALSE,
                           algorithm.control=list()) {
  # Make the network a proposal it can't refuse.
  # Here, formula is a typical ergm formula (i.e., nw ~ terms)
  # timestamps is a vector whose length is the same as the #rows of toggles.
  # Each toggle has a timestamp, and this function forces the network to make
  # all of the changes at each unique timestamp value (in increasing order)
  # and then it keeps track of the change statistics that result.
  # Thus, the final product is a matrix of change statistics in which the
  # number of rows is determined by the # of unique timestamps and the
  # number of columns is determined by the ERGM terms as usual.

  if(is.null(sim)){
    if(is.null(timestamps) | is.null(toggles)){
      stop("Both 'timestamps' and 'toggle' are required arguments if 'sim' ",
           "is not given.")
    }
  }else{
    timestamps <- sim$changed[,1]
    toggles <- sim$changed[,2:3]
  }
  ## Defaults :
  con <- list(maxedges=100000,
              maxchanges=1000000, 
              assume_consecutive_timestamps=TRUE,
              return_new_network=FALSE
              )
  con[(namc <- names(algorithm.control))] <- algorithm.control

  nw <- ergm.getnetwork(formula)
  m <- ergm.getmodel(formula, nw, drop=FALSE, initialfit=TRUE)
  Clist <- ergm.Cprepare(nw, m)
  Clist$obs <- summary(m$formula)
  ots <- order(timestamps)
  toggles <- matrix(toggles[ots,],ncol=2)
  timestamps <- timestamps[ots]
  mincol = apply(toggles,1,which.min)
  toggles[mincol==2,] <- toggles[mincol==2,2:1] # make sure col1 < col2
  maxedges <- con$return_new_network * max(con$maxedges, 5*Clist$nedges)
  obsstat <- summary(formula)  
  z <- .C("godfather_wrapper",
          as.integer(Clist$heads), as.integer(Clist$tails), 
          as.integer(Clist$nedges), as.integer(Clist$n),
          as.integer(Clist$dir), as.integer(Clist$bipartite),
          as.integer(Clist$nterms), 
          as.character(Clist$fnamestring),
          as.character(Clist$snamestring), 
          as.integer(length(timestamps)), as.integer(timestamps),
          as.integer(toggles[,1]), as.integer(toggles[,2]),
          as.double(Clist$inputs),
          s = double((1+length(unique(timestamps))) * Clist$nparam),
          newnwheads = integer(maxedges+1),
          newnwtails = integer(maxedges+1),
          as.integer(accumulate),
          as.integer(verbose),
          as.integer(maxedges), 
          PACKAGE="ergm")  
  stats <- matrix(z$s + obsstat, ncol=Clist$nparam, byrow=T)
  colnames(stats) <- m$coef.names
  uts <- unique(timestamps)
  if (con$assume_consecutive_timestamps) {
    if (!all(uts %% 1 == 0)) {
      print(paste("ergm.godfather cannnot assume consecutive timestamps",
                  "unless timestamps are integers"))
    } else {
      # Insert timesteps whenever an integer is skipped
      uts2 <- min(uts):max(uts)
      stats <- stats[c(1,1+cumsum(uts2 %in% uts)),] # Don't forget:  first row
      # of original 'stats' matrix is initial value, not matching any timestamp
      uts <- uts2
    }
  }
  out <- list(stats = stats, timestamps = c(NA, uts))
  if (con$return_new_network) { 
    cat("Creating new network...\n")
    out$newnetwork <- newnw.extract(nw,z)
  }
  out
}
