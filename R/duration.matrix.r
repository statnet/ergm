duration.matrix <- function(gsim) {
  if (!is.network((g0<-gsim$networks)) || class(gsim) != "network.series") {
    stop("This function requires that the argument be a network.series")
  }
  cha <- gsim$changed[,c(1,3,2)]
  cha <- cha[order(gsim$change[,1]),]
  N <- max(cha[,1])+1
  Nfem <- g0%n%"bipartite"
  Ntot <- g0%n%"n"
  edges <- as.edgelist(g0)[,2:1] 
  nedge <- nrow(edges)
  nchange <- nrow(cha)
  dissolve <- gsim$dissolve[order(gsim$dissolve[,1]), c(1,2,3)]
  ndissolve <- nrow(dissolve)
  # Note bizarre difference in ordering columns between $change and $dissolve.
  # dissolve is not actually used by the DurationMatrix function at all.

  allties <- .C("DurationMatrix", as.integer(nedge), as.integer(edges), 
                 as.integer(N), as.integer(Nfem), as.integer(Ntot), 
                 as.integer(nchange), as.integer(cha),
                 as.integer(ndissolve), as.integer(dissolve), 
                 duration = as.integer(rep(0,5*(nedge+nchange))))$duration
  allties <- matrix(allties, ncol=5)
  colnames(allties) <- c("Fem", "Male", "Start", "End", "Noncensored")
  allties <- allties[allties[,1]!=0,] # Get rid of unused rows
  allties[allties==-1] <- N # Edges that didn't end will be censored at N
  allties <- cbind(allties, duration=allties[,4]-allties[,3])
}

