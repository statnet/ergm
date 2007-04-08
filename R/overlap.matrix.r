overlap.matrix <- function(gsim, maxoverlaps=100000) {
  if (!is.network((g0<-gsim$networks)) || class(gsim) != "network.series") {
    stop("This function requires that the first argument ",
         "be a network (the original) and the second ",
         "argument be a network.series (the changes)")
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

  overlap <- .C("OverlapDurations", as.integer(nedge), as.integer(edges), 
                as.integer(N), as.integer(Nfem), as.integer(Ntot), 
                as.integer(nchange), as.integer(cha),
                as.integer(ndissolve), as.integer(dissolve), as.integer(maxoverlaps),
                omatrix = as.integer(rep(0,8*maxoverlaps)))$omatrix
  overlap <- matrix(overlap, ncol=8)
  colnames(overlap) <- c("Fem1", "Male1", "Fem2", "Male2", "start1", "start2", 
                         "endtime", "firsttoend")
  overlap <- overlap[overlap[,1]>0,] # Get rid of unused rows
  overlap[overlap==-1] <- N # Overlaps that didn't end will be censored at N
  duration <- overlap[,7] - apply(overlap[,5:6], 1, max) 
  fts <- apply(overlap[,c("start1","start2")], 1, 
                        function(x) (x[1]<=x[2]) + 2*(x[2]<=x[1]))
  fte <- overlap[,"firsttoend"]
  type <- (((fts == 1 & fte == 2) | (fts == 2 & fte == 1)) 
           + 2*(fts == 3 | fte == 3 | fte == fts) 
           + 3*(fte==0))
  type <- c("transitional","embedded","truncated")[type]
  data.frame(overlap, duration = duration, firsttostart=fts, type=type)
}



