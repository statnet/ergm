duration.matrix <- function(gsim,
                            actornames=if(is.bipartite(gsim$networks))
                            c("Male","Female") else c("Ego","Alter")) {
  if (!is.network((g0<-gsim$networks)) || class(gsim) != "network.series") {
    stop("This function requires that the argument be a network.series")
  }
  cha <- gsim$changed[,c(1,2,3)]
  N <- max(cha[,1])+1
  Nbip <- g0%n%"bipartite"
  if(is.null(Nbip)) Nbip<-0
  Ntot <- g0%n%"n"
  edges <- as.edgelist(g0)
  # Workaround --- if a network is undirected (or bipartite), force heads<tails.
  if(!is.directed(g0)) edges<-t(apply(edges,1,sort))
  nedge <- nrow(edges)
  nchange <- nrow(cha)

  allties <- .C("DurationMatrix", as.integer(nedge), as.integer(edges), 
                 as.integer(N), as.integer(Ntot), 
                 as.integer(nchange), as.integer(cha),
                 duration = as.integer(rep(0,5*(nedge+nchange))),
                 PACKAGE = "ergm")$duration
  allties <- matrix(allties, ncol=5)
  colnames(allties) <- c(actornames, "Start", "End", "Noncensored")
  allties <- allties[allties[,1]!=0,] # Get rid of unused rows
  allties[allties==-1] <- N # Edges that didn't end will be censored at N
  
  as.data.frame(allties <- cbind(allties, duration=allties[,4]-allties[,3]))
}
