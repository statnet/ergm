prevalence <- function(gsim, nsim=1, beta=0.1, randomseeds=FALSE) {
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
  prev <- rep(0,Ntot)
  infected <- rep(0,Ntot)
  infected[sample(1:Ntot)] <- 1
  prevalence <- .C("Prevalence", as.integer(Ntot), 
                  as.integer(nedge), as.integer(edges), 
                  as.integer(N), as.integer(Nfem), as.integer(sum(infected)), 
                  as.integer(Ntot), 
                  as.integer(nchange), as.integer(cha),
                  as.integer(ndissolve), as.integer(dissolve),
                  as.integer(randomseeds),
                  as.double(beta),
                  infected = as.integer(infected),
                  prev = as.integer(prev),
                  PACKAGE = "statnet")
  list(prevalence=prevalence$prev, infected=prevalence$infected)
}


