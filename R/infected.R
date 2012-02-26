###########################################################################
# The <prevalence> function calculates and returns the prevalence of ?? via
# <Prevalence.C> 
#
# --PARAMETERS--
#   gsim       : a bipartite network list
#   nsim       : the number of simulations to use; this is currently ignored
#                but should surely be passed to <Prevalence.C>; default=1
#   beta       : the rate of transmission; default=0.1
#   randomseeds: whether random seeds should be used (T or F); default=FALSE
#
# --RETURNED--
#   the prevalence as a list containing:
#      prevalence:
#      infected  :
#
#  without the R and C calls matching up, this is hard to determine
###########################################################################

prevalence <- function(gsim, nsim=1, beta=0.1, randomseeds=FALSE) {
  if (!is.network((g0<-gsim$networks)) || class(gsim) != "network.list") {
    stop("This function requires that the argument be a network.list")
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
                  PACKAGE = "ergm")
  list(prevalence=prevalence$prev, infected=prevalence$infected)
}


