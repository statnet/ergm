degree.mix.matrix <- function(gsim) {
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
  degmixmat <- .C("DegreeMixMatrix", as.integer(Ntot), 
                  as.integer(nedge), as.integer(edges), 
                  as.integer(N), as.integer(Nfem), as.integer(Ntot), 
                  as.integer(nchange), as.integer(cha),
                  as.integer(ndissolve), as.integer(dissolve),
                  dmm = as.integer(rep(0,4*Ntot)))$dmm
  degmixmat <- matrix(degmixmat, ncol=4)
  colnames(degmixmat) <- c("Deg0", "bothMono", "partnerConc", "Conc")
  degmixmat
}


