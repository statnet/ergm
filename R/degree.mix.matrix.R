##########################################################################
# The <degree.mix.matrix> function calculates a degree mixing matrix,
# via <DegreeMixMatrix.C>,  that shows the patterns of concurrency
# throughout a network list
#
# --PARAMETERS--
#   gsim: a bipartite network list, as returned by <simulate.stergm>
#
# --RETURNED--
#   degmixmat: an nx4 matrix, where the rows correspond to the n nodes
#              in each network of 'gsim' and the columns are:
#      "Deg0"       : the number of networks (time steps) in which the 
#                     node had no partnerships
#      "BothMono"   : the number of time steps in which the node was 
#                     mononogamous and the partner was monongamous
#      "partnerConc": the number of time steps in which the node was
#                     monogamous, but the partner had concurrent
#                     partnerships
#      "Conc"       : the number of time steps in which the node had 
#                     concurrent partnerships
#
###########################################################################

degree.mix.matrix <- function(gsim) {
  if (!is.network((g0<-gsim$networks)) || class(gsim) != "network.list") {
    stop("This function requires that the argument be a network.list")
  }
  
  if(!is.bipartite(gsim$networks))
    stop("This function only works for bipartite networks at the moment.")

  cha <- gsim$changed[,c(1,3,2)]
  cha <- cha[order(gsim$change[,1]),]
  N <- max(cha[,1])+1
  Nfem <- g0%n%"bipartite"
  Ntot <- g0%n%"n"
  edges <- as.edgelist(g0)[,2:1] 
  nedge <- nrow(edges)
  nchange <- nrow(cha)
  degmixmat <- .C("DegreeMixMatrix", as.integer(Ntot), 
                  as.integer(nedge), as.integer(edges), 
                  as.integer(N), as.integer(Nfem), as.integer(Ntot), 
                  as.integer(nchange), as.integer(cha),
                  dmm = as.integer(rep(0,4*Ntot)),
                  PACKAGE = "ergm")$dmm
  degmixmat <- matrix(degmixmat, ncol=4)
  colnames(degmixmat) <- c("Deg0", "bothMono", "partnerConc", "Conc")
  degmixmat
}


