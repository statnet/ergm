#####################################################################
# The <ergm.curved.statsmatrix> maps the stats matrix to its reduced
# form based on the etamap object created by <ergm.etamap>
#
# --PARAMETERS--
#   statsmatrix:  the matrix of sampled summary statistics
#   theta      :  the  model parameters producing
#                 'statsmatrix'
#   etamap     :  the theta-> eta mapping, as returned by <ergm.etamap> 
#                   
# --RETURNED--
#   a 2-element list containing
#     sm   :  the reduced form 'statsmatrix'
#     novar:  whether each row of the stats matrix has any variance (T or F)
#
# author: MSH  1/29/06
#
#######################################################################

"ergm.curved.statsmatrix" <- function(statsmatrix,theta,etamap){
  eta <- rep(0,etamap$etalength)
  ec <- etamap$canonical
  eta[ec[ec>0]] <- theta[ec>0]
  nstats <- sum(ec>0) + length(etamap$curved)

  if(length(etamap$curved)>0) {
    sm <- matrix(0, nrow=length(eta), ncol=nstats)
    namessm <- rep("", nstats)
    icurved <- 0
    istat <- 0
    for(i in 1:nstats) {
      istat <- istat + 1
      if(ec[istat] > 0){
#      sm[istat,i] <- eta[istat] 
       sm[istat,i] <- 1
       namessm[i] <- names(theta)[istat]
      }else{
       icurved <- icurved + 1
       cm <- etamap$curved[[icurved]]
       sm[cm$to,i] <- cm$map(theta[cm$from],length(cm$to),cm$cov)
#      (un)scale by linear coefficient
       sm[cm$to,i] <- sm[cm$to,i] / theta[cm$from][1]
       namessm[i] <- names(theta)[istat]
# MSH: Is this last line right? I have commented it out 7/19/06
#      And changed to "from"
#      istat <- istat + length(cm$to)
       istat <- istat + length(cm$from) - 1
      }
    }
    sm <- statsmatrix %*% sm
    colnames(sm) <- namessm
  }else{
    sm <- statsmatrix
  }
  novar <- apply(sm,2,var)<1e-6
  list(sm=sm, novar=novar)
}
