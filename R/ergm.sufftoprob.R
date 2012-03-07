#  File ergm/R/ergm.sufftoprob.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
########################################################################
# The <ergm.sufftoprob> function attaches a probability weight to each
# row of a given matrix or "mcmc" class; the resultant matrix may be 
# optionally 'compressed'
########################################################################
"ergm.sufftoprob"<- function(suff, compress=FALSE, probs=FALSE) {
  cnames <- dimnames(suff)[[2]]
  if(is.null(cnames)){
    cnames <- 1:ncol(suff)
  }
  if(!is.matrix(suff)){
    csuff <- table(suff)
    csuff <- cbind(as.numeric(names(csuff)),as.numeric(csuff))
    csuff[,2] <- csuff[,2]/sum(csuff[,2])
  }else{
    if(compress){
      if(probs){
        oprobs <- suff[,ncol(suff)]
        cnames <- cnames[-ncol(suff)]
        suff <- suff[,-ncol(suff)]
      }
      cbase <- apply(suff,2,min)
      suff <- sweep(suff, 2, cbase, "-")
      baseten <- apply(suff,2,max) + 1
      baseten <- c(rev(cumprod(rev(baseten[-1]))),1)
      if(probs){
        out <- tapply(oprobs,as.vector(round(suff %*% baseten,7)),sum)
      }else{
        suff <- as.vector(round(suff %*% baseten,7))
        out <- table(suff)
        rm(suff)
      }
      ress <- as.numeric(names(out)) 
      csuff <- matrix(0,ncol=length(baseten),nrow=length(ress))
      for(i in seq(along=baseten)){
        ri <- trunc(ress/baseten[i]+1e-10)
        csuff[,i] <- ri
        ress <- ress - baseten[i]*ri
      }
      csuff <- cbind(sweep(csuff,2,cbase,"+"),out/sum(out))
    }else{
      if(probs){
        csuff <- suff
        csuff[,ncol(suff)] <- csuff[,ncol(suff)]/sum(csuff[,ncol(suff)])
        cnames <- cnames[-length(cnames)]
      }else{
        csuff <- cbind(suff,1/nrow(suff))
      }
    }
  }
  dimnames(csuff) <- list(NULL, c(cnames,"prob"))
  csuff
}




