# convert a MCMC matrix of sampled
# sufficient statistics be "compressed" into a reduced
# form. In its natural state the matrix has one row for
# each sampled run. Only the
# unique values of the vector of sufficient statistics
# are retained and an additional column is added to
# the matrix with the proportion of the MCMC runs that
# returned that value.  The objective is to keep the
# size of the matrix small for very long runs to save
# memory, speed calculations and make it easier to read.
#
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




