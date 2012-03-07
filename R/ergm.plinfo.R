#  File ergm/R/ergm.plinfo.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
###############################################################################
# The <ergm.plinfo> function prepares via <plinfo_wrapper.C> 2 of the
# components needed for MPL estimation
###############################################################################

ergm.plinfo<-function(Clist, Clist.miss, m, fix=NULL, theta.offset=NULL)
{
  numobs <- Clist$ndyads

  z <- .C("plinfo_wrapper",
          as.integer(Clist$tails),    as.integer(Clist$heads),
          as.integer(Clist$nedges), 
          as.integer(Clist$n), 
          as.integer(Clist$dir),     as.integer(Clist$nstats), 
          as.character(Clist$fnamestring),as.character(Clist$snamestring),
	  as.double(Clist$inputs),        
	  y = double(numobs),  x = double(numobs*Clist$nstats),
          start=as.integer(1), end=as.integer(numobs),
          PACKAGE="ergm")

  xmat <- matrix(z$x, numobs, Clist$nstats, byrow=TRUE)
  zy <- z$y

  if(Clist.miss$nedges>0){
    temp <- matrix(0,ncol=Clist$n,nrow=Clist$n)
    base <- cbind(as.vector(col(temp)), as.vector(row(temp)))
    base <- base[base[, 2] > base[, 1], ]
    if(!Clist.miss$dir){
      base <- base[base[, 2] > base[, 1], ]
    }else{
      base <- base[base[, 2] != base[, 1], ]
    }
    ubase <- base[,1] + Clist$n*base[,2]
    dmiss <- !is.na(match(ubase, Clist.miss$heads+Clist.miss$tails*Clist$n))
    xmat <- matrix(xmat[!dmiss,], ncol=Clist$nstats, nrow=sum(!dmiss))
    zy <- zy[!dmiss]
  }
  
  colnames(xmat) <- m$coef.names
  
#
# Adjust for the offset
#
  if(is.null(fix) || all(!fix) ){
   fix <- rep(FALSE, length=ncol(xmat))
   theta.offset <- rep(0, length=ncol(xmat))
   names(theta.offset) <- m$coef.names
   foffset <- rep(0, length=nrow(xmat))
   fxmat <- xmat
  }else{
   foffset <- xmat %*% theta.offset
   fxmat <- xmat[,!fix]
  }

# mplefit <- glm(zy ~ .-1 + offset(foffset), data=data.frame(fxmat),
#       family=binomial)

  list(zy=zy, xmat=xmat)
}
