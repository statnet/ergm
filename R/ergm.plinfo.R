###############################################################################
# The <ergm.plinfo> function prepares via <plinfo_wrapper.C> 2 of the
# components needed for MPL estimation
#
# --PARAMETERS--
#   Clist       : a list of parameters used for fitting and returned
#                 by <ergm.Cprepare>
#   Clist.miss  : the corresponding 'Clist' for the network of missing
#                 edges returned by <ergm.design>
#   m           : the model, as returned by <ergm.getmodel>
#
#
# --IGNORED--
#   these are essentially ignored, since they are both used to calculate
#   items which are not returned:
#       theta.offset: a logical vector specifying which of the model
#                     coefficients are offset, i.e. fixed
#       fix         : appears to be have the same meaning as 'theta.offset'
#
#
# --RETURNED--
#   the pseudo likelihood info as a list containing
#     xmat: the matrix of change statistics??
#     zy  : the corresponding vector of responses, i.e. tie values
#
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

  list(zy=zy, xmat=xmat)
}
