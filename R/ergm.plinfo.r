ergm.plinfo<-function(Clist, mClist, m, fix=NULL, theta.offset=NULL)
{
  numobs <- Clist$ndyads

  z <- .C("plinfo_wrapper",
          as.double(Clist$heads),    as.double(Clist$tails),
          as.double(Clist$nedges),   as.double(Clist$n), 
          as.integer(Clist$dir),     as.integer(Clist$nparam), 
          as.character(Clist$fnamestring),as.character(Clist$snamestring),
	  as.double(Clist$inputs),        
	  y = double(numobs),  x = double(numobs*Clist$nparam),
          start=as.integer(1), end=as.integer(numobs),
          PACKAGE="statnet")

  xmat <- matrix(z$x, numobs, Clist$nparam, byrow=TRUE)
  zy <- z$y

  if(mClist$nedges>0){
    temp <- matrix(0,ncol=Clist$n,nrow=Clist$n)
    base <- cbind(as.vector(col(temp)), as.vector(row(temp)))
    base <- base[base[, 2] > base[, 1], ]
    if(!mClist$dir){
      base <- base[base[, 2] > base[, 1], ]
    }else{
      base <- base[base[, 2] != base[, 1], ]
    }
    ubase <- base[,1] + Clist$n*base[,2]
    dmiss <- !is.na(match(ubase, mClist$tails+mClist$heads*Clist$n))
    xmat <- matrix(xmat[!dmiss,], ncol=Clist$nparam, nrow=sum(!dmiss))
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
