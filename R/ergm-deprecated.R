#  File R/ergm-deprecated.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2018 Statnet Commons
#######################################################################
# Since deprecated functions are sprinkled all over the source code, this filename is prepended with "AAA" to ensure that it's evaluated before any others.

#' @name ergm-deprecated
#' @rdname ergm-deprecated
#' @title Functions that will no longer be supported in future releases of the package
#' @description Functions that have been superceed, were never documented, or will be removed from the package for other reasons
#' @param nw,arglist,...,fname,varnames,vartypes,defaultvalues,required,nw.bipartiteflag,requirement,extramessage,nw.directedflag,m,response,formula,object,new,from.new,x,model,proposal,eta0,control,verbose,update.nws,Clist,prev.run,burnin,samplesize,interval,maxedges,attrname,label,coord,jitter,thresh,usearrows,mode,displayisolates,interactive,xlab,ylab,xlim,ylim,pad,label.pad,displaylabels,boxed.labels,label.pos,label.bg,vertex.sides,vertex.rot,arrowhead.cex,label.cex,loop.cex,vertex.cex,edge.col,label.col,vertex.col,label.border,vertex.border,edge.lty,label.lty,vertex.lty,edge.lwd,label.lwd,edge.len,edge.curve,edge.steps,loop.steps,object.scale,uselen,usecurve,suppress.axes,vertices.last,layout.par,cex.main,cex.sub,seed,latent.control,colornames,latent
#' Arguments to deprecated functions.
#' @keywords misc
NULL

#' @rdname ergm-deprecated
#' @export colMeans.mcmc.list
colMeans.mcmc.list <- function(...){
  .Deprecated("statnet.common::colMeans.mcmc.list()")
  statnet.common::colMeans.mcmc.list(...)
}

#' @rdname ergm-deprecated
#' @export sweep.mcmc.list
sweep.mcmc.list <- function(...){
  .Deprecated("statnet.common::sweep.mcmc.list()")
  statnet.common::sweep.mcmc.list(...)
}

#' @rdname ergm-deprecated
#' @export lapply.mcmc.list
lapply.mcmc.list <- function(...){
  .Deprecated("statnet.common::lapply.mcmc.list()")
  statnet.common::lapply.mcmc.list(...)
}

.dep_method <- function(generic, class){
  fullname <- paste(generic,class,sep=".")
  me <- sys.call(1)
  parent <- sys.call(2)
  if(me[[1]]==fullname && parent[[1]]!=generic)
    .Deprecated(msg=paste0("You appear to be calling ", fullname,"() directly.", fullname,"() is a method, and will not be exported in a future version of ", sQuote("ergm"),". Use ", generic, "() instead, or getS3method() if absolutely necessary."))

}
