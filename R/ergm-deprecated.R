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
#' @param nw,arglist,...,fname,varnames,vartypes,defaultvalues,required,nw.bipartiteflag,requirement,extramessage,nw.directedflag,m,response,formula,object,new,from.new,x,model,proposal,eta0,control,verbose,update.nws,Clist,prev.run,burnin,samplesize,interval,maxedges,attrname,label,coord,jitter,thresh,usearrows,mode,displayisolates,interactive,xlab,ylab,xlim,ylim,pad,label.pad,displaylabels,boxed.labels,label.pos,label.bg,vertex.sides,vertex.rot,arrowhead.cex,label.cex,loop.cex,vertex.cex,edge.col,label.col,vertex.col,label.border,vertex.border,edge.lty,label.lty,vertex.lty,edge.lwd,label.lwd,edge.len,edge.curve,edge.steps,loop.steps,object.scale,uselen,usecurve,suppress.axes,vertices.last,layout.par,cex.main,cex.sub,seed,latent.control,colornames,latent,constraints,constraints.obs,prototype
#' Arguments to deprecated functions.
#' @keywords misc
NULL

#' @rdname ergm-deprecated
#' @export colMeans.mcmc.list
colMeans.mcmc.list <- function(...){
  .dep_once("statnet.common::colMeans.mcmc.list()")
  statnet.common::colMeans.mcmc.list(...)
}

#' @rdname ergm-deprecated
#' @export sweep.mcmc.list
sweep.mcmc.list <- function(...){
  .dep_once("statnet.common::sweep.mcmc.list()")
  statnet.common::sweep.mcmc.list(...)
}

.dep_method <- local({
  warned <- c()
  function(generic, class){
    fullname <- paste(generic,class,sep=".")
    if(! fullname%in%warned){
      me <- sys.call(-1)[[1]]
      if(length(me)>1 && me[[1]]=="::") me <- me[[3]]
      parent <- sys.call(-2)[[1]]
      if(length(parent)>1 && parent[[1]]=="::") parent <- parent[[3]]
      if(me==fullname && NVL(parent,"")!=generic){
        do.call(".Deprecated", list(msg=paste0("You appear to be calling ", fullname,"() directly. ", fullname,"() is a method, and will not be exported in a future version of ", sQuote("ergm"),". Use ", generic, "() instead, or getS3method() if absolutely necessary."), old=fullname))
        warned <<- c(warned, fullname)
      }
    }
  }
})

# Only evaluate deprecation warning once per function.
.dep_once <- local({
  warned <- c()
  function(...){
    me <- sys.call(-1)
    myname <- as.character(me[[1]])
    if(length(myname)>1 && myname[[1]]=="::") myname <- myname[[3]]
    if(! myname%in%warned){
      do.call(".Deprecated", modifyList(list(old=myname),list(...)))
      warned <<- c(warned, myname)
    }
  }
})

#' @rdname ergm-deprecated
#' @export get.miss.dyads
get.miss.dyads <- function(constraints, constraints.obs){
  .dep_once("as.rlebdm() ergm_conlist method")
  # Returns an network indicating which dyads are missing.
  NVL3(as.rlebdm(constraints, constraints.obs=constraints.obs, which="missing"),
       as.network(as.edgelist(.), matrix.type="edgelist",directed=TRUE),
       network.initialize(2))
}
