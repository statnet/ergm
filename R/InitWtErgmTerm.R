#  File R/InitWtErgmTerm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2018 Statnet Commons
#######################################################################

binary_wrap <- function(InitFun, nw, a, valued_args, ddd, namemap = identity, cnmap = identity){
  namemap <- as_mapper(namemap)
  cnmap <- as_mapper(cnmap)

  for(arg in valued_args) a[[arg]] <- NULL

  term <- do.call(InitFun,c(list(nw, a),ddd))
  term$name <- namemap(term$name)
  term$coef.names <- cnmap(term$coef.names)
  term
}

binary_dind_wrap <- function(name, nw, a, ..., cn=name){
  form<-match.arg(a$form,c("sum","nonzero"))
  binary_wrap(get(paste0("InitErgmTerm.",name), mode="function"), nw, a, "form", list(...), namemap=~paste(.,form,sep="_"), cnmap=~sub(cn,paste(cn,form,sep="."), .))
}

InitWtErgmTerm.absdiff <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=NULL, bipartite=NULL,
                      varnames = c("attrname","pow","form"),
                      vartypes = c("character","numeric","character"),
                      defaultvalues = list(NULL,1,"sum"),
                      required = c(TRUE,FALSE,FALSE))
  binary_dind_wrap("absdiff", nw, a, ...)
}

InitWtErgmTerm.absdiffcat <- function(nw, arglist, response, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=NULL, bipartite=NULL,
                      varnames = c("attrname","base","form"),
                      vartypes = c("character","numeric","character"),
                      defaultvalues = list(NULL,NULL,"sum"),
                      required = c(TRUE,FALSE,FALSE))
  binary_dind_wrap("absdiffcat", nw, a, ..., cn="absdiff")
}


InitWtErgmTerm.atleast<-function(nw, arglist, response, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("threshold"),
                      vartypes = c("numeric"),
                      defaultvalues = list(0),
                      required = c(FALSE))
  list(name="atleast",
       coef.names=paste("atleast",a$threshold,sep="."),
       inputs=a$threshold,
       dependence=FALSE,
       minval=0, maxval=network.dyadcount(nw,FALSE),
       emptynwstats=if(0>=a$threshold) network.dyadcount(nw,FALSE) else 0)
}

InitWtErgmTerm.atmost<-function(nw, arglist, response, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("threshold"),
                      vartypes = c("numeric"),
                      defaultvalues = list(0),
                      required = c(FALSE))
  list(name="atmost",
       coef.names=paste("atmost",a$threshold,sep="."),
       inputs=a$threshold,
       dependence=FALSE,
       minval=0, maxval=network.dyadcount(nw,FALSE),
       emptynwstats=if(0<=a$threshold) network.dyadcount(nw,FALSE) else 0)
}

InitWtErgmTerm.b1cov<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE, 
                      varnames = c("attrname","transform","transformname","form"),
                      vartypes = c("character","function","character","character"),
                      defaultvalues = list(NULL,function(x)x,"","sum"),
                      required = c(TRUE,FALSE,FALSE,FALSE))
  binary_dind_wrap("b1cov", nw, a, ...)
}

InitWtErgmTerm.b1factor<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                      varnames = c("attrname", "base", "levels", "form"),
                      vartypes = c("character", "numeric", "character,numeric,logical", "character"),
                      defaultvalues = list(NULL, 1, NULL, "sum"),
                      required = c(TRUE, FALSE, FALSE, FALSE))                                    
  binary_dind_wrap("b1factor", nw, a, ...)
}

InitWtErgmTerm.b2cov<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE, 
                      varnames = c("attrname","transform","transformname","form"),
                      vartypes = c("character","function","character","character"),
                      defaultvalues = list(NULL,function(x)x,"","sum"),
                      required = c(TRUE,FALSE,FALSE,FALSE))
  binary_dind_wrap("b2cov", nw, a, ...)
}

InitWtErgmTerm.b2factor<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                      varnames = c("attrname", "base", "levels", "form"),
                      vartypes = c("character", "numeric", "character,numeric,logical", "character"),
                      defaultvalues = list(NULL, 1, NULL, "sum"),
                      required = c(TRUE, FALSE, FALSE, FALSE))
  binary_dind_wrap("b2factor", nw, a, ...)
}

InitWtErgmTerm.diff <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=NULL, bipartite=NULL,
                      varnames = c("attrname","pow", "dir", "sign.action", "form"),
                      vartypes = c("character","numeric", "character", "character", "character"),
                      defaultvalues = list(NULL,1, "t-h", "identity", "sum"),
                      required = c(TRUE, FALSE, FALSE, FALSE, FALSE))
  binary_dind_wrap("diff", nw, a, ...)
}

InitWtErgmTerm.edgecov <- function(nw, arglist, response, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, 
                      varnames = c("x", "attrname", "form"),
                      vartypes = c("matrix,network,character", "character", "character"),
                      defaultvalues = list(NULL, NULL, "sum"),
                      required = c(TRUE, FALSE, FALSE))
  binary_dind_wrap("edgecov", nw, a, ...)
}


InitWtErgmTerm.equalto<-function(nw, arglist, response, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("value", "tolerance"),
                      vartypes = c("numeric", "numeric"),
                      defaultvalues = list(0, 0),
                      required = c(FALSE,FALSE))
  list(name="ininterval",
       coef.names=paste("equalto",a$value,"pm",a$tolerance,sep="."),
       inputs=with(a, c(value-tolerance, value+tolerance, FALSE, FALSE)),
       dependence=FALSE,
       minval=0, maxval=network.dyadcount(nw,FALSE),
       emptynwstats=if(abs(a$value)<=a$tolerance) network.dyadcount(nw,FALSE) else 0)
}


InitWtErgmTerm.ininterval<-function(nw, arglist, response, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("lower","upper","open"),
                      vartypes = c("numeric","numeric","logical,character"),
                      defaultvalues = list(-Inf,+Inf,c(TRUE,TRUE)),
                      required = c(FALSE,FALSE,FALSE))

  open <- switch(mode(a$open),
                 character = paste0(a$open, collapse=""),
                 logical = rep(a$open, length.out=2))

  OPENSPECS = c('()', '(]', '[)', '[]')
  if(is(open, "character")){
    if(! open%in%OPENSPECS) ergm_Init_abort("Interval openness specification via a string must be ", paste.and(OPENSPECS,'"','"',"or"),".")
    open <- c(substr(open,1,1)=="(",
              substr(open,2,2)==")")
  }

  list(name="ininterval",
       coef.names=paste("ininterval",if(open[1]) "(" else "[", a$lower,",",a$upper, if(open[2]) ")" else "]",sep=""),
       inputs=c(.deinf(a$lower),.deinf(a$upper),open),
       dependence=FALSE,
       minval=0, maxval=network.dyadcount(nw,FALSE),
       emptynwstats=if(
       ((open[1] & 0>a$lower) | (!open[1] & 0>=a$lower)) &
       ((open[2] & 0<a$upper) | (!open[2] & 0<=a$upper))
       ) network.dyadcount(nw,FALSE) else 0)
}

InitWtErgmTerm.greaterthan<-function(nw, arglist, response, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("threshold"),
                      vartypes = c("numeric"),
                      defaultvalues = list(0),
                      required = c(FALSE))
  list(name="greaterthan",
       coef.names=paste("greaterthan",a$threshold,sep="."),
       inputs=a$threshold,
       dependence=FALSE,
       minval=0, maxval=network.dyadcount(nw,FALSE),
       emptynwstats=if(0>a$threshold) network.dyadcount(nw,FALSE) else 0)
}

InitWtErgmTerm.smallerthan<-function(nw, arglist, response, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("threshold"),
                      vartypes = c("numeric"),
                      defaultvalues = list(0),
                      required = c(FALSE))
  list(name="smallerthan",
       coef.names=paste("smallerthan",a$threshold,sep="."),
       inputs=a$threshold,
       dependence=FALSE,
       minval=0, maxval=network.dyadcount(nw,FALSE),
       emptynwstats=if(0<a$threshold) network.dyadcount(nw,FALSE) else 0)
}


InitWtErgmTerm.sum<-function(nw, arglist, response, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("pow"),
                      vartypes = c("numeric"),
                      defaultvalues = list(1),
                      required = c(FALSE))
  if(a$pow==1){
    list(name="sum",
         coef.names="sum",
         inputs=NULL,
         dependence=FALSE)
  }else{
    list(name="sum_pow",
         coef.names=paste("sum",a$pow,sep=""),
         inputs=a$pow,
         dependence=FALSE)
  }
}

InitWtErgmTerm.nodecovar<-function (nw, arglist, response, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, bipartite = FALSE,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = NULL,
                      required = NULL)
  ### Process the arguments

  list(name="nodecovar",
       coef.names = "nodecovar",
       dependence = TRUE
       )
}

InitWtErgmTerm.nodesqrtcovar<-function (nw, arglist, response, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, bipartite = FALSE, nonnegative=TRUE, response=response,
                      varnames = c("center"),
                      vartypes = c("logical"),
                      defaultvalues = list(TRUE),
                      required = c(TRUE))
  ### Process the arguments

  name <- "nodesqrtcovar"

  if(a$center) name <- paste(name,"centered",sep="_")

  coef.name <- gsub("_",".",name)
  
  list(name=name,
       coef.names = coef.name,
       dependence = TRUE,
       minval = if(a$center) NULL else 0
       )
}

InitWtErgmTerm.nodeosqrtcovar<-function (nw, arglist, response, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed = TRUE, nonnegative=TRUE, response=response,
                      varnames = c(),
                      vartypes = c(),
                      defaultvalues = list(),
                      required = c())
  ### Process the arguments

  list(name="nodeosqrtcovar",
       coef.names = "nodeosqrtcovar",
       dependence = TRUE,
       # arithmetic mean >= geometric mean
       minval = 0
       )
}

InitWtErgmTerm.nodeisqrtcovar<-function (nw, arglist, response, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed = TRUE, nonnegative=TRUE, response=response,
                      varnames = c(),
                      vartypes = c(),
                      defaultvalues = list(FALSE),
                      required = c(FALSE))
  ### Process the arguments

  list(name="nodeisqrtcovar",
       coef.names = "nodeisqrtcovar",
       dependence = TRUE,
       # arithmetic mean >= geometric mean
       minval = 0
       )
}

InitWtErgmTerm.nodefactor<-function (nw, arglist, response, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, 
                      varnames = c("attrname", "base","form","levels"),
                      vartypes = c("character", "numeric","character", "character,numeric,logical"),
                      defaultvalues = list(NULL, 1, "sum", NULL),
                      required = c(TRUE, FALSE, FALSE, FALSE))
  binary_dind_wrap("nodefactor", nw, a, ...)
}


InitWtErgmTerm.sociality<-function (nw, arglist, response, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed = FALSE,
                      varnames = c("attrname", "base", "levels", "form"),
                      vartypes = c("character", "numeric", "character,numeric,logical", "character"),
                      defaultvalues = list(NULL, 1, NULL, "sum"),
                      required = c(FALSE, FALSE, FALSE, FALSE))
  binary_dind_wrap("sociality", nw, a, ...)
}


InitWtErgmTerm.nodeocovar<-function (nw, arglist, response, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed = TRUE,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = NULL,
                      required = NULL)
  ### Process the arguments

  list(name="nodeocovar",
       coef.names = "nodeocovar",
       dependence = TRUE
       )
}

InitWtErgmTerm.nodeofactor<-function (nw, arglist, response, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed = TRUE,
                      varnames = c("attrname", "base", "levels", "form"),
                      vartypes = c("character", "numeric", "character,numeric,logical", "character"),
                      defaultvalues = list(NULL, 1, NULL, "sum"),
                      required = c(TRUE, FALSE, FALSE, FALSE))
  binary_dind_wrap("nodeofactor", nw, a, ...)
}

InitWtErgmTerm.sender<-function (nw, arglist, response, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed = TRUE,
                      varnames = c("base", "form"),
                      vartypes = c("numeric", "character"),
                      defaultvalues = list(1, "sum"),
                      required = c(FALSE, FALSE))
  binary_dind_wrap("sender", nw, a, ...)
}

InitWtErgmTerm.nodeicovar<-function (nw, arglist, response, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed = TRUE,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = NULL,
                      required = NULL)
  ### Process the arguments

  list(name="nodeicovar",
       coef.names = "nodeicovar",
       dependence = TRUE
       )
}

InitWtErgmTerm.nodeifactor<-function (nw, arglist, response, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed = TRUE,
                      varnames = c("attrname", "base","levels","form"),
                      vartypes = c("character", "numeric","character,numeric,logical","character"),
                      defaultvalues = list(NULL, 1, NULL, "sum"),
                      required = c(TRUE, FALSE, FALSE, FALSE))
  binary_dind_wrap("nodeifactor", nw, a, ...)
}

InitWtErgmTerm.receiver<-function (nw, arglist, response, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed = TRUE,
                      varnames = c("base", "form"),
                      vartypes = c("numeric", "character"),
                      defaultvalues = list(1, "sum"),
                      required = c(FALSE, FALSE))
  binary_dind_wrap("receiver", nw, a, ...)
}


InitWtErgmTerm.nodematch<-InitWtErgmTerm.match<-function (nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, 
                      varnames = c("attrname", "diff", "keep", "levels", "form"),
                      vartypes = c("character", "logical", "numeric", "character,numeric,logical", "character"),
                      defaultvalues = list(NULL, FALSE, NULL, NULL, "sum"),
                      required = c(TRUE, FALSE, FALSE, FALSE, FALSE))
  binary_dind_wrap("nodematch", nw, a, ...)
}

InitWtErgmTerm.nodemix<-function (nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attrname", "base", "b1levels", "b2levels", "form"),
                      vartypes = c("character", "numeric", "character,numeric,logical", "character,numeric,logical", "character"),
                      defaultvalues = list(NULL, NULL, NULL, NULL, "sum"),
                      required = c(TRUE, FALSE, FALSE, FALSE, FALSE))
  binary_dind_wrap("nodemix", nw, a, ...)
}

InitWtErgmTerm.nodecov<-InitWtErgmTerm.nodemain<-function (nw, arglist, response, ...) {
  a <- check.ErgmTerm(nw, arglist,
                     varnames = c("attrname","transform","transformname","form"),
                     vartypes = c("character","function","character","character"),
                     defaultvalues = list(NULL,identity,"","sum"),
                     required = c(TRUE,FALSE,FALSE,FALSE))
  binary_dind_wrap("nodecov", nw, a, ...)
}

InitWtErgmTerm.nodeicov<-function (nw, arglist, response, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                     varnames = c("attrname","transform","transformname","form"),
                     vartypes = c("character","function","character","character"),
                     defaultvalues = list(NULL,identity,"","sum"),
                     required = c(TRUE,FALSE,FALSE,FALSE))
  binary_dind_wrap("nodeicov", nw, a, ...)
}


InitWtErgmTerm.nodeocov<-function (nw, arglist, response, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                     varnames = c("attrname","transform","transformname","form"),
                     vartypes = c("character","function","character","character"),
                     defaultvalues = list(NULL,identity,"","sum"),
                     required = c(TRUE,FALSE,FALSE,FALSE))
  binary_dind_wrap("nodeocov", nw, a, ...)
}

InitWtErgmTerm.edges<-InitWtErgmTerm.nonzero<-function(nw, arglist, response, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  list(name="nonzero",
       coef.names="nonzero",
       inputs=NULL,
       dependence=FALSE,
       minval=0, maxval=network.dyadcount(nw,FALSE))
}

InitWtErgmTerm.mutual<-function (nw, arglist, response, ...) {
### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=TRUE, bipartite=NULL,
                      varnames = c("form","threshold"),
                      vartypes = c("character","numeric"),
                      defaultvalues = list("min",0),
                      required = c(FALSE,FALSE))

  form <- match.arg(a$form,c("min","nabsdiff","threshold","product","geometric"))
  
  list(name=switch(form,min="mutual_wt_min",nabsdiff="mutual_wt_nabsdiff",threshold="mutual_wt_threshold",product="mutual_wt_product", geometric="mutual_wt_geom_mean"),
       coef.names=switch(form,min="mutual.min",nabsdiff="mutual.nabsdiff",threshold=paste("mutual",a$threshold,sep="."), product="mutual.product",geometric="mutual.geom.mean"),
       inputs=if(form=="threshold") a$threshold,
       dependence=TRUE,
       minval=switch(form,min=NULL,nabsdiff=NULL,threshold=0,product=NULL,geometric=0),
       maxval=switch(form,min=NULL,nabsdiff=0,threshold=NULL,product=NULL,geometric=NULL),
       emptynwstats=switch(form,min=0,nabsdiff=0,threshold=if(a$threshold<=0) network.dyadcount(nw, FALSE)/2 else 0,product=0,geometric=0)
       )
}

InitWtErgmTerm.transitiveties<-function (nw, arglist, response, ...) {
### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, bipartite=NULL,
                      varnames = c("threshold"),
                      vartypes = c("numeric"),
                      defaultvalues = list(0),
                      required = c(FALSE))
  
  list(name="transitiveweights_threshold",
       coef.names="transitiveties",
       inputs=c(a$threshold),
       dependence=TRUE,
       minval=0, maxval=network.dyadcount(nw,FALSE))  
}

InitWtErgmTerm.transitiveweights<-function (nw, arglist, response, ...) {
### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, bipartite=NULL, nonnegative=TRUE,
                      varnames = c("twopath","combine","affect"),
                      vartypes = c("character","character","character"),
                      defaultvalues = list("min","max","min"),
                      required = c(FALSE,FALSE,FALSE), response=response)
  twopaths<-c("min","geomean")
  twopath<-match.arg(a$twopath,twopaths)
  combines<-c("max","sum")
  combine<-match.arg(a$combine,combines)
  affects<-c("min","geomean")
  affect<-match.arg(a$affect,affects)

  list(name="transitiveweights",
       coef.names=paste("transitiveweights",twopath,combine,affect,sep="."),
       inputs=c(
         which(twopaths==twopath),
         which(combines==combine),
         which(affects==affect)
         ),
       dependence=TRUE,
       minval = 0)
}

InitWtErgmTerm.cyclicalties<-function (nw, arglist, response, ...) {
### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, bipartite=NULL,
                      varnames = c("threshold"),
                      vartypes = c("numeric"),
                      defaultvalues = list(0),
                      required = c(FALSE))
  list(name="cyclicweights_threshold",
       coef.names="cyclicalties",
       inputs=c(a$threshold),
       dependence=TRUE,
       minval=0, maxval=network.dyadcount(nw,FALSE))
  
}

InitWtErgmTerm.cyclicalweights<-function (nw, arglist, response, ...) {
### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, bipartite=NULL, nonnegative=TRUE,
                      varnames = c("twopath","combine","affect"),
                      vartypes = c("character","character","character"),
                      defaultvalues = list("min","max","min"),
                      required = c(FALSE,FALSE,FALSE), response=response)
  twopaths<-c("min","geomean")
  twopath<-match.arg(a$twopath,twopaths)
  combines<-c("max","sum")
  combine<-match.arg(a$combine,combines)
  affects<-c("min","geomean")
  affect<-match.arg(a$affect,affects)

  list(name="cyclicalweights",
       coef.names=paste("cyclicalweights",twopath,combine,affect,sep="."),
       inputs=c(
         which(twopaths==twopath),
         which(combines==combine),
         which(affects==affect)
         ),
       dependence=TRUE,
       minval = 0)
}

InitWtErgmTerm.mm<-function (nw, arglist, response, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attrs", "levels", "levels2", "form"),
                      vartypes = c(ERGM_VATTR_SPEC, ERGM_LEVELS_SPEC, ERGM_LEVELS_SPEC, "character"),
                      defaultvalues = list(NULL, NULL, NULL, "sum"),
                      required = c(TRUE, FALSE, FALSE, FALSE))
  binary_dind_wrap("mm", nw, a, ...)
}
