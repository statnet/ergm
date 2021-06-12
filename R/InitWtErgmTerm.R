#  File R/InitWtErgmTerm.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################

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
  a <- a[!attr(a,"missing")] # Remove missing arguments before propagating.
  form<-match.arg(a$form,c("sum","nonzero"))
  binary_wrap(get(paste0("InitErgmTerm.",name), mode="function"), nw, a, "form", list(...), namemap=~paste(.,form,sep="_"), cnmap=~sub(cn,paste(cn,form,sep="."), .))
}

InitWtErgmTerm.absdiff <- function(nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, directed=NULL, bipartite=NULL,
                        varnames = c("attrname","pow","form"),
                        vartypes = c("character","numeric","character"),
                        defaultvalues = list(NULL,1, "sum"),
                        required = c(TRUE,FALSE,FALSE))
  }else{
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, directed=NULL, bipartite=NULL,
                        varnames = c("attr","pow","form"),
                        vartypes = c(ERGM_VATTR_SPEC,"numeric","character"),
                        defaultvalues = list(NULL,1, "sum"),
                        required = c(TRUE,FALSE,FALSE))
  }

  binary_dind_wrap("absdiff", nw, a, ..., version=version)
}

InitWtErgmTerm.absdiffcat <- function(nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, directed=NULL, bipartite=NULL,
                        varnames = c("attrname","base","form"),
                        vartypes = c("character","numeric","character"),
                        defaultvalues = list(NULL,NULL, "sum"),
                        required = c(TRUE,FALSE,FALSE))
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=NULL, bipartite=NULL,
                        varnames = c("attr","base","levels","form"),
                        vartypes = c(ERGM_VATTR_SPEC,"numeric",ERGM_LEVELS_SPEC,"character"),
                        defaultvalues = list(NULL,NULL,NULL, "sum"),
                        required = c(TRUE,FALSE,FALSE,FALSE))
  }

  binary_dind_wrap("absdiffcat", nw, a, ..., cn="absdiff", version=version)
}


InitWtErgmTerm.atleast<-function(nw, arglist, ...) {
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
       emptynwstats=ifelse(0>=a$threshold, network.dyadcount(nw,FALSE), 0))
}

InitWtErgmTerm.atmost<-function(nw, arglist, ...) {
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
       emptynwstats=ifelse(0<=a$threshold, network.dyadcount(nw,FALSE), 0))
}

InitWtErgmTerm.b1cov<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE, 
                        varnames = c("attrname","transform","transformname","form"),
                        vartypes = c("character","function","character","character"),
                        defaultvalues = list(NULL,function(x)x,"", "sum"),
                        required = c(TRUE,FALSE,FALSE,FALSE))
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE, 
                        varnames = c("attr","form"),
                        vartypes = c(ERGM_VATTR_SPEC,"character"),
                        defaultvalues = list(NULL, "sum"),
                        required = c(TRUE,FALSE))
  }    
  binary_dind_wrap("b1cov", nw, a, ..., version=version)
}

InitWtErgmTerm.b1factor<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                        varnames = c("attrname", "base", "levels","form"),
                        vartypes = c("character", "numeric", "character,numeric,logical","character"),
                        defaultvalues = list(NULL, 1, NULL, "sum"),
                        required = c(TRUE, FALSE, FALSE,FALSE))
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                        varnames = c("attr", "base","levels","form"),
                        vartypes = c(ERGM_VATTR_SPEC, "numeric",ERGM_LEVELS_SPEC,"character"),
                        defaultvalues = list(NULL, 1, LEVELS_BASE1, "sum"),
                        required = c(TRUE, FALSE,FALSE,FALSE))
  }
                              
  binary_dind_wrap("b1factor", nw, a, ..., version=version)
}

InitWtErgmTerm.b1sociality<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                      varnames = c("nodes", "form"),
                      vartypes = c(ERGM_LEVELS_SPEC, "character"),
                      defaultvalues = list(LEVELS_BASE1, "sum"),
                      required = c(FALSE, FALSE))
  binary_dind_wrap("b1sociality", nw, a, ...)
                      
}

InitWtErgmTerm.b2cov<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                        varnames = c("attrname","transform","transformname","form"),
                        vartypes = c("character","function","character","character"),
                        defaultvalues = list(NULL,function(x)x,"", "sum"),
                        required = c(TRUE,FALSE,FALSE,FALSE))
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                        varnames = c("attr","form"),
                        vartypes = c(ERGM_VATTR_SPEC,"character"),
                        defaultvalues = list(NULL, "sum"),
                        required = c(TRUE,FALSE))
  }

  binary_dind_wrap("b2cov", nw, a, ..., version=version)
}

InitWtErgmTerm.b2factor<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                        varnames = c("attrname", "base", "levels","form"),
                        vartypes = c("character", "numeric", "character,numeric,logical","character"),
                        defaultvalues = list(NULL, 1, NULL, "sum"),
                        required = c(TRUE, FALSE, FALSE,FALSE))
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                        varnames = c("attr", "base", "levels","form"),
                        vartypes = c(ERGM_VATTR_SPEC, "numeric", ERGM_LEVELS_SPEC,"character"),
                        defaultvalues = list(NULL, 1, LEVELS_BASE1, "sum"),
                        required = c(TRUE, FALSE,FALSE,FALSE))
  }

  binary_dind_wrap("b2factor", nw, a, ..., version=version)
}


InitWtErgmTerm.b2sociality<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                      varnames = c("nodes", "form"),
                      vartypes = c(ERGM_LEVELS_SPEC, "character"),
                      defaultvalues = list(LEVELS_BASE1, "sum"),
                      required = c(FALSE, FALSE))
  binary_dind_wrap("b2sociality", nw, a, ...)
}



InitWtErgmTerm.diff <- function(nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, directed=NULL, bipartite=NULL,
                        varnames = c("attrname","pow", "dir", "sign.action","form"),
                        vartypes = c("character","numeric", "character", "character","character"),
                        defaultvalues = list(NULL,1, "t-h", "identity", "sum"),
                        required = c(TRUE, FALSE, FALSE, FALSE,FALSE))
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=NULL, bipartite=NULL,
                        varnames = c("attr","pow", "dir", "sign.action","form"),
                        vartypes = c(ERGM_VATTR_SPEC,"numeric", "character", "character","character"),
                        defaultvalues = list(NULL,1, "t-h", "identity", "sum"),
                        required = c(TRUE, FALSE, FALSE, FALSE,FALSE))
  }  

  binary_dind_wrap("diff", nw, a, ..., version=version)
}

InitWtErgmTerm.edgecov <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, 
                      varnames = c("x", "attrname", "form"),
                      vartypes = c("matrix,network,character", "character", "character"),
                      defaultvalues = list(NULL, NULL, "sum"),
                      required = c(TRUE, FALSE, FALSE))
  binary_dind_wrap("edgecov", nw, a, ...)
}


InitWtErgmTerm.equalto<-function(nw, arglist, ...) {
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


InitWtErgmTerm.ininterval<-function(nw, arglist, ...) {
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
       inputs=c(deInf(a$lower),deInf(a$upper),open),
       dependence=FALSE,
       minval=0, maxval=network.dyadcount(nw,FALSE),
       emptynwstats=if(
       ((open[1] & 0>a$lower) | (!open[1] & 0>=a$lower)) &
       ((open[2] & 0<a$upper) | (!open[2] & 0<=a$upper))
       ) network.dyadcount(nw,FALSE) else 0)
}

InitWtErgmTerm.greaterthan<-function(nw, arglist, ...) {
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
       emptynwstats=ifelse(0>a$threshold, network.dyadcount(nw,FALSE), 0))
}

InitWtErgmTerm.smallerthan<-function(nw, arglist, ...) {
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
       emptynwstats=ifelse(0<a$threshold, network.dyadcount(nw,FALSE), 0))
}


InitWtErgmTerm.sum<-function(nw, arglist, ...) {
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

InitWtErgmTerm.nodecovar<-function (nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed = FALSE,
                      varnames = c("center","transform"),
                      vartypes = c("logical","character"),
                      defaultvalues = list(FALSE,"identity"),
                      required = c(FALSE,FALSE))
  ### Process the arguments
  transcode <- match(match.arg(a$transform, c("identity","sqrt")), c("identity","sqrt"))-1  
  list(name="nodecovar",
       coef.names = "nodecovar",
       inputs = c(transcode,a$center), # Transformation codes: 0 for no transformation, 1 for square root.
       dependence = TRUE
       )
}

InitWtErgmTerm.nodeosqrtcovar<-function (nw, arglist, ...) {
  .Deprecated('nodeocovar(transform="sqrt")',
            old = "nodeosqrtcovar")
  arglist$transform <- "sqrt"
  InitWtErgmTerm.nodeocovar(nw, arglist, ...)
}

InitWtErgmTerm.nodeisqrtcovar<-function (nw, arglist, ...) {
  .Deprecated('nodeicovar(transform="sqrt")',
            old = "nodeisqrtcovar")
  arglist$transform <- "sqrt"
  InitWtErgmTerm.nodeicovar(nw, arglist, ...)
}

InitWtErgmTerm.nodesqrtcovar<-function (nw, arglist, ...) {
  .Deprecated('nodecovar(transform="sqrt")',
            old = "nodesqrtcovar")
  arglist$transform <- "sqrt"
  InitWtErgmTerm.nodecovar(nw, arglist, ...)
}

InitWtErgmTerm.nodefactor<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  ### Check the network and arguments to make sure they are appropriate.
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist,
                        varnames = c("attrname", "base", "levels","form"),
                        vartypes = c("character", "numeric", "character,numeric,logical","character"),
                        defaultvalues = list(NULL, 1, NULL, "sum"),
                        required = c(TRUE, FALSE, FALSE,FALSE))
  }else{
    a <- check.ErgmTerm(nw, arglist,
                        varnames = c("attr","base", "levels","form"),
                        vartypes = c(ERGM_VATTR_SPEC,"numeric", ERGM_LEVELS_SPEC,"character"),
                        defaultvalues = list(NULL,1, LEVELS_BASE1, "sum"),
                        required = c(TRUE, FALSE,FALSE,FALSE))
  }

  binary_dind_wrap("nodefactor", nw, a, ..., version=version)
}


InitWtErgmTerm.sociality<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  ### Check the network and arguments to make sure they are appropriate.
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=FALSE,
                        varnames = c("attrname", "base", "levels","form"),
                        vartypes = c("character", "numeric", "character,numeric,logical","character"),
                        defaultvalues = list(NULL, 1, NULL, "sum"),
                        required = c(FALSE, FALSE, FALSE,FALSE))
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=FALSE,
                        varnames = c("attr", "base", "levels", "nodes","form"),
                        vartypes = c(ERGM_VATTR_SPEC, "numeric", ERGM_LEVELS_SPEC, ERGM_LEVELS_SPEC,"character"),
                        defaultvalues = list(NULL, 1, NULL, LEVELS_BASE1, "sum"),
                        required = c(FALSE, FALSE, FALSE, FALSE, FALSE))  
  }

  binary_dind_wrap("sociality", nw, a, ..., version=version)
}


InitWtErgmTerm.nodeocovar<-function (nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed = TRUE,
                      varnames = c("center","transform"),
                      vartypes = c("logical","character"),
                      defaultvalues = list(FALSE,"identity"),
                      required = c(FALSE,FALSE))
  ### Process the arguments
  transcode <- match(match.arg(a$transform, c("identity","sqrt")), c("identity","sqrt"))-1  
  list(name="nodeocovar",
       coef.names = "nodeocovar",
       inputs = c(transcode,a$center), # Transformation codes: 0 for no transformation, 1 for square root.
       dependence = TRUE
       )
}

InitWtErgmTerm.nodeofactor<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  ### Check the network and arguments to make sure they are appropriate.
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=TRUE, 
                        varnames = c("attrname", "base", "levels","form"),
                        vartypes = c("character", "numeric", "character,numeric,logical","character"),
                        defaultvalues = list(NULL, 1, NULL, "sum"),
                        required = c(TRUE, FALSE, FALSE,FALSE))
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=TRUE, 
                        varnames = c("attr", "base", "levels","form"),
                        vartypes = c(ERGM_VATTR_SPEC, "numeric", ERGM_LEVELS_SPEC,"character"),
                        defaultvalues = list(NULL, 1, LEVELS_BASE1, "sum"),
                        required = c(TRUE, FALSE, FALSE,FALSE))
  }

  binary_dind_wrap("nodeofactor", nw, a, ..., version=version)
}

InitWtErgmTerm.sender<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  ### Check the network and arguments to make sure they are appropriate.
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                        varnames = c("base","form"),
                        vartypes = c("numeric","character"),
                        defaultvalues = list(1, "sum"),
                        required = c(FALSE,FALSE))
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                        varnames = c("base","nodes","form"),
                        vartypes = c("numeric",ERGM_LEVELS_SPEC,"character"),
                        defaultvalues = list(1,LEVELS_BASE1, "sum"),
                        required = c(FALSE,FALSE,FALSE))
  }  

  binary_dind_wrap("sender", nw, a, ..., version=version)
}

InitWtErgmTerm.nodeicovar<-function (nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed = TRUE,
                      varnames = c("center","transform"),
                      vartypes = c("logical","character"),
                      defaultvalues = list(FALSE,"identity"),
                      required = c(FALSE,FALSE))
  ### Process the arguments
  transcode <- match(match.arg(a$transform, c("identity","sqrt")), c("identity","sqrt"))-1  
  list(name="nodeicovar",
       coef.names = "nodeicovar",
       inputs = c(transcode,a$center), # Transformation codes: 0 for no transformation, 1 for square root.
       dependence = TRUE
       )
}

InitWtErgmTerm.nodeifactor<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  ### Check the network and arguments to make sure they are appropriate.
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=TRUE, 
                        varnames = c("attrname", "base", "levels","form"),
                        vartypes = c("character", "numeric", "character,numeric,logical","character"),
                        defaultvalues = list(NULL, 1, NULL, "sum"),
                        required = c(TRUE, FALSE, FALSE,FALSE))
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=TRUE, 
                        varnames = c("attr", "base","levels","form"),
                        vartypes = c(ERGM_VATTR_SPEC, "numeric",ERGM_LEVELS_SPEC,"character"),
                        defaultvalues = list(NULL, 1, LEVELS_BASE1, "sum"),
                        required = c(TRUE, FALSE,FALSE,FALSE))
  }  

  binary_dind_wrap("nodeifactor", nw, a, ..., version=version)
}

InitWtErgmTerm.receiver<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  ### Check the network and arguments to make sure they are appropriate.
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                        varnames = c("base","form"),
                        vartypes = c("numeric","character"),
                        defaultvalues = list(1, "sum"),
                        required = c(FALSE,FALSE))
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                        varnames = c("base","nodes","form"),
                        vartypes = c("numeric",ERGM_LEVELS_SPEC,"character"),
                        defaultvalues = list(1,LEVELS_BASE1, "sum"),
                        required = c(FALSE,FALSE,FALSE))
  }

  binary_dind_wrap("receiver", nw, a, ..., version=version)
}


InitWtErgmTerm.nodematch<-InitWtErgmTerm.match<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  ### Check the network and arguments to make sure they are appropriate.
  if(version <= as.package_version("3.9.4")){
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, 
                        varnames = c("attrname", "diff", "keep", "levels","form"),
                        vartypes = c("character", "logical", "numeric", "character,numeric,logical","character"),
                        defaultvalues = list(NULL, FALSE, NULL, NULL, "sum"),
                        required = c(TRUE, FALSE, FALSE, FALSE,FALSE))
  }else{
    a <- check.ErgmTerm(nw, arglist, 
                        varnames = c("attr", "diff", "keep", "levels","form"),
                        vartypes = c(ERGM_VATTR_SPEC, "logical", "numeric", ERGM_LEVELS_SPEC,"character"),
                        defaultvalues = list(NULL, FALSE, NULL, NULL, "sum"),
                        required = c(TRUE, FALSE, FALSE, FALSE,FALSE))
  }

  binary_dind_wrap("nodematch", nw, a, ..., version=version)
}

InitWtErgmTerm.nodemix<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist,
                        varnames = c("attrname", "base", "b1levels", "b2levels", "form"),
                        vartypes = c("character", "numeric", "character,numeric,logical", "character,numeric,logical", "character"),
                        defaultvalues = list(NULL, NULL, NULL, NULL, "sum"),
                        required = c(TRUE, FALSE, FALSE, FALSE, FALSE))
  }else{
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist,
                        varnames = c("attr", "base", "b1levels", "b2levels", "levels", "levels2", "form"),
                        vartypes = c(ERGM_VATTR_SPEC, "numeric", ERGM_LEVELS_SPEC, ERGM_LEVELS_SPEC, ERGM_LEVELS_SPEC, ERGM_LEVELS_SPEC, "character"),
                        defaultvalues = list(NULL, NULL, NULL, NULL, NULL, NULL, "sum"),
                        required = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))  
  }
  binary_dind_wrap("nodemix", nw, a, ..., version=version)
}

InitWtErgmTerm.nodecov<-InitWtErgmTerm.nodemain<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist,
                        varnames = c("attrname","transform","transformname","form"),
                        vartypes = c("character","function","character","character"),
                        defaultvalues = list(NULL,function(x)x,"", "sum"),
                        required = c(TRUE,FALSE,FALSE,FALSE))
  }else{
    a <- check.ErgmTerm(nw, arglist,
                        varnames = c("attr","form"),
                        vartypes = c(ERGM_VATTR_SPEC,"character"),
                        defaultvalues = list(NULL, "sum"),
                        required = c(TRUE,FALSE))
  }
  binary_dind_wrap("nodecov", nw, a, ..., version=version)
}

InitWtErgmTerm.nodeicov<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                        varnames = c("attrname","transform","transformname","form"),
                        vartypes = c("character","function","character","character"),
                        defaultvalues = list(NULL,identity,"", "sum"),
                        required = c(TRUE,FALSE,FALSE,FALSE))
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                        varnames = c("attr","form"),
                        vartypes = c(ERGM_VATTR_SPEC,"character"),
                        defaultvalues = list(NULL, "sum"),
                        required = c(TRUE,FALSE))
  }

  binary_dind_wrap("nodeicov", nw, a, ..., version=version)
}


InitWtErgmTerm.nodeocov<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=TRUE, 
                        varnames = c("attrname","transform","transformname","form"),
                        vartypes = c("character","function","character","character"),
                        defaultvalues = list(NULL,identity,"", "sum"),
                        required = c(TRUE,FALSE,FALSE,FALSE))
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=TRUE, 
                        varnames = c("attr","form"),
                        vartypes = c(ERGM_VATTR_SPEC,"character"),
                        defaultvalues = list(NULL, "sum"),
                        required = c(TRUE,FALSE))
  }

  binary_dind_wrap("nodeocov", nw, a, ..., version=version)
}


InitWtErgmTerm.edges<-InitWtErgmTerm.nonzero<-function(nw, arglist, ...) {
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

InitWtErgmTerm.mutual<-function (nw, arglist, ...) {
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

InitWtErgmTerm.transitiveties<-function (nw, arglist, ...) {
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

InitWtErgmTerm.transitiveweights<-function (nw, arglist, ...) {
### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, bipartite=NULL, nonnegative=TRUE,
                      varnames = c("twopath","combine","affect"),
                      vartypes = c("character","character","character"),
                      defaultvalues = list("min","max","min"),
                      required = c(FALSE,FALSE,FALSE))
  twopaths<-c("min","geomean")
  twopath<-match.arg(a$twopath,twopaths)
  combines<-c("max", "sum")
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

InitWtErgmTerm.cyclicalties<-function (nw, arglist, ...) {
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

InitWtErgmTerm.cyclicalweights<-function (nw, arglist, ...) {
### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, bipartite=NULL, nonnegative=TRUE,
                      varnames = c("twopath","combine","affect"),
                      vartypes = c("character","character","character"),
                      defaultvalues = list("min","max","min"),
                      required = c(FALSE,FALSE,FALSE))
  twopaths<-c("min","geomean")
  twopath<-match.arg(a$twopath,twopaths)
  combines<-c("max", "sum")
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

InitWtErgmTerm.mm<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attrs", "levels", "levels2", "form"),
                      vartypes = c(ERGM_VATTR_SPEC, ERGM_LEVELS_SPEC, ERGM_LEVELS_SPEC, "character"),
                      defaultvalues = list(NULL, NULL, NULL, "sum"),
                      required = c(TRUE, FALSE, FALSE, FALSE))
  binary_dind_wrap("mm", nw, a, ...)
}


################################################################################


