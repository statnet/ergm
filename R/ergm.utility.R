#  File ergm/R/ergm.utility.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
ostar2deg <- function(object, ninflast=TRUE){
 nnodes <- network.size(object$newnetwork)
 nodeg <- paste("odeg",1:(nnodes-1),sep="")
 nostar <- paste("ostar",1:(nnodes-1),sep="")
 if(ninflast){
  ostar <- rep(-10000,nnodes-1)
 }else{
  ostar <- rep(0,nnodes-1)
 }
 names(ostar) <- nostar
 mostar <- match(nostar,names(object$coef))
 for(i in 1:(nnodes-1)){
  if(!is.na(mostar[i])){
   ostar[i] <- object$coef[mostar[i]]
  }
 }
 odeg <- rep(0,nnodes-1)
 names(odeg) <- nodeg
 for(j in 1:(nnodes-1)){
  odeg[j] <- sum(choose(rep(j,j), 1:j)*ostar[nostar[1:j]])
 }
 odeg
}



is.invertible <- function(V, tol=1e-12)
{
    ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
    all(ev/max(ev) > tol)
}


is.ergm <- function(object)
{
    class(object)=="ergm"
}

###############################################################################
# The <degreedist> function computes and returns the degree distribution for
# a given network
###############################################################################

degreedist <- function(g, print=TRUE)
{
 if(!is.network(g)){
  stop("degreedist() requires a network object")
 }
 if(is.directed(g)){                                      
   mesp <- paste("c(",paste(0:(network.size(g)-1),collapse=","),")",sep="")
   outdegrees <- summary(as.formula(paste('g ~ odegree(',mesp,')',sep="")))
   indegrees <- summary(as.formula(paste('g ~ idegree(',mesp,')',sep="")))
   temp <- outdegrees > 0 | indegrees > 0
   outdegrees <- outdegrees[temp]
   indegrees <- indegrees[temp]
   if(!is.null(outdegrees) & print){print(outdegrees[outdegrees>0])}
   if(!is.null(indegrees) & print){print(indegrees[indegrees>0])}
   degrees <- rbind(indegrees, outdegrees)
 }else{
  if(is.bipartite(g)){
   nb1 <- get.network.attribute(g,"bipartite")
   nb2 <- network.size(g) - nb1
   mesp <- paste("c(",paste(0:nb2,collapse=","),")",sep="")
   b1degrees <- summary(as.formula(paste('g ~ b1degree(',mesp,')',sep="")))
   mesp <- paste("c(",paste(0:nb1,collapse=","),")",sep="")
   b2degrees <- summary(as.formula(paste('g ~ b2degree(',mesp,')',sep="")))
   names(b2degrees) <- 0:nb1
   if(!is.null(b2degrees) & print){
    cat("Bipartite mode 2 degree distribution:\n")
    if(any(b2degrees>0)){print(b2degrees[b2degrees>0])}
   }
   names(b1degrees) <- 0:nb2
   if(!is.null(b1degrees) & print){
    cat("Bipartite mode 1 degree distribution:\n")
    if(any(b1degrees>0)){print(b1degrees[b1degrees>0])}
   }
   degrees <- list(b2=b2degrees, b1=b1degrees)
  }else{              
   mesp <- paste("c(",paste(0:(network.size(g)-1),collapse=","),")",sep="")
   degrees <- summary(as.formula(paste('g ~ degree(',mesp,')',sep="")))
   degrees <- degrees[degrees > 0]
   if(!is.null(degrees) & print){print(degrees)}
  }
 }
 invisible(degrees)
}


###############################################################################
# The <degreedistfactor> function returns the cross table of the degree
# distribution for a network and a given factor
###############################################################################

degreedistfactor <- function(g,x)
{
 if(!is.network(g)){
  stop("degreedist() requires a network object")
 }
 x <- get.vertex.attribute(g,x)
 degrees <- as.edgelist(g)
 if(length(degrees)>0){
  if(is.directed(g)){
   outdegrees <- table(degrees[,1],x[degrees[,2]])
#  outdegrees <- c(rep(0, network.size(g)-nrow(outdegrees)), outdegrees)
   if(!is.null(outdegrees)){print(table(outdegrees[,1]))}
   if(!is.null(outdegrees)){print(table(outdegrees[,2]))}
   indegrees <- table(degrees[,2],x[degrees[,1]])
#  indegrees <- c(rep(0, network.size(g)-nrow(indegrees)), indegrees)
   if(!is.null(indegrees)){print(table(indegrees[,1]))}
   if(!is.null(indegrees)){print(table(indegrees[,2]))}
   degrees <- list(indegrees=indegrees, outdegrees=outdegrees)
#  degrees <- rbind(indegrees, outdegrees)
  }else{
   degrees <- table(degrees,x[degrees])
   degrees <- c(rep(0, network.size(g)-nrow(degrees)), degrees)
   if(!is.null(degrees)){print(table(degrees))}
  }
 }
 invisible(degrees)
}


#espartnerdist <- function(g, print=TRUE)
#{
# if(!is.network(g)){
#  stop("espartnerdist() requires a network object")
# }
## twopaths
# smatrix <- g[,]
##twopaths <- crossprod(smatrix)
# twopaths <- t(smatrix) %*% smatrix
## (transitive) three-triangle (= triangle for directed)
# twopaths[smatrix==0] <- 0
# degrees <- tabulate(twopaths[row(twopaths)<col(twopaths)]+1,
#                     nbins=network.size(g)-1)
# if(is.directed(g)){
#  degrees[1] <- sum(smatrix)-sum(degrees[-1])
# }else{
#  degrees[1] <- 0.5*sum(smatrix)-sum(degrees[-1])
# }
# names(degrees) <- paste(0:(network.size(g)-2))
# if(print){
#  cat("ESP (edgewise shared partner) distribution:\n")
#  if(!is.null(degrees)){print(degrees[degrees>0])}
# }
# invisible(degrees)
#}


espartnerdist <- function(g, print=TRUE)
{
 if(!is.network(g)){
  stop("espartnerdist() requires a network object")
 }
 mesp <- paste("c(",paste(0:(network.size(g)-2),collapse=","),")",sep="")
 degrees <- summary(as.formula(paste('g ~ esp(',mesp,')',sep="")))
#names(degrees) <- paste(0:(network.size(g)-2))
 if(print){
  cat("ESP (edgewise shared partner) distribution:\n")
  if(!is.null(degrees)){print(degrees[degrees>0])}
 }
 invisible(degrees)
}


#dspartnerdist <- function(g, print=TRUE)
#{
# if(!is.network(g)){
#  stop("dspartnerdist() requires a network object")
# }
## twopaths
# smatrix <- g[,]
##twopaths <- crossprod(smatrix)
# twopaths <- t(smatrix) %*% smatrix
## (transitive) three-triangle (= triangle for directed)
# degrees <- tabulate(twopaths[row(twopaths)<col(twopaths)]+1,
#                     nbins=network.size(g)-1)
# names(degrees) <- paste(0:(network.size(g)-2))
# if(print){
#  cat("DSP (dyadwise shared partner) distribution:\n")
#  if(!is.null(degrees)){print(degrees[degrees>0])}
# }
# invisible(degrees)
#}



dspartnerdist <- function(g, print=TRUE)
{
 if(!is.network(g)){
  stop("dspartnerdist() requires a network object")
 }
 mesp <- paste("c(",paste(0:(network.size(g)-2),collapse=","),")",sep="")
 degrees <- summary(as.formula(paste('g ~ dsp(',mesp,')',sep="")))
#names(degrees) <- paste(0:(network.size(g)-2))
 if(print){
  cat("DSP (dyadwise shared partner) distribution:\n")
  if(!is.null(degrees)){print(degrees[degrees>0])}
 }
 invisible(degrees)
}



twopathdist <- function(g, print=TRUE)
{
 if(!is.network(g)){
  stop("twopathdist() requires a network object")
 }
# twopaths
 smatrix <- as.sociomatrix(g)
#twopaths <- crossprod(smatrix)
 twopaths <- t(smatrix) %*% smatrix
# (transitive) three-triangle (= triangle for directed)
#twopaths[smatrix==0] <- 0
 degrees <- tabulate(twopaths[row(twopaths)<col(twopaths)]+1,
                     nbins=network.size(g)-1)
 if(is.directed(g)){
  degrees[1] <- sum(smatrix)-sum(degrees[-1])
 }else{
  degrees[1] <- 0.5*sum(smatrix)-sum(degrees[-1])
 }
 names(degrees) <- paste(0:(network.size(g)-2))
 if(print){
  cat("twopath distribution:\n")
  if(!is.null(degrees)){print(degrees[degrees>0])}
 }
 invisible(degrees)
}


"rspartnerdist" <- function (g, print = TRUE) 
{
    if (!is.network(g)) {
        stop("rspartnerdist() requires a network object")
    }
    ds <- dspartnerdist(g,print=FALSE)
    es <- espartnerdist(g,print=FALSE)
    vec <- 0: (network.size(g) - 2)
    upperlimit <- max(vec[ds>0]) + 1
    rs <- es[1:upperlimit]/ds[1:upperlimit]

    if (print) {
        cat("ESP/DSP ratio distribution:\n")
        print(rs)
    }
    invisible(rs)
}


summary.statsmatrix.ergm <- function(object, ...){
 c(summary(round(object,digits=8), ...),
   round(ergm.t.summary(object),5))
}


###############################################################################
# The <ergm.t.summary> function conducts a t test for comparing the mean of a
# given vector and a hypothesized mean
###############################################################################

ergm.t.summary <-
function(x, alternative = c("two.sided", "less", "greater"),
         mu = 0, var.equal = FALSE, conf.level = 0.95,
         ...)
{
    alternative <- match.arg(alternative)

    if(!missing(mu) && (length(mu) != 1 || is.na(mu)))
        stop("mu must be a single number")
    if(!missing(conf.level) &&
       (length(conf.level) != 1 || !is.finite(conf.level) ||
        conf.level < 0 || conf.level > 1))
        stop("conf.level must be a single number between 0 and 1")

    dname <- deparse(substitute(x))
    xok <- !is.na(x)
    yok <- NULL
    x <- x[xok]
    nx <- length(x)
    mx <- mean(x)
    rval <- c(NA, mx, NA, NA)
    names(rval) <- c("std. err.","std. t units.","p.value dev.","std. dev.")
    if(nx < 2){
#     stop("not enough x observations")
      return(rval)
    }
    vx <- var(x)
    estimate <- mx
    df <- nx-1
    stddev <- sqrt(vx)
    stderr <- sqrt(vx/nx)
    if(stderr < 10 *.Machine$double.eps * abs(mx)){
#       stop("data are essentially constant")
        return(rval)
    }
    tstat <- (mx-mu)/stderr
    method <- "One Sample t-test"
    names(estimate) <- "mean of x"
    if (alternative == "less") {
	pval <- pt(tstat, df)
    }
    else if (alternative == "greater") {
	pval <- pt(tstat, df, lower.tail = FALSE)
    }
    else {
	pval <- 2 * pt(-abs(tstat), df)
    }
    rval[1:4] <- c(stderr, tstat, pval, stddev)
    return(rval)
}

newnw.extract<-function(oldnw,z,output="network"){
  if("newedgelist" %in% names(z)){
    newedgelist<-z$newedgelist[,1:2,drop=FALSE]
  }else{
    nedges<-z$newnwtails[1]
    # *** don't forget - edgelists are cbind(tails, heads) now
    newedgelist <-
      if(nedges>0) cbind(z$newnwtails[2:(nedges+1)],z$newnwheads[2:(nedges+1)])
      else matrix(0, ncol=2, nrow=0)
  }
  
  network.update(oldnw,newedgelist,matrix.type="edgelist",output=output)
}

nvattr.copy.network <- function(to, from, ignore=c("bipartite","directed","hyper","loops","mnext","multiple","n")){
  for(a in list.vertex.attributes(from)){
    if(! a%in%ignore)
      to <- set.vertex.attribute(to, a, get.vertex.attribute(from, a, unlist=FALSE))
  }
  for(a in list.network.attributes(from)){
    if(! a%in%ignore)
      to <- set.network.attribute(to, a, get.network.attribute(from, a, unlist=FALSE))
  }
  to
}


statnet.edit <- function(name,package=c("statnet","ergm","network")){
  i <- 1
  while(i < length(package)){
   pkgpath <- .find.package(package[i],quiet=TRUE)
   if(length(pkgpath)>0){
    filepath <- file.path(pkgpath,name)
    if(file.exists(filepath)){
     i <- length(package)+1
     file.edit(filepath)
    }
   }
   i <- i + 1
  }
  if(i != length(package)+2){
   warning(paste("The file '",name,"' does not seem to exist in 'statnet'.",
           sep=""), call. = FALSE)
  }
  invisible(filepath)
}

## Compress a data frame by eliminating duplicate rows while keeping
## track of their frequency.
compress.data.frame<-function(x){
  x<-sort(x)
  firsts<-which(!duplicated(x))
  freqs<-diff(c(firsts,nrow(x)+1))
  x<-x[firsts,]
  list(rows=x,frequencies=freqs)
}

## Sorts rows of a data frame in lexicographic order.
sort.data.frame<-function(x, decreasing=FALSE, ...){
  x[do.call(order,c(sapply(seq_along(x),function(i)x[[i]],simplify=FALSE), decreasing=decreasing)),]
}

## Concatenate a character list with commas and ands in the right places.
paste.and <- function(x, oq='', cq=''){
  x <- paste(oq, x, cq, sep='')
  if(length(x)==0) return('')
  if(length(x)==1) return(x)
  if(length(x)==2) return(paste(x[1],'and',x[2]))
  if(length(x)>=3) return(paste(paste(x[-length(x)], collapse=", "),', and ',x[length(x)],sep=''))
}
