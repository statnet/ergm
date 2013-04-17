#  File R/ergm.utility.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
#######################################################################
#==============================================================
# This file contains the following 21 utility functions:
#      <ostar2deg>                  
#      <is.invertible>          <summary.statsmatrix.ergm>
#      <is.ergm>                <ergm.t.summary>
#      <is.latent>
#      <degreedist>             <is.latent.cluster>
#      <degreedistfactor>       <newnw.extract>
#      <espartnerdist>          <dspartnerdist>         
#      <rspartnerdist>         
#      <twopathdist>            <copy.named>
#      <compress.data.frame>    <sort.data.frame>
#==============================================================      



###############################################################################
# The <ostar2deg> function ??
#
# --PARAMETERS--
#   object  : an ergm object
#   ninflast: whether ??
#
# --RETURNED--
#   odeg: the vector of ??
#
###############################################################################

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
#
# --PARAMETERS--
#   g    : a network object
#   print: whether to print the degree distribution; default=TRUE
#
# --RETURNED--
#   degrees:
#      if directed  -- a matrix of the distributions of in and out degrees;
#                      this is row bound and only contains degrees for which
#                      one of the in or out distributions has a positive count
#      if bipartite -- a list containing the degree distributions of b1 and b2
#      otherwise    -- a vector of the positive values in the degree
#                      distribution
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
#
# --PARAMETERS--
#   g: a network
#   x: a nodal attribute, as a character string
#
# --RETURNED--
#   degrees:
#      if directed  -- a list containing 2 cross tables, the in degree
#                      distributions by 'x', and out degree dist by 'x'
#      otherwise    -- a table of the degree distribution by 'x'
#
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
#
# --PARAMETERS--
#   x          : a numeric vector
#   alternative: a string to indicate whether the test is two-sided or one-
#                sided to the left or right, as either "two.sided", "less",
#                or "greater"; default="two.sided"
#   mu         : the hypothesized mean; default = 0
#
# --IGNORED PARAMETERS--
#   var.equal : whether the variance of ?? is ??; default=FALSE
#   conf.level: the confidence level; default=0.95
#   ...       : ??
#
# --RETURNED--
#   rval: a vetor of the standard error, the t statistic, the p value, and the
#         standard deviation, consistent with 'alternative'; if the length of
#         x is <2, this vector will be predominently NA's
#
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

newnw.extract<-function(oldnw,z,output="network",response=NULL){
  if("newedgelist" %in% names(z)){
    newedgelist<-z$newedgelist[,1:2,drop=FALSE]
    if(!is.null(response))
       newnwweights<-z$newedgelist[,3]
  }else{
    nedges<-z$newnwtails[1]
    # *** don't forget - edgelists are cbind(tails, heads) now
    newedgelist <-
      if(nedges>0) cbind(z$newnwtails[2:(nedges+1)],z$newnwheads[2:(nedges+1)])
      else matrix(0, ncol=2, nrow=0)
    newnwweights <- z$newnwweights[2:(nedges+1)]
  }
  
  newnw<-network.update(oldnw,newedgelist,matrix.type="edgelist",output=output)
  if(!is.null(response)){
    # It's very important that the order of weights here is the same
    # as the one that network accepts.
    newnw<-set.edge.attribute(newnw,attrname=response,newnwweights,e=apply(newedgelist,1,function(e) get.edgeIDs(newnw,e[1],e[2])))
  }
  newnw
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


# Create a copy of a network of interest with certain guarantees about its internal representation:
# * tails < heads
# * no (tail,head) pair has more than one edge ID associated with it
standardize.network <- function(nw, preserve.eattr=TRUE){
  if(preserve.eattr){
    el <- rbind(as.edgelist(nw),as.edgelist(is.na(nw)))
    eids <- apply(el, 1, function(e) get.edgeIDs(nw, e[1], e[2], na.omit=FALSE))
    if(is.list(eids)){ # I.e., apply() wasn't able to simplify.
      bad.ei <- which(sapply(eids,length)>1)
      for(ei in bad.ei){
        dup.eids <- duplicated(nw$mel[eids[[ei]]])
        if(sum(!dup.eids)!=1) stop("Edge (",el[ei,1],",",el[ei,2],") has multiple IDs with distinct attributes. Cannot repair.")
        eids[[ei]] <- eids[[ei]][!dup.eids]
      }
      eids <- unlist(eids)
    }
    
    vals <- lapply(nw$mel,"[[","atl")[eids]
    names <- lapply(vals, names)
    el.na <- NULL
  }else{
    el <- rbind(as.edgelist(nw))
    vals <- NULL
    names <- NULL
    el.na <- as.edgelist(is.na(nw))
  }
  
  nw <- delete.edges(nw, seq_along(nw$mel))
  nw <- add.edges(nw, el[,1], el[,2], names.eval=names, vals.eval=vals)
  if(!is.null(el.na)) nw[el.na] <- NA
  nw
}

get.free.dyads <- function(constraints){
  y <- NULL
  for(con in constraints){
    if(!is.null(con$free.dyads)){
      y <- if(is.null(y)) standardize.network(con$free.dyads(),FALSE) else y & standardize.network(con$free.dyads(),FALSE)
    }
  }
  y
}

get.miss.dyads <- function(constraints, constraints.obs){
  free.dyads <- get.free.dyads(constraints)
  free.dyads.obs <- get.free.dyads(constraints.obs)
  
  if(is.null(free.dyads)){
    if(is.null(free.dyads.obs)) NULL
    else free.dyads.obs
  }else{
    if(is.null(free.dyads.obs)) standardize.network(invert.network(free.dyads),FALSE)
    else standardize.network(invert.network(free.dyads),FALSE) | free.dyads.obs
  }
}

.hash.el <- function(x){
  apply(x, 1, paste, collapse="\r")
}

invert.network <- function(nw){
  n <- network.size(nw)
  m <- nw %n% "bipartite"

  # Borrows from !.network:
  el <- cbind(rep(1:n, each = n), rep(1:n, n))
  if (!is.directed(nw)) 
    el <- el[el[, 1] <= el[, 2], , drop=FALSE]
  if (!has.loops(nw)) 
    el <- el[el[, 1] != el[, 2], , drop=FALSE]

  if(m) el <- el[(el[,1]<=m)!=(el[,2]<=m), , drop=FALSE]

  # el now contains the set of possible dyads

  el <- el[! .hash.el(el) %in% .hash.el(as.edgelist(nw)), , drop=FALSE]

  network.update(nw, el, matrix.type="edgelist")
  
}


# Return the name of the package containing function f visible from
# environment env.
which.package.InitFunction <- function(f, env = parent.frame()){
  f <- as.character(f)
  # Find the first entity named f in the search path, and get its name
  # attribute (if present).
  loc <- attr(findFunction(f, where=env)[[1]], "name")
  # If name attribute is not NULL and begins with "package:", return
  # the package name. Otherwise, return NULL.
  if(!is.null(loc) && grepl("^package:", loc)) sub("^package:", "", loc) else NULL
}
