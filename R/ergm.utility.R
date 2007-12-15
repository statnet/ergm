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
    ev <- eigen(V, sym = TRUE, only.values = TRUE)$values
    all(ev/max(ev) > tol)
}
is.ergm <- function(object)
{
    class(object)=="ergm"
}
is.matrixnetwork<-function(x){
 is.matrix(x)|is.network(x)
}
degreedist <- function(g, print=TRUE)
{
 if(!is.network(g)){
  stop("degreedist() requires a network object")
 }
 if(is.directed(g)){                                      
   mesp <- paste("c(",paste(0:(network.size(g)-1),collapse=","),")",sep="")
   outdegrees <- summary(as.formula(paste('g ~ odegree(',mesp,')',sep="")),drop=FALSE)
   indegrees <- summary(as.formula(paste('g ~ idegree(',mesp,')',sep="")),drop=FALSE)
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
   adegrees <- summary(as.formula(paste('g ~ adegree(',mesp,')',sep="")),drop=FALSE)
   mesp <- paste("c(",paste(0:nb1,collapse=","),")",sep="")
   edegrees <- summary(as.formula(paste('g ~ edegree(',mesp,')',sep="")),drop=FALSE)
   names(edegrees) <- 0:nb1
   if(!is.null(edegrees) & print){
    cat("Event degree distribution:\n")
    if(any(edegrees>0)){print(edegrees[edegrees>0])}
   }
   names(adegrees) <- 0:nb2
   if(!is.null(adegrees) & print){
    cat("Actor degree distribution:\n")
    if(any(adegrees>0)){print(adegrees[adegrees>0])}
   }
   degrees <- list(b2=edegrees, b1=adegrees)
  }else{              
   mesp <- paste("c(",paste(0:(network.size(g)-1),collapse=","),")",sep="")
   degrees <- summary(as.formula(paste('g ~ degree(',mesp,')',sep="")),drop=FALSE)
   degrees <- degrees[degrees > 0]
   if(!is.null(degrees) & print){print(degrees)}
  }
 }
 invisible(degrees)
}
degreedistfactor <- function(g,x)
{
 if(!is.network(g)){
  stop("degreedist() requires a network object")
 }
 x <- get.vertex.attribute(g,x)
 degrees <- as.matrix.network.edgelist(g)
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
 degrees <- summary(as.formula(paste('g ~ esp(',mesp,')',sep="")),drop=FALSE)
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
 degrees <- summary(as.formula(paste('g ~ dsp(',mesp,')',sep="")),drop=FALSE)
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
    ds <- dspartnerdist(g,print=F)
    es <- espartnerdist(g,print=F)
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
	pval <- pt(tstat, df, lower = FALSE)
    }
    else {
	pval <- 2 * pt(-abs(tstat), df)
    }
    rval[1:4] <- c(stderr, tstat, pval, stddev)
    return(rval)
}

# Rewritten version of mixingmatrix that does not use the nodemix model term nor
# return totals:
mixingmatrix <- function(nw, attrname) {
  if(!is.network(nw)){
    stop("mixingmatrix() requires a network object")
  }                                                                
  nodecov <- get.node.attr(nw, attrname)
  u<-sort(unique(nodecov))
  # nodecovnum <- match(nodecov, u)
  el <- as.matrix.network.edgelist(nw)
  type <- "directed"
  if (is.bipartite(nw)) { # must have heads < tails now
    if (is.directed(nw)) 
      cat("Warning!  Bipartite networks are currently\n",
          "automatically treated as undirected\n")
    type <- "bipartite"
    rowswitch <- apply(el, 1, function(x) x[1]>x[2])
    el[rowswitch, 1:2] <- el[rowswitch, 2:1]
  }
  From <- c(u, nodecov[el[,1]])
  To <- c(u, nodecov[el[,2]])
  tabu <- table(From, To)  # Add u,u diagonal to ensure each 
  # value is represented, then subtract it later
  diag(tabu) <- diag(tabu) - 1
  if(!is.directed(nw) && !is.bipartite(nw)){
    type <- "undirected"
    tabu <- tabu + t(tabu)
    diag(tabu) <- diag(tabu)/2
  }
  ans <- list(type=type, matrix=tabu)
  class(ans) <- "mixingmatrix"
  ans
}

print.mixingmatrix <- function(x, ...) {
  m <- x$mat
  rn <- rownames(m)
  cn <- colnames(m)  
  if (x$type == "undirected") {
    dimnames(m) <- list(rn, cn)
    cat("Note:  Marginal totals can be misleading\n",
        "for undirected mixing matrices.\n")
  } else {
    total <- apply(m,1,sum)
    m <- cbind(m,total)
    total <- apply(m,2,sum)
    m <- rbind(m,total)
    rn <- c(rn, "Total")
    cn <- c(cn, "Total")
    if (x$type == "bipartite")
      dimnames(m) <- list(Actors = rn,Events = cn)
    else
      dimnames(m) <- list(From = rn,To = cn)
  }
  print(m)
}


## Removed old mixingmatrix.via.edgelist function (DH, Oct. 17 2007)
## There is now only one mixingmatrix function, and it does not
## depend on nodemix.  (NB:  Nodemix sometimes depends on mixingmatrix,
## so it's good to avoid circularity)
#mixingmatrix.via.edgelist <- function(g,attrname)
#{
# if(!is.network(g)){
#  stop("mixingmatrix() requires a network object")
# }
# nodecov <- get.node.attr(g, attrname)
# u<-sort(unique(nodecov))
# if(any(is.na(nodecov))){u<-c(u,NA)}
# nodecov <- match(nodecov,u,nomatch=length(u)+1)
# if(any(is.na(u))){u <- paste(u)}
# ui<-seq(along=u)
##
##   Check for matches
##
# degrees <- as.matrix.network.edgelist(g)
# tabu <- table(c(nodecov[degrees[,1]],ui),c(nodecov[degrees[,2]],ui)) 
# diag(tabu) <- diag(tabu)-1
# if(!is.directed(g)  & !is.bipartite(g)){
#   tabu <- tabu + t(tabu)
#   total <- apply(tabu,1,sum)
#   tabu <- cbind(tabu,total)
#   tabu <- rbind(tabu,c(total,sum(total)))
#   tabu[row(tabu)>col(tabu)] <- 0
#   diag(tabu) <- diag(tabu)/2
#   tabu[nrow(tabu),] <- c(total,sum(total))
##  dimnames(tabu)[[1]][nrow(tabu)] <- "total"
#   dimnames(tabu) <- list(c(u,"total"),c(u,"total"))
# }else{
#   total <- apply(tabu,1,sum)
#   tabu <- cbind(tabu,total)
#   total <- apply(tabu,2,sum)
#   tabu <- rbind(tabu,total)
#   dimnames(tabu) <- list(c(u,"total"),c(u,"total"))
##  if(is.bipartite(g)){tabu <- t(tabu)}
# }
# tabu
#}


## Eliminated old mixingmatrix function (DH, Oct. 17 2007)
## because it depends on nodemix.  Replaced by 
## modified version of old mixingmatrix.via.edgelist
#mixingmatrix <- function(g,attrname)
#{
# if(!is.network(g)){
#  stop("mixingmatrix() requires a network object")
# }
# nodecov <- get.node.attr(g, attrname)
# u<-sort(unique(nodecov))
# su <- summary(as.formula(paste("g~nodemix('",attrname,"')",sep="")))
# if(!is.directed(g) & !is.bipartite(g)){
#   tabu <- matrix(0,ncol=length(u),nrow=length(u))
#   tabu[row(tabu)<=col(tabu)] <- su
#   tabu <- t(tabu)
#   tabu[row(tabu)<=col(tabu)] <- su   
##   tabu[row(tabu)>=col(tabu)] <- su
##   tabu <- tabu + t(tabu) - diag(diag(tabu))
#   total <- apply(tabu,1,sum)
#   tabu <- cbind(tabu,total)
#   tabu <- rbind(tabu,c(total,sum(total)))
#   tabu[row(tabu)>col(tabu)] <- 0
##  diag(tabu) <- diag(tabu)/2
#   tabu[nrow(tabu),] <- c(total,sum(total))
##  dimnames(tabu)[[1]][nrow(tabu)] <- "total"
#   dimnames(tabu) <- list(c(u,"total"),c(u,"total"))
# }else{
#   tabu <- matrix(su, ncol=length(u))
#   total <- apply(tabu,1,sum)
#   tabu <- cbind(tabu,total)
#   total <- apply(tabu,2,sum)
#   tabu <- rbind(tabu,total)
#   dimnames(tabu) <- list(c(u,"total"),c(u,"total"))
##  if(is.bipartite(g)){tabu <- t(tabu)}
# }
# tabu
#}



#
# Return TRUE iff object x is a latentfit object
# or a latent model
#
is.latent<-function(x){
  out <- is.list(x)
  if(out){
   out <- is.logical(x$latent)
   if(out){
    out <- x$latent
   }
  }
  out
}
#
# Return TRUE iff object x is a latentfit object
# or a latent model and a "latentcluster" model or fit
#
is.latent.cluster<-function(x){
  out <- FALSE
  if(is.latent(x)){
   out <- is.list(x)
   if(out){
    out <- is.logical(x$cluster)
    if(out){
     out <- x$cluster
    }
   }
  }
  out
}

newnw.extract<-function(oldnw,z){
  nedges<-z$newnwheads[1]
  newedgelist <-
    if(nedges>0) cbind(z$newnwheads[2:(nedges+1)],z$newnwtails[2:(nedges+1)])
    else matrix(0, ncol=2, nrow=0)
  
  network.update(oldnw,newedgelist,"edgelist")
}
