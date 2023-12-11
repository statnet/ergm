#  File R/ergm.utility.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2023 Statnet Commons
################################################################################

#' @rdname ergm
#' @importFrom methods is
#' @export
is.ergm <- function(object)
{
    is(object,"ergm")
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



#' Computes and Returns the Degree Distribution Information for a Given Network
#' 
#' The \code{degreedist} generic computes and returns the degree
#' distribution (number of vertices in the network with each degree
#' value) for a given network. This help page documents the
#' function. For help about [the ERGM sample space constraint with
#' that name][degreedist-ergmConstraint], try
#' `help("degreedist-constraint")`.
#' 
#' @param object a \code{network} object or some other object for
#'   which degree distribution is meaningful.
#' @param \dots Additional arguments to functions.
#' @return If directed, a matrix of the distributions of in and out
#'   degrees; this is row bound and only contains degrees for which
#'   one of the in or out distributions has a positive count.  If
#'   bipartite, a list containing the degree distributions of b1 and
#'   b2.  Otherwise, a vector of the positive values in the degree
#'   distribution
#' @examples
#' 
#' data(faux.mesa.high)
#' degreedist(faux.mesa.high)
#' 
#' @export
degreedist <- function(object, ...) UseMethod("degreedist")

#' @describeIn degreedist Method for [`network`] objects.
#' @param print logical, whether to print the degree distribution.
#' @export
degreedist.network <- function(object, print=TRUE, ...)
{
 g <- object
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

#' Copy network- and vertex-level attributes between two network objects
#' 
#' An internal ergm utility function to copy the network-level attributes and
#' vertex-level attributes from one \code{\link{network}} object to another,
#' ignoring some standard properties by default.
#' 
#' 
#' @param to the \code{\link{network}} that attributes should be copied to
#' @param from the \code{\link{network}} that attributes should be copied to
#' @param ignore vector of charcter names of network attributes that should not
#' be copied. Default is the standard list of network properties created by
#' \code{\link{network.initialize}}
#' @return returns the \code{to} network, with attributes copied from
#' \code{from}
#' @note does not check that networks are of the same size, etc
#' @seealso \code{\link{set.vertex.attribute}},
#' \code{\link{set.network.attribute}}
#' @keywords internal
#' @export nvattr.copy.network
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

#' @importFrom tibble tibble
single.impute.dyads <- function(nw, constraints=NULL, constraints.obs=NULL, min_informative=NULL, default_density=NULL, output=c("network","ergm_state"), verbose=FALSE){
  output <- match.arg(output)
  response <- nw %ergmlhs% "response"
  stopifnot(!is.null(constraints)||is.null(constraints.obs))

  if(!is.null(constraints)){
    imputable <- as.rlebdm(constraints, constraints.obs, "missing")
    nae <- NVL3(imputable, sum(.), 0)
    if(nae) na.el <- as.edgelist(imputable, output="tibble") # FIXME: Avoid creating edgelists.
  }else{
    nae <- network.naedgecount(nw)
    if(nae) na.el <- as.edgelist(is.na(nw), output="tibble")
  }
  if(nae==0){
    if(output=="network") return(nw)
    else return(ergm_state(nw))
  }

  if(verbose) message("Imputing ", nae, " dyads is required.")

  el2s <- function(el) apply(el, 1, paste, collapse=",")
  s2el <- function(s) matrix(as.integer(do.call(rbind,strsplit(s,","))),ncol=2)

  min_informative <- NVL3(min_informative, if(is.function(.)) .(nw) else ., 0)
  default_density <- if(is.function(default_density)) default_density(nw)

  if(!is.null(constraints)){ # Constraints
    informative <- as.rlebdm(constraints, constraints.obs, "informative")
    nonzeros <- as.rlebdm(nw)
    if(!is.valued(nw)){
      d <-
        if(sum(informative)<min_informative){
          message("Number of informative dyads is too low. Using default imputation density.")
          default_density
        }else sum(nonzeros & informative)/sum(informative)
      nimpute <- round(d*nae)
    }else{
      if(sum(informative)<min_informative){
        message("Number of informative dyads is too low. Imputing valued relations is not possible.")
        return(nw)
      }
      x <- as.edgelist(nw, attrname=response, output="tibble")
      x <- x[[3L]][! el2s(x[1:2])%in%el2s(na.el)]
      zeros <- sum(informative) - length(x)
    }
  }else{ # No Constraints
    if(!is.valued(nw)){
      d <-
        if(network.dyadcount(nw,na.omit=TRUE)<min_informative){
          message("Number of informative dyads is too low. Using default imputation density.")
          default_density
        }else network.edgecount(nw,na.omit=TRUE)/network.dyadcount(nw,na.omit=TRUE)
      nimpute <- round(d*nae)
    }else{
      if(network.dyadcount(nw,na.omit=TRUE)<min_informative){
        message("Number of informative dyads is too low. Imputing valued relations is not possible.")
        return(nw)
      }
      x <- as.edgelist(nw, attrname=response, output="tibble")[[3]]
      zeros <- network.dyadcount(nw,na.omit=TRUE) - length(x)
    }
  }
  
  if(!is.valued(nw)){
    if(verbose) message("Imputing ", nimpute, " edges at random.")
    i.new <- sample.int(nae,nimpute)
    if(output=="network"){
      y.cur <- nw[na.el]
      i.na <- which(is.na(y.cur))
      i.cur <- which(y.cur!=0)
      todel <- union(setdiff(i.cur, i.new), setdiff(i.na, i.new))
      toadd <- union(setdiff(i.new, i.cur), intersect(i.na, i.new))
      nw[na.el[c(todel,toadd),,drop=FALSE]] <- rep(0:1, c(length(todel),length(toadd)))
    }else{ # edgelist
      el <- s2el(union(setdiff(el2s(as.edgelist(nw, output="tibble")), el2s(na.el)), el2s(na.el[i.new,,drop=FALSE])))
      colnames(el) <- c(".tail",".head")
      nw <- ergm_state(el, nw=nw)
    }
  }else{
    if(output=="network"){
      nw[na.el,names.eval=response,add.edges=TRUE] <- sample(c(0,x),nae,replace=TRUE,prob=c(zeros,rep(1,length(x))))
    }else{ # edgelist
      el <- as.edgelist(nw, attrname=response, output="tibble")
      el <- el[!el2s(el[1:2])%in%el2s(na.el),,drop=FALSE]
      na.el <- cbind(na.el, sample(c(0,x),nae,replace=TRUE,prob=c(zeros,rep(1,length(x)))))
      colnames(na.el)[3] <- response
      el <- rbind(el, na.el)
      el <- el[el[[3L]]!=0,,drop=FALSE]
      colnames(el) <- c(".tail",".head",response)
      nw <- ergm_state(el, nw=nw)
    }
  }

  nw
}

## A is a matrix. V is a column vector that may contain Infs
## computes A %*% V, counting 0*Inf as 0
.multiply.with.inf <- function(A,V) {
  cbind(colSums(t(A)*c(V), na.rm=TRUE))
}

trim_env_const_formula <- function(x, keep=NULL){
  terms <- list_rhs.formula(x)
  needs_env <- FALSE
  basefn <- ls(baseenv())

  is_simple <- function(x, toplevel=FALSE){
    if(is.null(x) || is.character(x) || is.numeric(x) || is.logical(x) || is.complex(x) || identical(x, as.name("."))) TRUE
    else if(is.call(x)){
      fn <- as.character(x[[1]])
      if(!toplevel && ! fn%in%basefn) FALSE
      else all(sapply(as.list(x)[-1], is_simple))
    } else FALSE
  }

  for(trm in terms){
    if(is.call(trm) && !is_simple(trm, TRUE)){
        needs_env <- TRUE
        break
    }
  }

  if(needs_env) x else trim_env(x, keep)
}

### A patched version of statnet.common::sginv() that uses
### eigendecomposition rather than the SVD for the case when symmetric
### non-negative-definite (snnd) is TRUE.
###
### TODO: Delete after incorporation into statnet.common.
sginv <- function(X, tol = sqrt(.Machine$double.eps), ..., snnd = TRUE){
  if(!snnd) statnet.common::sginv(X, ..., snnd)

  d <- diag(as.matrix(X))
  d <- ifelse(d==0, 1, 1/d)

  d <- sqrt(d)
  if(anyNA(d)) stop("Matrix a has negative elements on the diagonal, and snnd=TRUE.")
  dd <- rep(d, each = length(d)) * d
  X <- X * dd

  ## Perform inverse via eigendecomposition, removing too-small eigenvalues.
  e <- eigen(X, symmetric=TRUE)
  keep <- e$values > max(tol * e$values[1L], 0)
  e$vectors[, keep, drop=FALSE] %*% ((1/e$values[keep]) * t(e$vectors[, keep, drop=FALSE])) * dd
}

#' Evaluate a quadratic form with a possibly singular matrix using
#' eigendecomposition after scaling to correlation.
#'
#' Let \eqn{A} be the matrix of interest, and let \eqn{D} is a
#' diagonal matrix whose diagonal is same as that of \eqn{A}.
#'
#' Let \eqn{R = D^{-1/2} A D^{-1/2}}. Then \eqn{A = D^{1/2} R D^{1/2}} and
#' \eqn{A^{-1} = D^{-1/2} R^{-1} D^{-1/2}}.
#'
#' Decompose \eqn{R = P L P'} for \eqn{L} diagonal matrix of eigenvalues
#' and \eqn{P} orthogonal. Then \eqn{R^{-1} = P L^{-1} P'}.
#'
#' Substituting,
#' \deqn{x' A^{-1} x  = x' D^{-1/2} P L^{-1} P' D^{-1/2} x = h' L^{-1} h}
#' for \eqn{h = P' D^{-1/2} x}.
#'
#' @noRd
xTAx_seigen <- function(x, A, tol=sqrt(.Machine$double.eps), ...) {
  d <- diag(as.matrix(A))
  d <- ifelse(d<=0, 0, 1/d)

  d <- sqrt(d)
  dd <- rep(d, each = length(d)) * d
  A <- A * dd

  e <- eigen(A)

  keep <- e$values > max(tol * e$values[1L], 0)

  h <- drop(crossprod(x*d, e$vectors[,keep,drop=FALSE]))

  structure(sum(h*h/e$values[keep]), rank = sum(keep), nullity = sum(!keep))
}
