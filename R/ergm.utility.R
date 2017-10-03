#  File R/ergm.utility.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
#==============================================================
# This file contains the following 21 utility functions:
#      <ostar2deg>                  
#      <is.invertible>          <summary.statsmatrix.ergm>
#      <is.ergm>                <ergm.t.summary>
#      <is.latent>
#      <degreedist>             <is.latent.cluster>
#      <newnw.extract>
#      <espartnerdist>          <dspartnerdist>         
#      <rspartnerdist>         
#      <twopathdist>            <copy.named>
#      <compress.data.frame>    <sort.data.frame>
#      <catchToList>
#==============================================================      

#' @rdname ergm
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
#' The \code{degreedist} generic computes and returns the degree distribution
#' (number of vertices in the network with each degree value) for a given
#' network.
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

# generate a network object from the edgelist output of the mcmc sample


#' Internal function to create a new network from the ergm MCMC sample output
#' 
#' An internal function to generate a new \code{\link{network}} object using
#' the output (lists of toggled heads and tail vertices) from an ERGM MCMC or
#' SAN process.
#' 
#' 
#' @param oldnw a network object (presumably input to the ergm process) from
#' which the network- and vertex-level attributes will be copied
#' @param z a list having either a component named \code{newedgelist} or two
#' components \code{newtails} and \code{newheads} containing the ids of the
#' head and tails vertices of the edges.  Optionall \code{newweights}
#' containing edgewights.
#' @param output passed to \code{\link{network.update}}, which claims not to
#' use it
#' @param response optional character string giving the name of the edge
#' attribute where the edge values (weight/count) should be stored.
#' @return a \code{\link{network}} object with properties copied from
#' \code{oldnw} and edges corresponding to the lists of tails and head vertex
#' ids in \code{z}
#' @note This is an internal ergm function, it most cases with edgelists to be
#' converted to networks it will probably be simpler to use
#' \code{\link{network.edgelist}}
#' @seealso \code{\link{network.edgelist}}, \code{\link{network.update}}
#' @export newnw.extract
newnw.extract<-function(oldnw,z,output="network",response=NULL){
  # if z has a newedgelist attached, use it
  if("newedgelist" %in% names(z)){
    newedgelist<-z$newedgelist[,1:2,drop=FALSE]
    if(!is.null(response))
       newnwweights<-z$newedgelist[,3]
  }else{
    # expect that z will have seperate lists of heads and tails
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



#' Copy a network object enforcing ergm-appropriate guarantees about its
#' internal representation
#' 
#' Create a copy of a \code{\link{network}} of interest with certain guarantees
#' about its internal representation: \itemize{ \item for every edge, the id of
#' the 'tails' vertex is < id of the 'heads' vertex if the network is
#' undirected \item no (tail,head) id pair has more than one edge ID associated
#' with it (no multiplex edges) }
#' 
#' This function is needed because the \code{\link{network}} object can support
#' added non-directed edges in arbitrary order, as well as multiplex edges and
#' hypergraphs (even if the network is not marked as such), which are not
#' supported in the underlying ergm processes. Uses \code{\link{as.edgelist}}
#' internally to make the conversion.
#' 
#' @param nw a \code{\link{network}} object to be copied
#' @param preserve.eattr logical; should the edge attributes be preserved
#' during the copying process (presumably slower)
#' @return returns an ergm-appropriate network object.
#' @note This function may be time expensive.  Also, the transformation is
#' performed by deleting edges in initial network and re-adding them.
#' @keywords internal
#' @examples
#' 
#' test<-network.initialize(5,directed=FALSE)
#' test[2,1]<-1  #ergm wont like this
#' test$mel[[1]] # peek at internal representation
#' 
#' test2<-standardize.network(test) # enforce!
#' test2$mel[[2]]  # 1 and 2 have traded places
#' 
#' 
#' @export standardize.network
standardize.network <- function(nw, preserve.eattr=TRUE){
  if(preserve.eattr){
    el <- rbind(as.edgelist(nw),as.edgelist(is.na(nw)))
    eids <- lapply(seq_len(nrow(el)), function(i) get.edgeIDs(nw, el[i,1], el[i,2], na.omit=FALSE))

    bad.ei <- which(sapply(eids,length)>1)
    for(ei in bad.ei){
      dup.eids <- duplicated(nw$mel[eids[[ei]]])
      if(sum(!dup.eids)!=1) stop("Edge (",el[ei,1],",",el[ei,2],") has multiple IDs with distinct attributes. Cannot repair.")
      eids[[ei]] <- eids[[ei]][!dup.eids]
    }
    eids <- unlist(eids)
    
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


.hash.el <- function(x){
  apply(x, 1, paste, collapse="\r")
}

locate.InitFunction <- function(name, prefix, errname=NULL, env = parent.frame()){
  if(is.call(name)) name <- name[[1]]
  name <- as.character(name)
  fname <- paste(prefix,name,sep=".")
  
  f <- try(get(fname, mode='function', envir=env), silent=TRUE)
  if(inherits(f, "try-error")){
    m <- getAnywhere(fname)
    if(length(m$objs)){
      ## Prioritise visible over not:
      if(any(m$visible)){
        m <- lapply(m[-1], "[", m$visible)
      }
      if(length(m$objs)>1) warning("Name ",fname," matched by multiple objects; using the first one on the list.")
      f <- m$objs[[1]] 
    }else{
      if(!is.null(errname)) stop(errname,' ', sQuote(name), " initialization function ", sQuote(fname), " not found.") else f <- NULL
    }
  }
  f
}

# Return the name of the package containing function f visible from
# environment env.
which.package.InitFunction <- function(f, env = parent.frame()){
  f <- as.character(f)
  # Find the first entity named f in the search path, and get its name
  # attribute (if present).
  found <- methods::findFunction(f, where=env)
  if (length(found) > 0) {
    loc <- attr(found[[1]], "name")
    # If name attribute is not NULL and begins with "package:", return
    # the package name. Otherwise, return NULL.
    if(!is.null(loc) && grepl("^package:", loc)) sub("^package:", "", loc) else NULL
  } else {
    # can't find the function normally; the package might have been imported
    # instead of attached
    found <- get(f, envir=env)
    environmentName(environment(found))
  }
  
}

single.impute.dyads <- function(nw, response=NULL){
    nae <- network.naedgecount(nw)
    if(nae==0) return(nw)
    
    na.el <- as.edgelist(is.na(nw))

    if(is.null(response)){
        d <- network.edgecount(nw,na.omit=TRUE)/network.dyadcount(nw,na.omit=TRUE)
        nimpute <- round(d*nae)
        nw[na.el] <- 0
        nw[na.el[sample.int(nae,nimpute),,drop=FALSE]] <- 1
    }else{
        x <- as.edgelist(nw,attrname=response)[,3]
        zeros <- network.dyadcount(nw,na.omit=TRUE)-length(x)
        nw[na.el] <- 0
        nw[na.el,names.eval=response,add.edges=TRUE] <- sample(c(0,x),nae,replace=TRUE,prob=c(zeros,rep(1,length(x))))
    }

    nw
}

# Given a vector, truncate all infinite (or, really, bigger in
# magnitude than replace=) values with replace= with the appropriate
# sign. Leave NAs and NANs alone.
.deinf <- function(x, replace=1/.Machine$double.eps) ifelse(is.nan(x) | abs(x)<replace, x, sign(x)*replace)



# executes expression, returns the result in a list with any warnings and errors
.catchToList <- function(expr) {
  val <- NULL
  myWarnings <- NULL
  wHandler <- function(w) {
    myWarnings <<- c(myWarnings, w$message)
    invokeRestart("muffleWarning")
  }
  myError <- NULL
  eHandler <- function(e) {
    myError <<- e$message
    NULL
  }
  val <- tryCatch(withCallingHandlers(expr, warning = wHandler), error = eHandler)
  list(value = val, warnings = myWarnings, error=myError)
} 

