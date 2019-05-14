#  File R/ergm.utility.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
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
#' that name][degreedist-constraint], try
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

#' @describeIn ergm-deprecated Use the [`pending_update_network`] "API".
#' @export newnw.extract
newnw.extract<-function(oldnw,z=NULL,output="network",response=NULL){
  .Deprecate_once('pending_network_update "API"')
  if(is(oldnw,"pending_update_network") && is.null(z)){
    class(oldnw) <- "network"
    z <- oldnw%n%".update"
    delete.network.attribute(oldnw, ".update")
  }

  newedgelist <- .extract_z_edgelist(z, response)
  
  newnw<-network.update(oldnw,newedgelist,matrix.type="edgelist",output=output)
  if(!is.null(response)){
    newnwweights <- newedgelist[,3]
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

#' @rdname ergm-deprecated
#' @export standardize.network
standardize.network <- function(nw, preserve.eattr=TRUE){
  .Deprecate_once(msg=paste0(sQuote("standardize.network"), " has been obviated by improvements to ", sQuote("network"), "."))
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

single.impute.dyads <- function(nw, response=NULL, constraints=NULL, constraints.obs=NULL, min_informative=NULL, default_density=NULL, output=c("network","pending_update_network"), verbose=FALSE){
  output <- match.arg(output)
  stopifnot(!is.null(constraints)||is.null(constraints.obs))
  
  if(!is.null(constraints)){
    imputable <- as.rlebdm(constraints, constraints.obs, "missing")
    nae <- NVL3(imputable, sum(.), 0)
    if(nae) na.el <- as.edgelist(imputable) # FIXME: Avoid creating edgelists.
  }else{
    nae <- network.naedgecount(nw)
    if(nae) na.el <- as.edgelist(is.na(nw))
  }
  if(nae==0) return(nw)

  if(verbose) message("Imputing ", nae, " dyads is required.")

  el2s <- function(el) apply(el, 1, paste, collapse=",")
  s2el <- function(s) matrix(as.integer(do.call(rbind,strsplit(s,","))),ncol=2)

  min_informative <- NVL3(min_informative, if(is.function(.)) .(nw) else ., 0)
  default_density <- if(is.function(default_density)) default_density(nw)

  if(!is.null(constraints)){ # Constraints
    informative <- as.rlebdm(constraints, constraints.obs, "informative")
    nonzeros <- as.rlebdm(nw)
    if(is.null(response)){
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
      x <- as.edgelist(nw,attrname=response)
      x.el <- x[,1:2,drop=FALSE]
      x <- x.el[! el2s(x.el)%in%el2s(na.el),3]
      zeros <- sum(informative) - length(x)
    }
  }else{ # No Constraints
    if(is.null(response)){
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
      x <- as.edgelist(nw,attrname=response)[,3]
      zeros <- network.dyadcount(nw,na.omit=TRUE)-length(x)
    }
  }
  
  if(is.null(response)){
    if(verbose) message("Imputing ", nimpute, " edges at random.")
    i.new <- sample.int(nae,nimpute)
    if(output=="network"){
      y.cur <- nw[na.el]
      i.na <- which(is.na(y.cur))
      i.cur <- which(y.cur!=0)
      todel <- union(setdiff(i.cur, i.new), setdiff(i.na, i.new))
      toadd <- union(setdiff(i.new, i.cur), intersect(i.na, i.new))
      nw[na.el[c(todel,toadd),,drop=FALSE]] <- rep(0:1, c(length(todel),length(toadd)))
    }else{ # pending_update_network
      el <- s2el(union(setdiff(el2s(as.edgelist(nw)), el2s(na.el)), el2s(na.el[i.new,,drop=FALSE])))
      nw <- pending_update_network(nw, list(newedgelist = el))
    }
  }else{
    if(output=="network"){
      nw[na.el,names.eval=response,add.edges=TRUE] <- sample(c(0,x),nae,replace=TRUE,prob=c(zeros,rep(1,length(x))))
    }else{ # pending_update_network
      el <- as.edgelist(nw, attrname=response)
      el <- el[!el2s(el[,-3,drop=FALSE])%in%el2s(na.el),,drop=FALSE]
      el <- rbind(el, cbind(na.el, sample(c(0,x),nae,replace=TRUE,prob=c(zeros,rep(1,length(x))))))
      el <- el[el[,3]!=0,,drop=FALSE]
      nw <- pending_update_network(nw, list(newedgelist = el), response=response)
    }
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


# TODO: Move to statnet.common and export after next release.
#' Evaluate an expression, restarting on error
#'
#' A pair of functions paralleling [eval()] and [evalq()] that make
#' multiple attempts at evaluating an expression, retrying on error up
#' to a specified number of attempts, and optionally evaluating
#' another expression before restarting.
#'
#' @param expr an expression to be retried; note the difference
#'   between [eval()] and [evalq()].
#' @param retries number of retries to make; defaults to
#'   `"eval.retries"` option, or 5.
#' @param before_retry if given, an expression that will be evaluated
#'   before each retry if the initial attempt fails; it is evaluated
#'   in the same environment and with the same quoting semantics as
#'   `expr`, but its errors are not handled.
#' @param envir,enclos see [eval()].
#' @param verbose Whether to output retries.
#'
#' @details Note if `expr` returns a `"try-error"` object (returned by
#'   [try()]), it will be treated as an error. This behavior may
#'   change in the future.
#'
#' @return Results of evaluating `expr`, including side-effects such
#'   as variable assignments, if successful in `retries` retries.
#'
#' @examples
#' x <- 0
#' persistEvalQ({if((x<-x+1)<3) stop("x < 3") else x}, before_retry = {cat("Will try incrementing...\n")})
#'
#' x <- 0
#' e <- quote(if((x<-x+1)<3) stop("x < 3") else x)
#' persistEval(e, before_retry = quote(cat("Will try incrementing...\n")))
#' @noRd
.persistEval <- function(expr, retries=NVL(getOption("eval.retries"), 5), before_retry,
                        envir = parent.frame(),
                        enclos = if (is.list(envir) ||
                                     is.pairlist(envir)) parent.frame() else baseenv(), verbose=FALSE){
  for(attempt in seq_len(retries)){
    out <- try(eval(expr, envir=envir, enclos=enclos), silent=TRUE)
    if(!is(out, "try-error")) return(out)
    else{
      if(!missing(before_retry)) eval(before_retry, envir=envir, enclos=enclos)
      if(verbose) inform(paste0("Retrying: retry ", attempt,"."))
    }
  }
  out <- eval(expr, envir=envir, enclos=enclos)
}

## #' @rdname persistEval
#' @noRd
.persistEvalQ <- function(expr, retries=NVL(getOption("eval.retries"), 5), before_retry,
                         envir = parent.frame(),
                         enclos = if (is.list(envir) ||
                                      is.pairlist(envir)) parent.frame() else baseenv(), verbose=FALSE){
  expr <- substitute(expr)
  before_retry <- substitute(before_retry)
  envir <- force(envir)
  enclos <- force(enclos)

  .persistEval(expr=expr, retries=retries, before_retry=before_retry, envir=envir, enclos=enclos, verbose=verbose)
}
