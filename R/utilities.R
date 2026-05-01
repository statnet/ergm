#  File R/ergm.utility.R in package ergm, part of the Statnet suite of packages
#  for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################

#' @rdname ergm
#' @importFrom methods is
#' @export
is.ergm <- function(object)
{
    is(object,"ergm")
}


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
   nb1 <- b1.size(g)
   nb2 <- b2.size(g)
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
#' vertex-level attributes from one [`network`] object to another,
#' ignoring some standard properties by default.
#' 
#' 
#' @param to the [`network`] that attributes should be copied to
#' @param from the [`network`] that attributes should be copied to
#' @param ignore vector of charcter names of network attributes that should not
#' be copied. Default is the standard list of network properties created by
#' [network.initialize()]
#' @return returns the \code{to} network, with attributes copied from
#' \code{from}
#' @note does not check that networks are of the same size, etc
#' @seealso [set.vertex.attribute()],
#' [set.network.attribute()]
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

  if(!is.dyad.independent(constraints) || !is.dyad.independent(constraints.obs)){
    message("Model and/or observational constraints are not dyad-independent. Dyad imputation cannot be used. Please ensure your LHS network satisfies all constraints.")
    if(output=="network") return(nw)
    else return(ergm_state(nw))
  }

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

## A is a matrix. x is a column vector that may contain Infs, NaNs,
## and NAs. This function computes A %*% x, treating 0*Inf, 0*NA, and
## 0*NaN as 0.
mat_by_coef <- function(A, x) {
  .Call("mat_by_coef", A, x)
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

is.const.sample <- function(x){
  if(is.mcmc.list(x)) x <- as.matrix(x)
  isTRUE(all.equal(apply(x,2,stats::sd), rep(0,ncol(x)), check.names=FALSE))
}

#' Obtain the ABI version of a library package (e.g., \pkg{ergm}) against which a client package was compiled.
#'
#' @param client name of the \R package that is `LinkingTo` the library package; defaults to \pkg{ergm} itself.
#' @param lib name of the \R package that is being linked to; defaults to \pkg{ergm} itself.
#'
#' @note This function is [memoise()]d.
#'
#' @return an integer containing the ABI version, typically `major*1e6 + minor`, or `NA` if none was stored.
#'
#' @noRd
packageBuiltABI <- memoise(function(client = "ergm", lib = "ergm"){
  NVL(.Call("GetBuiltABIVersion_wrapper", as.character(client), as.character(lib)),
      NA_integer_)
})

message_once <- once(message)

#' ABI version checking
#'
#' Check if package `client` was built under the same ABI version of
#' `lib` as the currently loaded version of `lib` and perform the
#' specified action if not.
#'
#' @param action What to do in case of a mismatch: \describe{
#'
#' \item{`"stop"`, `"abort"`}{stop with an error}
#'
#' \item{`"warning"`}{warn and proceed}
#'
#' \item{`"message"`, `"inform"`}{print a message and proceed}
#'
#' \item{`"silent"`}{return the value without side-effects}
#'
#' \item{`"disable"`}{skip the check, always returning `TRUE`}
#' }
#' Partial matching is supported.
#'
#' @return `TRUE` if the ABI matches or is not exported by `lib`,
#'   `FALSE` if exported by `lib` but not `client` or does not match.
#'
#' @noRd
check_ABI <- once(function(client = "ergm", lib  = "ergm", action = getOption("ergm.ABI.action")){
  action <- match.arg(action, c("stop", "abort", "warning", "message", "inform", "silent", "disable"))
  if(action == "disable") return(TRUE)

  if(client == lib) return(TRUE)

  if(NVL3(packageDescription(client)$LinkingTo, grepl(paste0("([^a-zA-Z0-9.]|^)", lib, "([^a-zA-Z0-9.]|$)"), .), FALSE)){

    asVer <- function(x) if(is.na(x)) NA_integer_ else package_version(paste(c(x %/% 1e6, x %% 1e6), collapse = "."))
    clientABI <- asVer(packageBuiltABI(client, lib))
    libABI <- asVer(packageBuiltABI(lib, lib))

    if(is.na(libABI)) return(TRUE)

    if(action %in% c("inform", "message")) {
      msg <- if(is.na(clientABI)) sprintf("NOTE: Possible C application binary interface (ABI) version mismatch between %s and %s detected. If you experience mysterious crashes and errors, rebuilding %s may help. You can silence this message by setting %s. (See %s for possible values.)",
                                          sQuote(client), sQuote(lib), sQuote(client), sQuote("options(ergm.ABI.action = ...)"), sQuote("options?ergm"))
             else if(clientABI != libABI) sprintf("NOTE: Possible C application binary interface (ABI) version mismatch in %s detected (compiled against %s %s, running against %s). If you experience mysterious crashes and errors, rebuilding %s may help. You can silence this message by setting %s. (See %s for possible values.)",
                                                  sQuote(client), sQuote(lib), clientABI, libABI, sQuote(client), sQuote("options(ergm.ABI.action = ...)"), sQuote("options?ergm"))
      msg <- NVL3(msg, paste(strwrap(., indent = 0, exdent = 2), collapse = "\n"))
    } else {
      msg <- if(is.na(clientABI)) sprintf("Package %s was compiled with an older version of %s that did not save its C application binary interface (ABI) version information at the time. Inconsistent ABI versions may cause malfunctions ranging from incorrect results to R crashing. Please rebuild the package against the current %s version. If you believe message to be spurious, you can override by changing %s (See %s for possible values.) and proceed at your own risk.",
                                          sQuote(client), sQuote(lib), sQuote(lib), sQuote("options(ergm.ABI.action = ...)"), sQuote("options?ergm"))
             else if(clientABI != libABI) sprintf("Package %s was compiled with %s that had one C application binary interface (ABI), but it is being loaded with %s that has a different C ABI version. (Respectively, %s and %s; note that the ABI version is not the same as the package version and may lag behind it.) Inconsistent versions may result in malfunctions ranging from incorrect results to R crashing. Please rebuild the package with the current %s version. If you believe message to be spurious, you can override by changing %s (See %s for possible values.) and proceed at your own risk.",
                                                  sQuote(client), sQuote(lib), sQuote(lib), clientABI, libABI, sQuote(lib), sQuote("options(ergm.ABI.action = ...)"), sQuote("options?ergm"))
      msg <- NVL3(msg, paste(c("", strwrap(., indent = 2, exdent = 2)), collapse = "\n"))
    }


    if(!is.null(msg)){
      switch(action,
             abort =,
             stop = stop(msg, call. = FALSE),
             warning = warning_once(msg, call. = FALSE, immediate. = TRUE),
             inform =,
             message = message_once(msg))
      FALSE
    }else TRUE
  }else TRUE
})

# Useful shorthands
logit <- qlogis
expit <- plogis

# TODO: Document and move to statnet.common.

`%[f]%` <- function(x, y) if (is.function(y)) x[y(x)] else x[y]

`%[f]%<-` <- function(x, y, value) replace(x, y, value)

`%[!f]%` <- function(x, y) if (is.function(y)) x[!y(x)] else x[!y]

`%[!f]%<-` <- function(x, y, value) {
  replace(x,
          if (is.function(y)) !y(x)
          else if (is.logical(y)) !y
          else -y,
          value)
}

# Evaluate RHS substituting LHS for ., then use the result to index
# LHS.
`%[.]%` <- function(x, y) {
  x[eval.parent(do.call(substitute, list(substitute(y), list(. = x))))]
}

# Note: This tries to guess whether to subset the RHS by comparing
# lengths. This may not be a good heuristic in all cases, so should
# not be relied on.
`%[.]%<-` <- function(x, y, value) {
  i <- eval.parent(do.call(substitute, list(substitute(y), list(. = x))))
  if (length(value) > length(seq_along(x)[i])) value <- value[i]
  x[i] <- value
  x
}

`%[.,]%<-` <- function(x, y, value) {
  i <- eval.parent(do.call(substitute, list(substitute(y), list(. = x))))
  if (length(value) > length(seq_len(dim(x)[1L])[i])) value <- value[i]
  x[i, ] <- value
  x
}

`%[,.]%<-` <- function(x, y, value) {
  i <- eval.parent(do.call(substitute, list(substitute(y), list(. = x))))
  if (length(value) > length(seq_len(dim(x)[2L])[i])) value <- value[i]
  x[, i] <- value
  x
}


`%[.|.]%<-` <- function(x, y, value) {
  i <- eval.parent(do.call(substitute, list(substitute(y), list(. = x))))
  if (length(value) > length(seq_len(dim(x)[2L])[i])) value <- value[i]
  x[, i] <- x[i, ] <- value
  x
}

`%[.,.]%<-` <- function(x, y, value) {
  i <- eval.parent(do.call(substitute, list(substitute(y), list(. = x))))
  if (!is.null(dim(value)) && dim(value)[2L] > length(seq_len(dim(x)[2L])[i]))
    value <- value[i, i]
  x[i, i] <- value
  x
}

rowcolnames <- function(x, ...) {
  rownames <- rownames(x, ...)
  colnames <- colnames(x, ...)
  stopifnot(identical(rownames, colnames))

  rownames
}

`rowcolnames<-` <- function(x, value) {
  colnames(x) <- rownames(x) <- value
  x
}

## This should be in 'network'.

b1.size <- function(x) if(is.bipartite(x)) x %n% "bipartite" else FALSE
b2.size <- function(x) if(is.bipartite(x)) network.size(x) - b1.size(x) else FALSE

## For handling skipping statistics.

is.NA <- function(x) is.na(x) & !is.nan(x)

## Enclose variable names with nonalphanumeric characters in them in
## backticks.
quote_var_name <- function(x) ifelse(grepl("[^a-zA-Z0-9._]", x), paste0("`", x, "`"), x)

##

strwrppst <- function(x, ..., parsep = "") {
  paste0(
    sapply(lapply(x, strwrap, ...),
           paste0, "\n", collapse = ""),
    collapse = parsep
  )
}
