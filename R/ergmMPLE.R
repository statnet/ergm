#  File R/ergmMPLE.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################

#' ERGM Predictors and response for logistic regression calculation of MPLE
#' 
#' Return the predictor matrix, response vector, and vector of weights that can
#' be used to calculate the MPLE for an ERGM.
#' 
#' The MPLE for an ERGM is calculated by first finding the matrix of change
#' statistics.  Each row of this matrix is associated with a particular pair
#' (ordered or unordered, depending on whether the network is directed or
#' undirected) of nodes, and the row equals the change in the vector of network
#' statistics (as defined in \code{formula}) when that pair is toggled from a 0
#' (no edge) to a 1 (edge), holding all the rest of the network fixed.  The
#' MPLE results if we perform a logistic regression in which the predictor
#' matrix is the matrix of change statistics and the response vector is the
#' observed network (i.e., each entry is either 0 or 1, depending on whether
#' the corresponding edge exists or not).
#' 
#' Using \code{output="matrix"}, note that the result of the fit may be
#' obtained from the \code{\link{glm}} function, as shown in the examples
#' below.
#' 
#' @param formula,constraints,obs.constraints An ERGM formula and
#'   (optional) constraint specification formulas. See \code{\link{ergm}}.
#' 
#' @param output Character, partially matched. See Value.
#'
#' @templateVar mycontrol control.ergm
#' @template control
#' @template verbose
#' @template expand.bipartite
#' @template basis
#'
#' @param \dots Additional arguments, to be passed to lower-level functions.
#' @return
#' 
#' If \code{output=="matrix"} (the default), then only the response, predictor,
#' and weights are returned; thus, the MPLE may be found by hand or the vector
#' of change statistics may be used in some other way. To save space, the
#' algorithm will automatically search for any duplicated rows in the predictor
#' matrix (and corresponding response values). \code{ergmMPLE} function will
#' return a list with three elements, \code{response}, \code{predictor}, and
#' \code{weights}, respectively the response vector, the predictor matrix, and
#' a vector of weights, which are really counts that tell how many times each
#' corresponding response, predictor pair is repeated.
#'
#' If \code{output=="dyadlist"}, as `"matrix"`, but rather than
#' coalescing the duplicated rows, every relation in the network that
#' is not fixed and is observed will have its own row in `predictor`
#' and element in `response` and `weights`, and `predictor` matrix
#' will have two additional rows at the start, `tail` and `head`,
#' indicating to which dyad the row and the corresponding elements
#' pertain.
#' 
#' If \code{output=="array"}, a list with similarly named three elements is
#' returned, but \code{response} is formatted into a sociomatrix;
#' \code{predictor} is a 3-dimensional array of with cell
#' \code{predictor[t,h,k]} containing the change score of term \code{k} for
#' dyad (\code{t},\code{h}); and \code{weights} is also formatted into a
#' sociomatrix, with an element being 1 if it is to be added into the
#' pseudolikelihood and 0 if it is not.
#' 
#' In particular, for a unipartite network, cells corresponding to self-loops,
#' i.e., \code{predictor[i,i,k]} will be \code{NA} and \code{weights[i,i]} will
#' be 0; and for a unipartite undirected network, lower triangle of each
#' \code{predictor[,,k]} matrix will be set to \code{NA}, with the lower
#' triangle of \code{weights} being set to 0.
#'
#' To all of the above output types, `attr(., "etamap")` is attached
#' containing the [mapping and offset information][ergm.eta].
#' 
#' If \code{output=="fit"}, then \code{ergmMPLE} simply calls the
#' \code{\link{ergm}} function with the \code{estimate="MPLE"} option set,
#' returning an object of class \code{ergm} that gives the fitted
#' pseudolikelihood model.
#' @seealso \code{\link{ergm}}, \code{\link{glm}}
#' @keywords regression models
#' @examples
#' 
#' data(faux.mesa.high)
#' formula <- faux.mesa.high ~ edges + nodematch("Sex") + nodefactor("Grade")
#' mplesetup <- ergmMPLE(formula)
#' 
#' # Obtain MPLE coefficients "by hand":
#' coef(glm(mplesetup$response ~ . - 1, data = data.frame(mplesetup$predictor),
#'          weights = mplesetup$weights, family="binomial"))
#' 
#' # Check that the coefficients agree with the output of the ergm function:
#' coef(ergmMPLE(formula, output="fit"))
#' 
#' # We can also format the predictor matrix into an array:
#' mplearray <- ergmMPLE(formula, output="array")
#' 
#' # The resulting matrices are big, so only print the first 8 actors:
#' mplearray$response[1:8,1:8]
#' mplearray$predictor[1:8,1:8,]
#' mplearray$weights[1:8,1:8]
#'
#' # Constraints are handled:
#' faux.mesa.high%v%"block" <- seq_len(network.size(faux.mesa.high)) %/% 4
#' mplearray <- ergmMPLE(faux.mesa.high~edges, constraints=~blockdiag("block"), output="array")
#' mplearray$response[1:8,1:8]
#' mplearray$predictor[1:8,1:8,]
#' mplearray$weights[1:8,1:8]
#'
#' # Or, a dyad list:
#' faux.mesa.high%v%"block" <- seq_len(network.size(faux.mesa.high)) %/% 4
#' mplearray <- ergmMPLE(faux.mesa.high~edges, constraints=~blockdiag("block"), output="dyadlist")
#' mplearray$response[1:8]
#' mplearray$predictor[1:8,]
#' mplearray$weights[1:8]
#' 
#' # Curved terms produce predictors on the canonical scale:
#' formula2 <- faux.mesa.high ~ gwesp
#' mplearray <- ergmMPLE(formula2, output="array")
#' # The resulting matrices are big, so only print the first 5 actors:
#' mplearray$response[1:5,1:5]
#' mplearray$predictor[1:5,1:5,1:3]
#' mplearray$weights[1:5,1:5]
#' @export ergmMPLE
ergmMPLE <- function(formula, constraints=~., obs.constraints=~-observed, output=c("matrix", "array", "dyadlist", "fit"), expand.bipartite=FALSE, control=control.ergm(),
                     verbose=FALSE, ..., basis=ergm.getnetwork(formula)){
  check.control.class("ergm", "ergmMPLE")
  handle.control.toplevel("ergm", ...)

  output <- match.arg(output)
  if (output=="fit") {
    return(
      ergm(formula, estimate="MPLE", control=control, verbose=verbose, constraints=constraints, obs.constraints=obs.constraints, basis=basis, ...)
    )
  }

  nw <- basis

  if(output %in% c("array", "dyadlist")) formula <- nonsimp_update.formula(formula, .~indices+.)

  # Construct the model
  model <- ergm_model(formula, nw, ..., term.options=control$term.options)

  # Handle the observation process constraints.
  tmp <- .handle.auto.constraints(nw, constraints, obs.constraints)
  nw <- tmp$nw; constraints <- tmp$constraints; constraints.obs <- tmp$constraints.obs
  
  if("constraints" %in% names(control$MCMC.prop.args)){
    conlist <- prune.ergm_conlist(control$MCMC.prop.args$constraints)
    class(conlist) <- "ergm_conlist"
  }else{
    conlist <- ergm_conlist(constraints, nw, term.options=control$term.options)
  }

  if("constraints" %in% names(control$obs.MCMC.prop.args)){
    conlist.obs <- prune.ergm_conlist(control$obs.MCMC.prop.args$constraints)
    class(conlist.obs) <- "ergm_conlist"
  }else{
    conlist.obs <- ergm_conlist(constraints.obs, nw, term.options=control$term.options)
  }

  fd <- as.rlebdm(conlist, conlist.obs, which="informative")

  # Get the MPLE predictors
  pl <- ergm.pl(nw, fd, model, verbose=verbose, control=control, ignore.offset=TRUE,...)

  structure(
    switch(output,
         matrix = list(response = pl$zy, predictor = pl$xmat.full,
           weights = pl$wend),
         dyadlist = {
           o <- order(pl$xmat.full[,"tail"], pl$xmat.full[,"head"])
           list(response = pl$zy[o], predictor = pl$xmat.full[o,,drop=FALSE],
                weights = pl$wend[o])
         },
         array = {
           # If expand.bipartite==TRUE, then no special treatment for bipartite networks is needed.
           bip <- if(!expand.bipartite) NVL(nw %n% "bipartite", 0) else 0

           vn <- if(all(is.na(nw %v% "vertex.names"))) 1:network.size(nw) else nw %v% "vertex.names"
           t.names <- if(bip) vn[seq_len(bip)] else vn
           h.names <- if(bip) vn[-seq_len(bip)] else vn
           term.names <- colnames(pl$xmat.full)[-(1:2),drop=FALSE]

           if(bip) pl$xmat.full[,2] <- pl$xmat.full[,2] - bip
           
           xa <- array(NA, dim = c(length(t.names), length(h.names), ncol(pl$xmat.full)-2), dimnames = list(tail = t.names, head = h.names, term = term.names))
           
           for(k in seq_along(term.names))
             xa[cbind(pl$xmat.full[,1:2,drop=FALSE],k)] <- pl$xmat.full[,k+2]

           ym <- replace(matrix(NA, length(t.names), length(h.names), dimnames = list(tail = t.names, head = h.names)),
                         pl$xmat.full[,1:2,drop=FALSE],
                         pl$zy)

           wm <- replace(matrix(0, length(t.names), length(h.names), dimnames = list(tail = t.names, head = h.names)),
                         pl$xmat.full[,1:2,drop=FALSE],
                         pl$wend)

           list(response = ym, predictor = xa, weights = wm)
         }
         ),
    etamap = model$etamap
  )
}

