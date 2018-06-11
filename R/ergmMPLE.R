#  File R/ergmMPLE.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2018 Statnet Commons
#######################################################################
############################################################################
# The <ergmMPLE> function has different behavior based on whether the given
# formula should be fit or not. If so, <ergm> is called. If not, the elements
# needed for logit regression are computed and returned using <ergm.pl>
#
# --PARAMETERS--
#   formula : a formula as 'nw ~ term(s)'
#   fitmodel: whether to fit the model given by 'formula' (T or F);
#             default=FALSE
#   control : a list of parameters to control the fitting process; this
#             is ignored if 'fitmodel'=FALSE
#   verbose : whether the <ergm> or <ergm.pl> functions should be verbose
#             (T or F); default=FALSE
#   ..      : additional parameters that will be passed onto <ergm> or
#             <ergm.pl>
#
# --RETURNED--
#   if output
#     ="fit"  -- an ergm object, as returned by <ergm>
#     ="matrix" -- a list with 3 components:
#                response : the vector of dyad values; this is tabulated
#                           according to 'weights' 
#                predictor: the design matrix of change stats;  this is 
#                           tabulated according to 'weights' 
#                weights  : the weights for each entry/row of 'response'/
#                           'predictor'
#     ="array" -- a list with 3 components:
#                response : the sociomatrix
#                predictor: an array of change stats with dimensions
#                           as follows: (tail, head, term)
#                weights  : a sociomatrix with weights: typically 0 if a
#                           dyad was not visited, 1 if it was
#                Note that for undirected unipartite networks, the lower
#                triangle for each term's predictor matrix is set to NA
#                and its weight to 0.
#
###########################################################################



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
#' When \code{output="array"}, the \code{MPLE.max.dyad.types} control parameter
#' must be greater than \code{network.dyadcount(.)} of the response network, or
#' not all elements of the array that ought to be filled in will be.
#' 
#' @param formula An ERGM formula. See \code{\link{ergm}}.
#' @param fitmodel Deprecated. Use \code{output="fit"} instead.
#' @param output Character, partially matched. See Value.
#' @param as.initialfit Logical. Specifies whether terms are initialized with
#' argument \code{initialfit==TRUE} (the default). Generally, if \code{TRUE},
#' all curved ERGM terms will be treated as having their curved parameters
#' fixed. See Example.
#' @param control A list of control parameters for tuning the fitting of an
#' ERGM.  Most of these parameters are irrelevant in this context.  See
#' \code{\link{control.ergm}} for details about all of the control parameters.
#' @param verbose Logical; if \code{TRUE}, the program will print out some
#' additional information.
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
#' glm(mplesetup$response ~ . - 1, data = data.frame(mplesetup$predictor), 
#'     weights = mplesetup$weights, family="binomial")$coefficients
#' 
#' # Check that the coefficients agree with the output of the ergm function:
#' ergmMPLE(formula, output="fit")$coef
#' 
#' # We can also format the predictor matrix into an array:
#' mplearray <- ergmMPLE(formula, output="array")
#' 
#' # The resulting matrices are big, so only print the first 5 actors:
#' mplearray$response[1:5,1:5]
#' mplearray$predictor[1:5,1:5,]
#' mplearray$weights[1:5,1:5]
#' 
#' formula2 <- faux.mesa.high ~ gwesp(0.5,fix=FALSE)
#' 
#' # The term is treated as fixed: only the gwesp term is returned:
#' colnames(ergmMPLE(formula2, as.initialfit=TRUE)$predictor)
#' 
#' # The term is treated as curved: individual esp# terms are returned:
#' colnames(ergmMPLE(formula2, as.initialfit=FALSE)$predictor)
#' @export ergmMPLE
ergmMPLE <- function(formula, fitmodel=FALSE, output=c("matrix", "array", "fit"), as.initialfit = TRUE, control=control.ergm(),
                     verbose=FALSE, ...){
  if(!missing(fitmodel)){
      warning("Argument fitmodel= to ergmMPLE() has been deprecated and will be removed in a future version. Use output=\"fit\" instead.")
      if(fitmodel) output <- "fit"
  }
  check.control.class("ergm", "ergmMPLE")
  control.toplevel(...,myname="ergm")
  output <- match.arg(output)
  if (output=="fit") {
    return(ergm(formula, estimate="MPLE", control=control, verbose=verbose, ...))
  }

  if(output == "array") formula <- nonsimp_update.formula(formula, .~indices+.)
  
  nw <- ergm.getnetwork(formula)
  model <- ergm_model(formula, nw, initialfit=as.initialfit, term.options=control$term.options)
  fd <- ergm.design(nw, verbose=verbose)
  pl <- ergm.pl(nw, fd, model, verbose=verbose, control=control,...)

  switch(output,
         matrix = list(response = pl$zy, predictor = pl$xmat, 
           weights = pl$wend),
         array = {
           vn <- if(all(is.na(nw %v% "vertex.names"))) 1:network.size(nw) else nw %v% "vertex.names"
           t.names <- if(is.bipartite(nw)) vn[seq_len(nw %n% "bipartite")] else vn
           h.names <- if(is.bipartite(nw)) vn[-seq_len(nw %n% "bipartite")] else vn
           term.names <- colnames(pl$xmat)[-(1:2),drop=FALSE]
           
           xa <- array(NA, dim = c(length(t.names), length(h.names), ncol(pl$xmat)-2), dimnames = list(tail = t.names, head = h.names, term = term.names))
           
           for(k in seq_along(term.names))
             xa[cbind(pl$xmat[,1:2,drop=FALSE],k)] <- pl$xmat[,k+2]
           
           ym <- as.matrix(nw, matrix.type="adjacency")

           wm <- matrix(0, nrow(ym), ncol(ym))
           wm[cbind(pl$xmat[,1:2,drop=FALSE])] <- pl$wend

           list(response = ym, predictor = xa, weights = wm)
         }
         )
}

