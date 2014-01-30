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

ergmMPLE <- function(formula, fitmodel=FALSE, output=c("matrix", "array", "fit"), as.initialfit = TRUE, control=control.ergm(),
                     verbose=FALSE, ...){
  if(!missing(fitmodel)){
      warning("Argument fitmodel= to ergmMPLE() has been deprecated and will be removed in a future version. Use output=\"fit\" instead.")
      if(fitmodel) output <- "fit"
  }
  check.control.class("ergm")
  output <- match.arg(output)
  if (output=="fit") {
    return(ergm(formula, estimate="MPLE", control=control, verbose=verbose, ...))
  }

  if(output == "array") formula <- ergm.update.formula(formula, .~indices+.)
  
  nw <- ergm.getnetwork(formula)
  model <- ergm.getmodel(formula, nw, initialfit=as.initialfit)
  Clist <- ergm.Cprepare(nw, model)
  Clist.miss <- ergm.design(nw, model, verbose=verbose)
  pl <- ergm.pl(Clist, Clist.miss, model, verbose=verbose, control=control,...)

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

