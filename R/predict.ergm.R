#' ERGM-based conditional tie probabilities
#'
#' This function calculates **conditional** tie probabilities for all dyads in a
#' given network based on an ERGM model.
#'
#' @param object fitted ERGM model object
#' @param net network object
#' @param ... other arguments passed to/from other methods
#'
#' @return 
#' N x N matrix of conditional tie probabilities, where N is the number of nodes
#' in `net`.
#' 
#' @method predict ergm
#'
#' @export

predict.ergm <- function(object, net, ...) {
  # Argument checking
  stopifnot(inherits(object, "ergm"))
  stopifnot(inherits(net, "network"))
  if( is.directed(object$network) != is.directed(net) )
    stop(
      "the model was fit to ",
      if(is.directed(object$network)) "directed" else "undirected",
      " network while `net` is ",
      if(is.directed(net)) "directed" else "undirected"
    )
  # Fix/enformulate curved ERGMs if necessary
  if(is.curved(object)) {
    fixed <- fix.curved(object)
    coeff <- fixed$theta
    form <- fixed$formula
  } else {
    coeff <- coef(object)
    form <- object$formula
  }
  n <-  network.size(net)
  # Substitute network object in formula's environment
  e <- rlang::env(net = net)
  rlang::f_env(form) <- e
  rlang::f_lhs(form) <- quote(net)
  # Calculate sufficient stats in `net`
  z <- summary(form)
  result <- matrix(0, n, n)
  # Loop over all dyads and calculate expits on toggle
  if(is.directed(net)) {
    iseq <- 1:n
  } else {
    iseq <- seq(1, n-1)
  }
  for(i in iseq) {
    if(is.directed(net)) {
      jseq <- setdiff(1:n, i)
    } else {
      jseq <- seq(i + 1, n)
    }
    for(j in jseq) {
      # Toggle the tie
      environment(form)$net[i,j] <- 1 - environment(form)$net[i,j]
      zAlt <- summary(form)
      # Changestats opposite for 1->0 toggles
      if(as.logical(environment(form)$net[i,j])) {
        d <- zAlt - z
      } else {
        d <- z - zAlt
      }
      # Expit
      result[i,j] <- (1.0 / (1.0 + exp(-sum(coeff * d))))
      # Symmetrize for undirected network
      if(!is.directed(net)) result[j,i] <- result[i,j]
      # Toggle back
      environment(form)$net[i,j] <- 1 - environment(form)$net[i,j]
    }
  }
  result
}
