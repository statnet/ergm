#' ERGM-based conditional tie probabilities
#' 
#' @description 
#' Calculate model-predicted **conditional** tie probabilities for dyads given
#' the rest of the graph. Currently there are two methods implemented:
#' - Method for formula objects requires (1) an ERGM model formula with an existing
#' network object on the left hand side and terms on the right hand side, and
#' (2) and a vector of corresponding parameter values.
#' - Method for `ergm` objects, as returned by [ergm()], takes both the formula
#' and parameter values from the fitted model object.
#' 
#' @description 
#' Both methods can limit calculations to specific set of dyads of interest.
#'
#' @param object a formula or a fitted ERGM model object
#' @param theta numeric vector of ERGM model parameter values
#' @param output character, type of object returned. Defaults to `"data.frame"`.
#'   See section Value below.
#' @param ... other arguments passed to/from other methods. For [ergm.formula()]
#'   the arguments are passed to [ergmMPLE()]
#'
#' @return 
#' Type of object returned depends on the argument `output`. If
#' `output="data.frame"` the function will return a data frame with columns:
#' 
#' - `tail`, `head` -- indices of nodes identifying a dyad
#' - `p` -- predicted conditional tie probability
#' 
#' If `output="matrix"` the function will return an "adjacency matrix" with the
#' predicted conditional tie probabilities.
#' 
#' @method predict formula
#'
#' @export
#' @examples 
#' # A three-node empty directed network
#' net <- network.initialize(3, directed=TRUE)
#' 
#' # In homogeneous Bernoulli model with odds of a tie of 1/5 all ties are
#' # equally likely
#' predict(net ~ edges, log(1/5))
#' 
#' # Let's add a tie so that `net` has 1 tie out of possible 6 (so odds of 1/5)
#' net[1,2] <- 1
#' 
#' # Fit the model
#' fit <- ergm(net ~ edges)
#' 
#' # The p's should be identical
#' predict(fit)
#' 

predict.formula <- function(object, theta, output=c("data.frame", "matrix"), ...) {
  stopifnot(is.numeric(theta))
  output <- match.arg(output)
  predmat <- ergmMPLE(
    update(object, . ~ . + indices),
    output = "matrix",
    ...
  )$predictor
  stopifnot(length(theta) == (ncol(predmat)-2))
  p <- 1 / (1 + exp(predmat[,seq(1, length(theta)), drop=FALSE] %*% theta))
  switch(
    output,
    data.frame = data.frame(predmat[,c("tail", "head"), drop=FALSE], p),
    matrix = .NotYetImplemented()
  )
}

#' @rdname predict.formula
#' @method predict ergm
#' @export
predict.ergm <- function(object, ...) {
  predict.formula(
    object = object$formula,
    theta = ergm.eta(object$coef, object$etamap),
    ...
  )
}



if(FALSE) {
  logit <- function(p) log(p/(1-p))
  expit <- function(x) 1 / (1 + exp(-x))
  r1 <- predict_ergm(fit, flomarriage)
  
  data("flo")
  fit <- ergm(flomarriage ~ edges + gwesp(0.25, fixed = TRUE))
  # term: indicies
  z <- ergmMPLE( 
    update(fit$formula, . ~ . + indices), 
    output="matrix"
  )
  str(z)
  coef(fit)
  v <- ergm.eta(fit$coef, fit$etamap)
  z$predictor[,1:2] %*% v -> m
  
  p <- predict_formula(flomarriage ~ edges + gwesp(0.25, fixed = TRUE), c(-1, 0.5))
  
}


# Old brute force method. Slow but works. Kept for now for testing purposes.
brutal_predict_ergm <- function(object, net, ...) {
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
    form <- fit$formula
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
