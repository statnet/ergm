#  File R/predict.ergm.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################
#' ERGM-based tie probabilities
#' 
#' @description 
#' Calculate model-predicted **conditional** and **unconditional** tie
#' probabilities for dyads in the given network. Conditional probabilities of a
#' dyad given the state of all the remaining dyads in the graph are computed
#' exactly. Unconditional probabilities are computed through simulating networks
#' using the given model. Currently there are two methods implemented:
#' - Method for formula objects requires (1) an ERGM model formula with an existing
#' network object on the left hand side and model terms on the right hand side, and
#' (2) a vector of corresponding parameter values.
#' - Method for `ergm` objects, as returned by [ergm()], takes both the formula
#' and parameter values from the fitted model object.
#' 
#' @description 
#' Both methods can limit calculations to specific set of dyads of interest.
#'
#' @param object a formula or a fitted ERGM model object
#' @param theta numeric vector of ERGM model parameter values
#' @param conditional logical whether to compute conditional or unconditional
#'   predicted probabilities
#' @param nsim integer, number of simulated networks used for computing
#'   unconditional probabilities. Defaults to 100.
#' @param type character element, one of `"response"` (default) or `"link"` -
#'   whether the returned predictions are on the probability scale or on the
#'   scale of linear predictor. This is similar to `type` argument of [predict.glm()].
#' @param output character, type of object returned. Defaults to `"data.frame"`.
#'   See section Value below.
#' @param ... other arguments passed to/from other methods. For the `predict.formula` method, if
#'   `conditional=TRUE` arguments are passed to [ergmMPLE()]. If `conditional=FALSE` arguments
#'   are passed to [simulate_formula()].
#'
#' @return 
#' Type of object returned depends on the argument `output`. If
#' `output="data.frame"` the function will return a data frame with columns:
#' 
#' - `tail`, `head` -- indices of nodes identifying a dyad
#' - `p` -- predicted conditional tie probability
#' 
#' If `output="matrix"` the function will return an "adjacency matrix" with the
#' predicted probabilities. Diagonal values are 0s.
#' 
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
#' @export
predict.formula <- function(object, theta,
                            conditional = TRUE,
                            type=c("response", "link"),
                            nsim = 100,
                            output=c("data.frame", "matrix"), ...) {
  stopifnot(is.numeric(theta))
  theta <- statnet.common::deInf(theta)
  output <- match.arg(output)
  type <- match.arg(type)
  stopifnot(nsim >= 2)
  
  # Transform extended ergmMPLE() output to matrix with 0s on the diagonal
  .df_to_matrix <- function(d) {
    N <- max(predmat[,c("tail", "head")])
    res <- replace(matrix(NA, N, N), as.matrix(d[,c("tail", "head")]), d[,"p"])
    diag(res) <- 0
    res
  }

  # Matrix to data.frame
  .matrix_to_df <- function(m, name=".value") {
    unames <- sort(unique(unlist(dimnames(m))))
    d <- as.data.frame(as.table(m), stringsAsFactors=FALSE)
    names(d) <- c("tail", "head", name)
    tail <- d$tail <- match(d$tail, unames)
    head <- d$head <- match(d$head, unames)
    d[tail!=head, , drop=FALSE]
  }

  # Simulated unconditional Ps
  if(!conditional) {
    if(type != "response") 
      stop("type='link' for unconditional probabilities is not supported")
    predm <- predict_ergm_unconditional(object=object, coef=theta, nsim=nsim, ...)
    return(
      switch(
        output,
        data.frame = .matrix_to_df(predm, name="p"),
        matrix = predm
      )
    )
  }
  
  predmat <- ergmMPLE(
    statnet.common::nonsimp_update.formula(object, . ~ indices + . ),
    output = "matrix", # reduced to number of informative dyads in ergm.pl
    ...
  )$predictor
  stopifnot(length(theta) == (ncol(predmat)-2))
  # Compute conditional Ps and cbind to ergmMPLE() output
  predmat <- cbind(predmat, p=drop(switch(
    type,
    link = predmat[,-(1:2), drop=FALSE] %*% theta,
    response = 1 / (1 + exp( - predmat[,-(1:2), drop=FALSE] %*% theta))
  ) ) )
  # Format output
  switch(
    output,
    data.frame = as.data.frame(predmat[,c("tail", "head", "p")]),
    matrix = {
      # Get vertex names
      vnames <- ergm.getnetwork(object) %v% "vertex.names"
      structure(.df_to_matrix(predmat), dimnames = list(vnames, vnames))
    }
  )
}

predict_ergm_unconditional <- function(object, coef, nsim=100, output="network", ...) {
  netlist <- simulate_formula(object=object, coef=coef, nsim=nsim, output=output, ...)
  mats <- vapply(
    netlist, as.matrix, 
    matrix(0, ncol=network.size(netlist[[1]]), nrow=network.size(netlist[[1]])),
    matrix.type="adjacency"
  )
  apply(mats, 1:2, mean)
}


#' @rdname predict.formula
#' @export
predict.ergm <- function(object, ...) {
  if(is.valued(object)) stop("Prediction for valued ERGMs is not implemented at this time.")
  predict.formula(
    object = object$formula,
    theta = ergm.eta(coef(object), object$etamap),
    ...
  )
}
