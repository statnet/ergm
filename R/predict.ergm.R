#  File R/predict.ergm.R in package ergm, part of the Statnet suite of packages
#  for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
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
#' @param eta numeric vector of ERGM model *canonical* parameter values
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
#' @template basis
#' @param theta deprecated name for the parameter vector, renamed to
#'   `eta` since the canonical parameters are expected.
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
predict.formula <- function(object, eta,
                            conditional = TRUE,
                            type=c("response", "link"),
                            nsim = 100,
                            output = c("data.frame", "matrix"), ...,
                            basis = ergm.getnetwork(object), theta = NULL) {

  ## TODO: Remove the following after May 2026 and ergm 4.10.
  if (!is.null(theta)) {
    .Deprecate_once(
      msg = paste0(sQuote("predict.formula()"), "'s parameter argument is now ",
                   sQuote("eta"), ".")
    )
    eta <- theta
  }

  stopifnot(is.numeric(eta))
  eta <- deInf(eta)
  output <- match.arg(output)
  type <- match.arg(type)
  stopifnot(nsim >= 2)

  # Simulated unconditional Ps
  if(!conditional) {
    if(type != "response") 
      stop("type='link' for unconditional probabilities is not supported")
    predm <- predict_ergm_unconditional(object, eta, nsim, ..., basis = basis)
    return(
      switch(
        output,
        data.frame = with(arr_to_coo(predm, FALSE, na.rm = TRUE),
                          data.frame(coord, x)) |>
          setNames(c("tail", "head", "p")) |>
          subset(tail != head),
        matrix = predm
      )
    )
  }

  # Compute conditional Ps
  predmat <- ergmMPLE(object, output = "dyadlist", ..., basis = basis)$predictor
  stopifnot(length(eta) == (ncol(predmat) - 2))
  link <- drop(predmat[, -(1:2), drop = FALSE] %*% eta)
  pred <- data.frame(predmat[, 1:2, drop = FALSE],
                     p = switch(type, link = link, response = expit(link)))
  # Format output
  switch(
    output,
    data.frame = pred,
    matrix = {
      vnames <- basis %v% "vertex.names"
      arr_from_coo(pred$p, pred[1:2], dimnames = list(vnames, vnames)) |>
        set_diag(0)
    }
  )
}

predict_ergm_unconditional <- function(object, coef, nsim = 100, ...) {
  output <- function(s, ...) {
    vn <- s$nw0 %v% "vertex.names"
    arr_from_coo(TRUE, s$el, x0 = FALSE, dimnames = list(vn, vn))
  }
  simulate_formula(object = object, coef = coef, nsim = nsim,
                   output = output, ...) |>
    Reduce(`+`, x = _) / nsim
}


#' @rdname predict.formula
#' @export
predict.ergm <- function(object, ...) {
  if(is.valued(object)) stop("Prediction for valued ERGMs is not implemented at this time.")
  predict.formula(
    object = object$formula,
    eta = ergm.eta(coef(object), object$etamap),
    basis = NVL(object$newnetwork, object$network),
    ...
  )
}
