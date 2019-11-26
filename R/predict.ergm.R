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
#' @param type character element, one of `"response"` (default) or `"link"`
#'   whether the returned predictions are on the probability scale or on the
#'   scale of linear predictor. This is similar to [predict.glm()].
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

predict.formula <- function(object, theta,
                            type=c("response", "link"),
                            output=c("data.frame", "matrix"), ...) {
  stopifnot(is.numeric(theta))
  output <- match.arg(output)
  type <- match.arg(type)
  
  # Transform extended ergmMPLE() output to matrix with 0s on the diagonal
  .df_to_matrix <- function(d) {
    res <- tapply(predmat[,"p"], list(predmat[,"tail"], predmat[,"head"]), identity)
    diag(res) <- 0
    res
  }
  
  predmat <- ergmMPLE(
    update(object, . ~ . + indices),
    output = "matrix",
    ...
  )$predictor
  stopifnot(length(theta) == (ncol(predmat)-2))
  # Compute conditional Ps and cbind to ergmMPLE() output
  predmat <- cbind(predmat, p=drop(switch(
    type,
    link = predmat[,seq(1, length(theta)), drop=FALSE] %*% theta,
    response = 1 / (1 + exp( - predmat[,seq(1, length(theta)), drop=FALSE] %*% theta))
  ) ) )
  # Format output
  switch(
    output,
    data.frame = as.data.frame(predmat[,c("tail", "head", "p")]),
    matrix = .df_to_matrix(predmat)
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
