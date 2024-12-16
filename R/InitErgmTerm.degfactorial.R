#' @templateVar name degm1lfactorial
#' @title Sum of log-factorials of degrees of non-isolates minus 1
#' @description This term adds one network statistic to the model,
#'   equalling \eqn{\sum_{i=1}^n \mathbb{I}(d_i>0) \log(d_i - 1)!}.
#'
#' @usage
#' # binary: degm1lfactorial
#'
#' @template ergmTerm-general
#'
#' @concept undirected
InitErgmTerm.degm1lfactorial <- function(nw, arglist, ..., version=packageVersion("ergm")) {
  a <- check.ErgmTerm(nw, arglist)
  list(name="sum_lfactorial_degreem1", coef.names="ldegm1lfactorial", minval = 0)
}


#' @templateVar name b2degm1lfactorial
#' @title Sum of log-factorials of degrees of Mode-2 non-isolates minus 1
#' @description This term adds one network statistic to the model,
#'   equalling \eqn{\sum_{i=b+1}^n \mathbb{I}(d_i>0) \log(d_i - 1)!}.
#'
#' @usage
#' # binary: b2degm1lfactorial
#'
#' @template ergmTerm-general
#'
#' @concept undirected
#' @concept bipartite
InitErgmTerm.b2degm1lfactorial <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, bipartite = TRUE)
  list(name="sum_lfactorial_idegreem1", coef.names="b2degm1lfactorial", minval = 0)
}


#' @templateVar name b1degm1lfactorial
#' @title Sum of log-factorials of degrees of Mode-1 non-isolates minus 1
#' @description This term adds one network statistic to the model,
#'   equalling \eqn{\sum_{i=1}^b \mathbb{I}(d_i>0) \log(d_i - 1)!}.
#'
#' @usage
#' # binary: b1degm1lfactorial
#'
#' @template ergmTerm-general
#'
#' @concept undirected
#' @concept bipartite
InitErgmTerm.b1degm1lfactorial <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, bipartite = TRUE)
  list(name="sum_lfactorial_odegreem1", coef.names="b1degm1lfactorial", minval = 0)
}


#' @templateVar name idegm1lfactorial
#' @title Sum of log-factorials of indegrees of non-in-isolates minus 1
#' @description This term adds one network statistic to the model,
#'   equalling \eqn{\sum_{i=1}^n \mathbb{I}(id_i>0) \log(id_i - 1)!}.
#'
#' @usage
#' # binary: idegm1lfactorial
#'
#' @template ergmTerm-general
#'
#' @concept directed
InitErgmTerm.idegm1lfactorial <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed = TRUE)
  list(name="sum_lfactorial_idegreem1", coef.names="idegm1lfactorial", minval = 0)
}


#' @templateVar name odegm1lfactorial
#' @title Sum of log-factorials of outdegrees of non-out-isolates minus 1
#' @description This term adds one network statistic to the model,
#'   equalling \eqn{\sum_{i=1}^n \mathbb{I}(od_i>0) \log(od_i - 1)!}.
#'
#' @usage
#' # binary: odegm1lfactorial
#'
#' @template ergmTerm-general
#'
#' @concept directed
InitErgmTerm.odegm1lfactorial <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed = TRUE)
  list(name="sum_lfactorial_odegreem1", coef.names="odegm1lfactorial", minval = 0)
}
