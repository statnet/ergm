#  File R/network.list.R in package ergm, part of the Statnet suite of packages
#  for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
#' A convenience container for a list of [`network`] objects, output
#' by [simulate.ergm()] among others.
#'
#' @param object,x a `list` of networks or a `network.list` object.
#' @param ... for `network.list`, additional attributes to be set on
#'   the network list; for others, arguments passed down to
#'   lower-level functions.
#'
#' @note Functions from the [simulate.ergm()] family can also return
#'   lists of lists of networks. In this case, they have an additional
#'   class `"network.list.list"`. At this time, it only affects
#'   printing.
#'
#' @aliases network.list.list
#' @export network.list
network.list <- function(object,...){
  if(any(!sapply(object, is.network))) stop("network.list() takes a list of networks as its first argument.")
  ns <- ...names()
  for(i in seq_len(...length())){
    attr(object, ns[i]) <- ...elt(i)
  }
  class(object) <- c("network.list","list")
  object
}

#' @describeIn network.list A [print()] method for network lists.
#' @export
print.network.list <- function(x, stats.print=FALSE, ...) {
  summary.network.list(x, stats.print=stats.print, ...)
  invisible(x)
}

#' @describeIn network.list A [summary()] method for network lists.
#' 
#' @param stats.print Logical: If TRUE, print network statistics.
#' @param net.print Logical: If TRUE, print network overviews.
#' @param net.summary Logical: If TRUE, print network summaries.
#' @seealso [simulate.ergm()]
#' @examples
#' 
#' # Draw from a Bernoulli model with 16 nodes
#' # and tie probability 0.1
#' #
#' g.use <- network(16, density=0.1, directed=FALSE)
#' #
#' # Starting from this network let's draw 3 realizations
#' # of a model with edges and 2-star terms
#' #
#' g.sim <- simulate(~edges+kstar(2), nsim=3, coef=c(-1.8, 0.03),
#'                basis=g.use, control=control.simulate(
#'                  MCMC.burnin=100000,
#'                  MCMC.interval=1000))
#' print(g.sim)
#' summary(g.sim)
#' 
#' @export
summary.network.list <- function (object, stats.print=TRUE, 
                       net.print=FALSE, net.summary=FALSE, ...){

  if (inherits(object, "network.list.list"))
    cat("List of lists of ", length(object), "*", length(object[[1]]),
        " Networks\n")
  else cat("List of ", length(object), " Networks\n")
  attrmap<-list(formula="Model: ",
                reference="Reference: ",
                constraints="Constraints: ",
                coef="Parameters:\n",
                stats="Stored network statistics:\n")
  if(!stats.print) attrmap$stats <- NULL
  for(a in names(attrmap)){
    if(a %in% names(attributes(object))){
      s<-attrmap[[a]]
      cat(s)
      if(substr(s,nchar(s),nchar(s))=="\n"){
        print(attr(object,a))
        cat("\n")
      }else{
        cat(format(attr(object,a)),"\n")
      }
    }
  }
  if(net.print){
    cat("Network overviews:\n")
    o <- object
    attributes(o)<-list()
    print(o)
  }
  if(net.summary){
    cat("Network summaries:\n")
    print(lapply(object,summary,...))
  }
  object
}





#' @rdname network.list
#'
#' @param check_attr Logical: should the attributes of the combined network
#'   lists be checked for consistency. If `TRUE` inconsistencies result in
#'   errors.
#'
#' @importFrom purrr map
#' @export
#'
#' @examples
#' # Simulate some more
#' g.sim2 <- simulate(~edges+kstar(2), nsim=3, coef=c(-1.8, 0.03),
#'   basis=g.use, control=control.simulate(
#'   MCMC.burnin=100000,
#'   MCMC.interval=1000))
#'
#' # Merge the simulations
#' g.simall <- c(g.sim, g.sim2)
#' length(g.simall) # 6
#'

c.network.list <- function(..., check_attr = TRUE) {
  dots <- list(...)

  # Merge network lists without attributes
  rval <- do.call("c", map(dots, unclass))

  # Check attributes
  if (check_attr) {
    # Names of attributes to check with `all.equal()`
    attr_names <- c("coefficients", "control", "response",
                    "formula", "constraints", "reference")
    for (an in attr_names) {
      al <- map(dots, ~ attr(.x, an))
      ok <- all_identical(al, all.equal)
      if (!ok) stop("network lists do not have equal values on attribute ", an)
    }

    # Check if "stats" have identical columns
    l_stats <- map(dots, ~ attr(.x, "stats"))
    ok <- all_identical(
      map(l_stats, function(x) colnames(x)),
      .p = identical
    )
    if (!ok) stop("network lists do not have identical columns of ",
                  sQuote("stats"), " attribute")
  }
  # Return the list of networks with attributes merged or taken from the
  # first object
  structure(
    rval,
    class = "network.list",
    coefficients = attr(dots[[1]], "coefficients"),
    control = attr(dots[[1]], "control"),
    response = attr(dots[[1]], "response"),
    stats = structure(
      do.call("rbind", map(dots, ~attr(.x, "stats"))),
      monitored = do.call("c", map(dots, ~ attr(attr(.x, "stats"), "monitored")))
    ),
    formula = attr(dots[[1]], "formula"),
    constraints = attr(dots[[1]], "constraints"),
    reference = attr(dots[[1]], "reference")
  )
}
