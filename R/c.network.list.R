#' @rdname network.list
#'
#' @param check_attr logical; should the attributes of the combined network
#'   lists be checked for consistency. If `TRUE` inconsistencies result in
#'   errors.
#'
#' @importFrom purrr map
#' @export

c.network.list <- function(..., check_attr = TRUE) {
  dots <- list(...)

  # Merge network lists without attributes
  lapply(dots, function(x) {
    attributes(x) <- NULL
    x
  }) -> l_networks
  rval <- do.call("c", l_networks)

  # Check attributes
  if(check_attr) {
    # Names of attributes to check with `all.equal()`
    attr_names <- c("coefficients", "control", "response",
                    "formula", "constraints", "reference")
    for(an in attr_names) {
      al <- map(dots, ~ attr(.x, an))
      ok <- all_identical(al, all.equal)
      if(!ok) stop(paste0("network lists do not have equal values on attribute ", an))
    }

    # Check if "stats" have identical columns
    l_stats <- map(dots, ~ attr(.x, "stats"))
    ok <- all_identical(
      lapply(l_stats, function(x) colnames(x)),
      fun = identical
    )
    if(!ok) stop("network lists do not have identical columns of 'stats' attribute")
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
