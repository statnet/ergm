#' @param by,homophily the optional argument `by` specifies a vertex attribute (see Specifying
#'   Vertex attributes and Levels (`?nodal_attributes`) for details).
#'   If this is specified and `homophily` is `TRUE` ,
#'   then degrees are calculated using the subnetwork consisting of only
#'   edges whose endpoints have the same value of the `by` attribute.
#'   If `by` is specified and
#'   `homophily` is `FALSE` (the default), then separate degree range
#'   statistics are calculated for nodes having each separate
#'   value of the attribute.
