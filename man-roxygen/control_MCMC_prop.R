#' @param  <%=MCMCType%>.prop.weights Specifies the proposal
#' distribution used in the <%=MCMCType%> Metropolis-Hastings algorithm.  Possible
#' choices depending on selected \code{reference} and \code{constraints}
#' arguments of the [ergm()] function, but often include \code{"TNT"}
#' and \code{"random"}, and the \code{"default"} is to use the one with the
#' highest priority available.
#' 
#' The \code{TNT} (tie / no tie) option puts roughly equal weight on selecting
#' a dyad with or without a tie as a candidate for toggling, whereas the
#' \code{random} option puts equal weight on all possible dyads, though the
#' interpretation of \code{random} may change according to the constraints in
#' place.  When no constraints are in place, the default is TNT, which appears
#' to improve Markov chain mixing particularly for networks with a low edge
#' density, as is typical of many realistic social networks.
#' 
#' @param <%=MCMCType%>.prop.args An alternative, direct way of
#' specifying additional arguments to proposal.
