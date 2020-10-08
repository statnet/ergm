
#' @param <%=MCMCType%>.prop Specifies the proposal (directly) and/or
#'   a series of "hints" about the structure of the model being
#'   sampled. The specification is in the form of a one-sided formula
#'   with hints separated by `+` operations. If the LHS exists and is
#'   a string, the proposal to be used is selected directly.
#'
#'   A common and default "hint" is `~TNT` (tie / no tie), indicating
#'   that the network is sparse and that the sample should put roughly
#'   equal weight on selecting a dyad with or without a tie as a
#'   candidate for toggling.
#'
#' @param  <%=MCMCType%>.prop.weights Specifies the proposal
#' distribution used in the <%=MCMCType%> Metropolis-Hastings algorithm.  Possible
#' choices depending on selected \code{reference} and \code{constraints}
#' arguments of the [ergm()] function, but often include \code{"TNT"}
#' and \code{"random"}, and the \code{"default"} is to use the one with the
#' highest priority available.
#' 
#' @param <%=MCMCType%>.prop.args An alternative, direct way of
#' specifying additional arguments to proposal.
