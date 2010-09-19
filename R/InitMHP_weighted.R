InitMHP.PseudoPoisson <- function(arguments, nw, model) {
  MHproposal <- list(name = "PseudoPoisson", args=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartitePseudoPoisson"
  }
  MHproposal
}
