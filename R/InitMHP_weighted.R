InitMHP.Poisson <- function(arguments, nw, model) {
  MHproposal <- list(name = "Poisson", args=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartitePoisson"
  }
  MHproposal
}


InitMHP.DescRank <- function(arguments, nw, model) {
  MHproposal <- list(name = "CompleteOrdering", args=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "CompleteOrderingBipartite"
  }
  MHproposal
}
