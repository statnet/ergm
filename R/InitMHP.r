#  DH:  This file is still a work in progress, but for now it should
#  be set up so as not to break anything!

InitMHP.randomtoggle <- function(arguments, nw, model) {
  MHproposal <- list(name = "randomtoggle", args=NULL, package="statnet")
  if(is.bipartite(nw)){
    MHproposal$name <- "Bipartiterandomtoggle"
  }
  MHproposal
}

InitMHP.godfather <- function(arguments, nw, model) {
  # This is a proposal you can't refuse.  It is useful for checking on
  # the change statistics that result from a prescribed sequence of 
  # edge toggles.
  # It is assumed that "arguments" is a kx3 matrix,
  # where the first column is the times and the next two give a list
  # of edges.  This proposal orders the times, then performs one set of
  # toggles for each unique value of time.
  arguments <- arguments[order(arguments[,1]),]  
  MHproposal <- list(name = "nonrandom", args=arguments, package="statnet")  
  MHproposal
}

InitMHP.TNT <- function(arguments, nw, model) {
  MHproposal <- list(name = "TNT", args=NULL, package="statnet")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteTNT"
  }
  MHproposal
}

InitMHP.ConstantEdges <- function(arguments, nw, model) {
  MHproposal <- list(name = "ConstantEdges", args=NULL, package="statnet")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteConstantEdges"
  }
# Check for redundant terms
  if("edges" %in% model$coef.names){
    stop("The model contains an 'edges' term and a proposal that", 
         "holds the edges fixed. One of them is redundant. Please", 
         "restate the model.")
  }
  MHproposal
}

InitMHP.HammingConstantEdges <- function(arguments, nw, model) {
  MHproposal <- list(name = "HammingConstantEdges", args=NULL, package="statnet")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteHammingConstantEdges"
  }
  MHproposal
}

InitMHP.Hamming <- function(arguments, nw, model) {
  MHproposal <- list(name = "Hamming", args=NULL, package="statnet")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteHamming"
  }
  MHproposal
}

InitMHP.formation <- function(arguments, nw, model) {
  MHproposal <- list(name = "formation", args=NULL, package="statnet")
  if(is.bipartite(nw)){
    MHproposal$name <- "Bipartiteformation"
  }
  MHproposal
}

InitMHP.formationTNT <- function(arguments, nw, model) {
  MHproposal <- list(name = "FormationTNT", args=NULL, package="statnet")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteFormationTNT"
  }
  MHproposal
}



