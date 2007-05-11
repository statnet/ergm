#  DH:  This file is still a work in progress, but for now it should
#  be set up so as not to break anything!

InitMHP.randomtoggle <- function(arguments, nw, model) {
  MHproposal <- list(name = "randomtoggle", args=NULL, package="statnet")
  if(is.bipartite(nw)){
    MHproposal$name <- "Bipartiterandomtoggle"
  }
  MHproposal
}

InitMHP.TNT <- function(arguments, nw, model) {
  MHproposal <- list(name = "TNT", args=NULL, package="statnet")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteTNT"
  }
  MHproposal
}

InitMHP.CondDegree <- function(arguments, nw, model) {
  MHproposal <- list(name = "CondDegree", args=NULL, package="statnet")
# if(is.bipartite(nw)){
#   MHproposal$name <- "BipartiteTNT"
# }
  MHproposal
}

InitMHP.ConstantEdges <- function(arguments, nw, model) {
  MHproposal <- list(name = "ConstantEdges", args=NULL, package="statnet")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteConstantEdges"
  }
# Check for redundant terms
  if("edges" %in% model$coef.names){
    stop("The model contains an 'edges' term and a proposal that ", 
         "holds the edges fixed. One of them is redundant. Please ", 
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

InitMHP.randomtoggleNonObserved <- function(arguments, nw, model) {
  MHproposal <- list(name = "randomtoggleNonObserved", args=NULL, package="statnet")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiterandomtoggleNonObserved"
  }
  MHproposal
}



