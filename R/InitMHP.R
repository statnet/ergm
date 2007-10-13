#  DH:  This file is still a work in progress, but for now it should
#  be set up so as not to break anything!

InitMHP.randomtoggle <- function(arguments, nw, model) {
  MHproposal <- list(name = "randomtoggle", args=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "Bipartiterandomtoggle"
  }
  MHproposal
}

InitMHP.TNT <- function(arguments, nw, model) {
  MHproposal <- list(name = "TNT", args=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteTNT"
  }
  MHproposal
}

InitMHP.CondDegree <- function(arguments, nw, model) {
  MHproposal <- list(name = "CondDegree", args=NULL, package="ergm")
  if (is.directed(nw)) {
    print("Warning:  Using the 'degree' constraint with a directed network",
          "is currently perilous.  We recommend that you use 'outdegree' or",
          "'indegree' instead.")
  }
  MHproposal
}

InitMHP.CondOutDegree <- function(arguments, nw, model) {
  MHproposal <- list(name = "CondOutDegree", args=NULL, package="ergm")
  if (!is.directed(nw)) {
    print("Warning:  The 'outdegree' constraint does not work with an",
          "undirected network.  Switching to 'degree' constraint.")
    return(InitMHP.CondDegree(arguments, nw, model))
  }
  MHproposal
}

InitMHP.CondInDegree <- function(arguments, nw, model) {
  MHproposal <- list(name = "CondInDegree", args=NULL, package="ergm")
  if (!is.directed(nw)) {
    print("Warning:  The 'indegree' constraint does not work with an",
          "undirected network.  Switching to 'degree' constraint.")
    return(InitMHP.CondDegree(arguments, nw, model))
  }
  MHproposal
}

InitMHP.ConstantEdges <- function(arguments, nw, model) {
  MHproposal <- list(name = "ConstantEdges", args=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteConstantEdges"
  }
# Check for redundant terms
  if("edges" %in% model$coef.names){
    cat("Warning: The model contains an 'edges' term and a proposal that\n", 
         "holds the edges fixed. The 'edges' term will be ignored.\n")
  }
  MHproposal
}

InitMHP.HammingConstantEdges <- function(arguments, nw, model) {
  MHproposal <- list(name = "HammingConstantEdges", args=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteHammingConstantEdges"
  }
  MHproposal
}

InitMHP.Hamming <- function(arguments, nw, model) {
  MHproposal <- list(name = "Hamming", args=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteHamming"
  }
  MHproposal
}

InitMHP.formation <- function(arguments, nw, model) {
  MHproposal <- list(name = "Formation", args=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteFormation"
  }
  MHproposal
}

InitMHP.formationTNT <- function(arguments, nw, model) {
  MHproposal <- list(name = "FormationTNT", args=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteFormationTNT"
  }
  MHproposal
}

InitMHP.dissolution <- function(arguments, nw, model) {
  MHproposal <- list(name = "Dissolution", args=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteDissolution"
  }
  MHproposal
}

InitMHP.randomtoggleNonObserved <- function(arguments, nw, model) {
  MHproposal <- list(name = "randomtoggleNonObserved", args=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiterandomtoggleNonObserved"
  }
  MHproposal
}

InitMHP.nobetweengroupties <- function(arguments, nw, model) {
  x <- get.node.attr(nw, arguments, "InitMHP.nobetweengroupties")
  if(any(is.na(x)) || any(table(x)==1)) {
    stop("nobetweengroups may not be used with a nodal covariate containing ",
         "NAs or nonrepeated values")
  }
  a <- sort(x)
  b <- table(a)
  d <- unique(a)
  e <- unlist(sapply(d, grep, x))
  f <- b*(b-1)
  args <- c(length(b), b, e)
  MHproposal <- list(name="nobetweengroupties", args = args, package="ergm")
  MHproposal
}



