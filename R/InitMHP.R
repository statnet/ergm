#  File ergm/R/InitMHP.R
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
# Copyright 2003 Mark S. Handcock, University of Washington
#                David R. Hunter, Penn State University
#                Carter T. Butts, University of California - Irvine
#                Steven M. Goodreau, University of Washington
#                Martina Morris, University of Washington
# Copyright 2007 The statnet Development Team
######################################################################
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

InitMHP.TNT10 <- function(arguments, nw, model) {
  MHproposal <- list(name = "TNT10", args=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteTNT"
  }
  MHproposal
}

InitMHP.CondDegree <- function(arguments, nw, model) {
  MHproposal <- list(name = "CondDegree", args=NULL, package="ergm")
  if (is.directed(nw)) {
    cat("Warning:  Using the 'degree' constraint with a directed network\n",
          "is currently perilous.  We recommend that you use 'outdegree' or\n",
          "'indegree' instead.\n")
  }
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteCondDegHexadToggles"
  }
  MHproposal
}

InitMHP.CondDegreeDist <- function(arguments, nw, model) {
  MHproposal <- list(name = "CondDegreeDist", args=NULL, package="ergm")
  if (is.directed(nw)) {
    cat("Warning:  Using the 'degreedist' constraint with a directed network\n",
          "is currently perilous.  We recommend that you use 'outdegree' or\n",
          "'indegree' instead.\n")
  }
  if(is.bipartite(nw)){
     MHproposal$name <- "BipartiteCondDegreeDist"
  }
  MHproposal
}

InitMHP.CondInDegreeDist <- function(arguments, nw, model) {
  MHproposal <- list(name = "CondInDegreeDist", args=NULL, package="ergm")
  if (!is.directed(nw)) {
    cat("Warning:  Using the 'indegreedist' constraint with an undirected network\n",
          "is currently perilous.  We recommend that you use 'degreedist'\n",
          " instead.\n")
  }
  if(is.bipartite(nw)){
     MHproposal$name <- "BipartiteCondDegreeDist"
  }
  MHproposal
}

InitMHP.CondOutDegreeDist <- function(arguments, nw, model) {
  MHproposal <- list(name = "CondOutDegreeDist", args=NULL, package="ergm")
  if (!is.directed(nw)) {
    cat("Warning:  Using the 'outdegreedist' constraint with an undirected network\n",
          "is currently perilous.  We recommend that you use 'degreedist'\n",
          " instead.\n")
  }
  if(is.bipartite(nw)){
     MHproposal$name <- "BipartiteCondDegreeDist"
  }
  MHproposal
}

#InitMHP.CondOutDegree <- function(arguments, nw, model) {
#  MHproposal <- list(name = "CondOutDegree", args=NULL, package="ergm")
#  if (!is.directed(nw)) {
#    cat("Warning:  The 'outdegree' constraint does not work with an\n",
#          "undirected network.  Switching to 'degree' constraint.\n")
#    return(InitMHP.CondDegree(arguments, nw, model))
#  }
#  MHproposal
#}

#InitMHP.CondInDegree <- function(arguments, nw, model) {
#  MHproposal <- list(name = "CondInDegree", args=NULL, package="ergm")
#  if (!is.directed(nw)) {
#    cat("Warning:  The 'indegree' constraint does not work with an\n",
#          "undirected network.  Switching to 'degree' constraint.\n")
#    return(InitMHP.CondDegree(arguments, nw, model))
#  }
#  MHproposal
#}

InitMHP.ConstantEdges <- function(arguments, nw, model) {
  MHproposal <- list(name = "ConstantEdges", args=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteConstantEdges"
  }
# Check for redundant terms
  e <- c("edges", "meandeg", "density")
  if(any(m<-(e %in% model$coef.names))){    
    cat(paste("Warning: The model contains the", e[m], 
              "term and the proposal constraints\nhold", e[m],
              "constant.  This term will be ignored.\n"))
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

InitMHP.HammingTNT <- function(arguments, nw, model) {
  MHproposal <- list(name = "HammingTNT", args=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteHammingTNT"
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



