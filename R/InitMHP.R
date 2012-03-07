#  File ergm/R/InitMHP.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
########################################################################
# Each of the <InitMHP.X> functions initializes and returns a
# MHproposal list; when appropriate, proposal types are checked against
# covariates and network types for 1 of 2 side effects: to print warning
# messages or to halt execution (only <InitMHP.nobetweengroupties> can
# halt execution)
############################################################################


InitMHP.randomtoggle <- function(arguments, nw) {
  MHproposal <- list(name = "randomtoggle", inputs=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "Bipartiterandomtoggle"
  }
  MHproposal
}

InitMHP.TNT <- function(arguments, nw) {
  MHproposal <- list(name = "TNT", inputs=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteTNT"
  }
  MHproposal
}

InitMHP.TNT10 <- function(arguments, nw) {
  MHproposal <- list(name = "TNT10", inputs=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteTNT"
  }
  MHproposal
}

InitMHP.CondDegreeSimple <- function(arguments, nw) {
  MHproposal <- list(name = "CondDegreeSimple", inputs=NULL, package="ergm")
  MHproposal
}

InitMHP.CondDegree <- function(arguments, nw) {
  MHproposal <- list(name = "CondDegree", inputs=NULL, package="ergm")
  if (is.directed(nw)) {
    cat("Warning:  Using the 'degree' constraint with a directed network\n",
          "is currently perilous.  We recommend that you use 'outdegree' or\n",
          "'indegree' instead.\n")
  }
  if(is.bipartite(nw)){
    MHproposal$name <- "CondDegreeSimpleTetrad"
  }
  MHproposal
}


InitMHP.CondDegreeTetrad <- function(arguments, nw) {
  MHproposal <- list(name = "CondDegreeTetradToggles", inputs=NULL, package="ergm")
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

InitMHP.CondDegreeHexad <- function(arguments, nw) {
  MHproposal <- list(name = "CondDegreeHexadToggles", inputs=NULL, package="ergm")
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

InitMHP.CondDegreeDist <- function(arguments, nw) {
  MHproposal <- list(name = "CondDegreeDist", inputs=NULL, package="ergm")
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

InitMHP.CondInDegreeDist <- function(arguments, nw) {
  MHproposal <- list(name = "CondInDegreeDist", inputs=NULL, package="ergm")
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

InitMHP.CondOutDegreeDist <- function(arguments, nw) {
  MHproposal <- list(name = "CondOutDegreeDist", inputs=NULL, package="ergm")
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

#InitMHP.CondOutDegree <- function(arguments, nw) {
#  MHproposal <- list(name = "CondOutDegree", inputs=NULL, package="ergm")
#  if (!is.directed(nw)) {
#    cat("Warning:  The 'outdegree' constraint does not work with an\n",
#          "undirected network.  Switching to 'degree' constraint.\n")
#    return(InitMHP.CondDegree(arguments, nw))
#  }
#  MHproposal
#}

#InitMHP.CondInDegree <- function(arguments, nw) {
#  MHproposal <- list(name = "CondInDegree", inputs=NULL, package="ergm")
#  if (!is.directed(nw)) {
#    cat("Warning:  The 'indegree' constraint does not work with an\n",
#          "undirected network.  Switching to 'degree' constraint.\n")
#    return(InitMHP.CondDegree(arguments, nw))
#  }
#  MHproposal
#}

InitMHP.ConstantEdges <- function(arguments, nw) {
  MHproposal <- list(name = "ConstantEdges", inputs=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteConstantEdges"
  }

  MHproposal
}

InitMHP.HammingConstantEdges <- function(arguments, nw) {
  MHproposal <- list(name = "HammingConstantEdges", inputs=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteHammingConstantEdges"
  }
  MHproposal
}

InitMHP.HammingTNT <- function(arguments, nw) {
  MHproposal <- list(name = "HammingTNT", inputs=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteHammingTNT"
  }
  MHproposal
}

InitMHP.randomtoggleNonObserved <- function(arguments, nw) {
  if(network.naedgecount(nw)==0){
   stop("The passed network does not have any non-observed dyads.\n Hence constraining to the observed will hold the network fixed at this network.\n Either the network or the constraint need to be altered.")
  }
  MHproposal <- list(name = "randomtoggleNonObserved", inputs=ergm.Cprepare.miss(nw), package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiterandomtoggleNonObserved"
  }
  MHproposal
}


# This one does not have a C function.
InitMHP.nobetweengroupties <- function(arguments, nw) {
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
  inputs <- c(length(b), b, e)
  MHproposal <- list(name="nobetweengroupties", inputs = inputs, package="ergm")
  MHproposal
}



