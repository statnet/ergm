#======================================================================
# The <InitMHP> file contains the following 14 functions for
# initializing the MHproposal object; each is prepended with 'InitMHP.'
#       <randomtoggle>       <CondInDegreeDist>
#       <TNT>                <CondOutDegreeDist>
#       <TNT10>              <ConstantEdges>
#       <CondDegree>         <HammingConstantEdges>
#       <CondDegreeTetrad>   <HammingTNT>
#       <CondDegreeHexad>    <randomtoggleNonObserved>      
#       <CondDegreeDist>     <nobetweengroupties>
#======================================================================


########################################################################
# Each of the <InitMHP.X> functions initializes and returns a
# MHproposal list; when appropriate, proposal types are checked against
# covariates and network types for 1 of 2 side effects: to print warning
# messages or to halt execution (only <InitMHP.nobetweengroupties> can
# halt execution)
#
# --PARAMETERS--
#   arguments: is ignored by all but <InitMHP.nobetweengroupties>,
#              where 'arguments' is used to get the nodal attributes
#              via <get.node.attr>
#   nw       : the network given by the model
#   model    : the model for 'nw', as returned by <ergm.getmodel>
#
# --RETURNED--
#   MHproposal: a list containing:
#        name   : the name of the proposal
#        args   : NULL for all but <InitMHP.nobetweengroupties>,
#                 where 'args' is the vector of 3 concatenated components:
#                  1) the number of unique attribute values for the
#                     given attribute, say a1, a2, a3.. aM
#                  2) the table of attribute value frequencies for a1 - aM
#                  3) a vector of node id's; those having a1, then those
#                     having a2, ... then those having aM
#        package: is "ergm"
#
############################################################################



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

InitMHP.CondDegreeTetrad <- function(arguments, nw, model) {
  MHproposal <- list(name = "CondDegreeTetradToggles", args=NULL, package="ergm")
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

InitMHP.CondDegreeHexad <- function(arguments, nw, model) {
  MHproposal <- list(name = "CondDegreeHexadToggles", args=NULL, package="ergm")
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



