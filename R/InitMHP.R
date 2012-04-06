#===========================================================================
# The <InitMHP> file contains the following 24 functions for
# initializing the MHproposal object; each is prepended with 'InitMHP.'
#       <randomtoggle>      <CondOutDegreeDist> 
#       <TNT>               <ConstantEdges>    
#       <TNT10>             <CondInDegree>      
#       <CondDegree>        <CondOutDegree>     <HammingTNT>   
#       <CondDegreeTetrad>         <HammingConstantEdges>
#       <CondDegreeHexad>            <randomtoggleNonObserved>
#       <CondDegreeDist>          <nobetweengroupties>
#       <CondInDegreeDist>  
#============================================================================


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
#
# --RETURNED--
#   MHproposal: a list containing:
#        name   : the name of the proposal
#        inputs : a vector to be passed to the proposal
#        package: is "ergm"
#
############################################################################
InitMHP.randomtoggle <- function(arguments, nw) {
  MHproposal <- list(name = "randomtoggle", inputs=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "Bipartiterandomtoggle"
  }
  MHproposal
}
#ergm.MHP.table("c", "Bernoulli", "",  0, "random", "randomtoggle")
#ergm.MHP.table("c", "Bernoulli", "bd",  0, "random", "randomtoggle")

InitMHP.TNT <- function(arguments, nw) {
  MHproposal <- list(name = "TNT", inputs=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteTNT"
  }
  MHproposal
}
#ergm.MHP.table("c", "Bernoulli", "",  1, "TNT", "TNT")
#ergm.MHP.table("c", "Bernoulli", "bd",  1, "TNT", "TNT")

InitMHP.TNT10 <- function(arguments, nw) {
  MHproposal <- list(name = "TNT10", inputs=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteTNT"
  }
  MHproposal
}
#ergm.MHP.table("c", "Bernoulli", "", -1, "TNT10", "TNT10")

InitMHP.CondDegree <- function(arguments, nw) {
  MHproposal <- list(name = "CondDegree", inputs=NULL, package="ergm")
  MHproposal
}
InitMHP.CondDegreeMix <- function(arguments, nw) {
  MHproposal <- list(name = "CondDegreeMix",
    inputs=get.vertex.attribute(nw,arguments$constraints$degreesmix$attrib),
    package="ergm")
  MHproposal
}
#ergm.MHP.table("c", "Bernoulli", "degrees",  0, "random", "CondDegree")
#ergm.MHP.table("c", "Bernoulli", "idegrees+odegrees",  0, "random", "CondDegree")
#ergm.MHP.table("c", "Bernoulli", "b1degrees+b2degrees",  0, "random", "CondDegree")

InitMHP.CondOutDegree <- function(arguments, nw) {
  MHproposal <- list(name = "CondOutDegree", inputs=NULL, package="ergm")
  if (!is.directed(nw)) # Really, this should never trigger, since the InitConstraint function should check.
    stop("The CondOutDegree proposal function does not work with an",
          "undirected network.")
  
  MHproposal
}
#ergm.MHP.table("c", "Bernoulli", "odegrees",  0, "random", "CondOutDegree")

InitMHP.CondInDegree <- function(arguments, nw) {
  MHproposal <- list(name = "CondInDegree", inputs=NULL, package="ergm")
  if (!is.directed(nw)) # Really, this should never trigger, since the InitConstraint function should check.
    stop("The CondInDegree proposal function does not work with an",
          "undirected network.")
  MHproposal
}
#ergm.MHP.table("c", "Bernoulli", "idegrees",  0, "random", "CondInDegree")

InitMHP.CondB1Degree <- function(arguments, nw) {
  MHproposal <- list(name = "CondB1Degree", inputs=NULL, package="ergm")
  if (!is.bipartite(nw)) # Really, this should never trigger, since the InitConstraint function should check.
    stop("The CondB1Degree proposal function does not work with a non-bipartite network.")
  
  MHproposal
}
#ergm.MHP.table("c", "Bernoulli", "b1degrees",  0, "random", "CondB1Degree")

InitMHP.CondB2Degree <- function(arguments, nw) {
  MHproposal <- list(name = "CondB2Degree", inputs=NULL, package="ergm")
  if (!is.bipartite(nw)) # Really, this should never trigger, since the InitConstraint function should check.
    stop("The CondB2Degree proposal function does not work with a non-bipartite network.")
  MHproposal
}
#ergm.MHP.table("c", "Bernoulli", "b2degrees",  0, "random", "CondB2Degree")

InitMHP.CondDegreeDist <- function(arguments, nw) {
  MHproposal <- list(name = "CondDegreeDist", inputs=NULL, package="ergm")
  if (is.directed(nw)) {
    cat("Warning:  Using the 'degreedist' constraint with a directed network\n",
          "is currently perilous.  We recommend that you use 'outdegree' or\n",
          "'idegrees' instead.\n")
  }
  if(is.bipartite(nw)){
     MHproposal$name <- "BipartiteCondDegreeDist"
  }
  MHproposal
}
#ergm.MHP.table("c", "Bernoulli", "degreedist",  0, "random", "CondDegreeDist")

InitMHP.CondInDegreeDist <- function(arguments, nw) {
  MHproposal <- list(name = "CondInDegreeDist", inputs=NULL, package="ergm")
  if (!is.directed(nw)) {
    cat("Warning:  Using the 'idegreedist' constraint with an undirected network\n",
          "is currently perilous.  We recommend that you use 'degreedist'\n",
          " instead.\n")
  }
  if(is.bipartite(nw)){
     MHproposal$name <- "BipartiteCondDegreeDist"
  }
  MHproposal
}
#ergm.MHP.table("c", "Bernoulli", "idegreedist",  0, "random", "CondInDegreeDist")

InitMHP.CondOutDegreeDist <- function(arguments, nw) {
  MHproposal <- list(name = "CondOutDegreeDist", inputs=NULL, package="ergm")
  if (!is.directed(nw)) {
    cat("Warning:  Using the 'odegreedist' constraint with an undirected network\n",
          "is currently perilous.  We recommend that you use 'degreedist'\n",
          " instead.\n")
  }
  if(is.bipartite(nw)){
     MHproposal$name <- "BipartiteCondDegreeDist"
  }
  MHproposal
}
#ergm.MHP.table("c", "Bernoulli", "odegreedist",  0, "random", "CondOutDegreeDist")

InitMHP.ConstantEdges <- function(arguments, nw) {
  MHproposal <- list(name = "ConstantEdges", inputs=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteConstantEdges"
  }
  MHproposal
}
#ergm.MHP.table("c", "Bernoulli", "bd+edges",  0, "random", "ConstantEdges")
#ergm.MHP.table("c", "Bernoulli", "edges",  0, "random", "ConstantEdges")

InitMHP.HammingConstantEdges <- function(arguments, nw) {
  MHproposal <- list(name = "HammingConstantEdges", inputs=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteHammingConstantEdges"
  }
  MHproposal
}
#ergm.MHP.table("c", "Bernoulli", "edges+hamming",  0, "random", "HammingConstantEdges")

InitMHP.HammingTNT <- function(arguments, nw) {
  MHproposal <- list(name = "HammingTNT", inputs=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteHammingTNT"
  }
  MHproposal
}
#ergm.MHP.table("c", "Bernoulli", "hamming",  0, "random", "HammingTNT")

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
#ergm.MHP.table("c", "Bernoulli", "bd+observed",  0, "random", "randomtoggleNonObserved")
#ergm.MHP.table("c", "Bernoulli", "observed",  0, "random", "randomtoggleNonObserved")

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



