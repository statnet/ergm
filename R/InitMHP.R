#  File R/InitMHP.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
#===========================================================================
# The <InitMHP> file contains the following 24 functions for
# initializing the MHproposal object; each is prepended with 'InitMHP.'
#       <randomtoggle>      <CondOutDegreeDist> 
#       <TNT>               <ConstantEdges>    
#       <CondInDegree>      <CondDegree>
#       <CondOutDegree>     <HammingTNT>   
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
  MHproposal <- list(name = "randomtoggle", inputs=NULL)
  MHproposal
}

InitMHP.TNT <- function(arguments, nw) {
  MHproposal <- list(name = "TNT", inputs=NULL)
  MHproposal
}

InitMHP.CondDegree <- function(arguments, nw) {
  MHproposal <- list(name = "CondDegree", inputs=NULL)
  MHproposal
}
InitMHP.CondDegreeMix <- function(arguments, nw) {
  MHproposal <- list(name = "CondDegreeMix",
    inputs=get.vertex.attribute(nw,arguments$constraints$degreesmix$attrib))
  MHproposal
}

InitMHP.CondOutDegree <- function(arguments, nw) {
  MHproposal <- list(name = "CondOutDegree", inputs=NULL)
  if (!is.directed(nw)) # Really, this should never trigger, since the InitConstraint function should check.
    stop("The CondOutDegree proposal function does not work with an",
          "undirected network.")
  
  MHproposal
}

InitMHP.CondInDegree <- function(arguments, nw) {
  MHproposal <- list(name = "CondInDegree", inputs=NULL)
  if (!is.directed(nw)) # Really, this should never trigger, since the InitConstraint function should check.
    stop("The CondInDegree proposal function does not work with an",
          "undirected network.")
  MHproposal
}

InitMHP.CondB1Degree <- function(arguments, nw) {
  MHproposal <- list(name = "CondB1Degree", inputs=NULL)
  if (!is.bipartite(nw)) # Really, this should never trigger, since the InitConstraint function should check.
    stop("The CondB1Degree proposal function does not work with a non-bipartite network.")
  
  MHproposal
}

InitMHP.CondB2Degree <- function(arguments, nw) {
  MHproposal <- list(name = "CondB2Degree", inputs=NULL)
  if (!is.bipartite(nw)) # Really, this should never trigger, since the InitConstraint function should check.
    stop("The CondB2Degree proposal function does not work with a non-bipartite network.")
  MHproposal
}

InitMHP.CondDegreeDist <- function(arguments, nw) {
  MHproposal <- list(name = "CondDegreeDist", inputs=NULL)
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

InitMHP.CondInDegreeDist <- function(arguments, nw) {
  MHproposal <- list(name = "CondInDegreeDist", inputs=NULL)
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

InitMHP.CondOutDegreeDist <- function(arguments, nw) {
  MHproposal <- list(name = "CondOutDegreeDist", inputs=NULL)
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

InitMHP.ConstantEdges <- function(arguments, nw) {
  MHproposal <- list(name = "ConstantEdges", inputs=NULL)
  MHproposal
}

InitMHP.HammingConstantEdges <- function(arguments, nw) {
  MHproposal <- list(name = "HammingConstantEdges", inputs=NULL)
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteHammingConstantEdges"
  }
  MHproposal
}

InitMHP.HammingTNT <- function(arguments, nw) {
  MHproposal <- list(name = "HammingTNT", inputs=NULL)
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteHammingTNT"
  }
  MHproposal
}

InitMHP.randomtoggleNonObserved <- function(arguments, nw) {
  if(network.naedgecount(nw)==0){
   stop("The passed network does not have any non-observed dyads.\n Hence constraining to the observed will hold the network fixed at this network.\n Either the network or the constraint need to be altered.")
  }
  MHproposal <- list(name = "randomtoggleList", inputs=ergm.Cprepare.miss(nw))
  MHproposal
}

InitMHP.NonObservedTNT <- function(arguments, nw) {
  if(network.naedgecount(nw)==0){
   stop("The passed network does not have any non-observed dyads.\n Hence constraining to the observed will hold the network fixed at this network.\n Either the network or the constraint need to be altered.")
  }
  MHproposal <- list(name = "listTNT", inputs=ergm.Cprepare.miss(nw))
  MHproposal
}


InitMHP.fixedas <- function(arguments, nw){
	y0<-arguments$constraints$fixedas$free.dyads()
	## Given the list of toggleable dyads, no formation-specific proposal function is needed:
	MHproposal <- list(name = "randomtoggleList", inputs=c(ergm.Cprepare.el(y0)), pkgname="ergm")
	
	MHproposal
	
}

InitMHP.fixedasTNT <- function(arguments, nw){
	y0<-arguments$constraints$fixedas$free.dyads()
	## Given the list of toggleable dyads, no formation-specific proposal function is needed:
	MHproposal <- list(name = "listTNT", inputs=c(ergm.Cprepare.el(y0)), pkgname="ergm")
	
	MHproposal
	
}

InitMHP.fixallbut <- function(arguments, nw){
	y0<-arguments$constraints$fixallbut$free.dyads()
	## Given the list of toggleable dyads, no formation-specific proposal function is needed:
	MHproposal <- list(name = "randomtoggleList", inputs=c(ergm.Cprepare.el(y0)), pkgname="ergm")
	
	MHproposal
	
}


InitMHP.fixallbutTNT <- function(arguments, nw){
	y0<-arguments$constraints$fixallbut$free.dyads()
	## Given the list of toggleable dyads, no formation-specific proposal function is needed:
	MHproposal <- list(name = "listTNT", inputs=c(ergm.Cprepare.el(y0)), pkgname="ergm")
	
	MHproposal
	
}




