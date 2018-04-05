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
# initializing the proposal object; each is prepended with 'InitMHP.'
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
# proposal list; when appropriate, proposal types are checked against
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
#   proposal: a list containing:
#        name   : the name of the proposal
#        inputs : a vector to be passed to the proposal
#        package: is "ergm"
#
############################################################################
InitMHP.randomtoggle <- function(arguments, nw) {
  proposal <- list(name = "randomtoggle", inputs=NULL)
  proposal
}

InitMHP.TNT <- function(arguments, nw) {
  proposal <- list(name = "TNT", inputs=NULL)
  proposal
}

InitMHP.CondDegree <- function(arguments, nw) {
  proposal <- list(name = "CondDegree", inputs=NULL)
  proposal
}
InitMHP.CondDegreeMix <- function(arguments, nw) {
  proposal <- list(name = "CondDegreeMix",
    inputs=get.vertex.attribute(nw,arguments$constraints$degreesmix$attrib))
  proposal
}

InitMHP.CondOutDegree <- function(arguments, nw) {
  proposal <- list(name = "CondOutDegree", inputs=NULL)
  if (!is.directed(nw)) # Really, this should never trigger, since the InitConstraint function should check.
    stop("The CondOutDegree proposal function does not work with an",
          "undirected network.")
  
  proposal
}

InitMHP.CondInDegree <- function(arguments, nw) {
  proposal <- list(name = "CondInDegree", inputs=NULL)
  if (!is.directed(nw)) # Really, this should never trigger, since the InitConstraint function should check.
    stop("The CondInDegree proposal function does not work with an",
          "undirected network.")
  proposal
}

InitMHP.CondB1Degree <- function(arguments, nw) {
  proposal <- list(name = "CondB1Degree", inputs=NULL)
  if (!is.bipartite(nw)) # Really, this should never trigger, since the InitConstraint function should check.
    stop("The CondB1Degree proposal function does not work with a non-bipartite network.")
  
  proposal
}

InitMHP.CondB2Degree <- function(arguments, nw) {
  proposal <- list(name = "CondB2Degree", inputs=NULL)
  if (!is.bipartite(nw)) # Really, this should never trigger, since the InitConstraint function should check.
    stop("The CondB2Degree proposal function does not work with a non-bipartite network.")
  proposal
}

InitMHP.CondDegreeDist <- function(arguments, nw) {
  proposal <- list(name = "CondDegreeDist", inputs=NULL)
  if (is.directed(nw)) {
    message("Warning:  Using the 'degreedist' constraint with a directed network ",
          "is currently perilous.  We recommend that you use 'outdegree' or ",
          "'idegrees' instead.")
  }
  if(is.bipartite(nw)){
     proposal$name <- "BipartiteCondDegreeDist"
  }
  proposal
}

InitMHP.CondInDegreeDist <- function(arguments, nw) {
  proposal <- list(name = "CondInDegreeDist", inputs=NULL)
  if (!is.directed(nw)) {
    message("Warning:  Using the 'idegreedist' constraint with an undirected network ",
          "is currently perilous.  We recommend that you use 'degreedist' ",
          " instead.")
  }
  if(is.bipartite(nw)){
     proposal$name <- "BipartiteCondDegreeDist"
  }
  proposal
}

InitMHP.CondOutDegreeDist <- function(arguments, nw) {
  proposal <- list(name = "CondOutDegreeDist", inputs=NULL)
  if (!is.directed(nw)) {
    message("Warning:  Using the 'odegreedist' constraint with an undirected network n",
          "is currently perilous.  We recommend that you use 'degreedist' ",
          " instead.")
  }
  if(is.bipartite(nw)){
     proposal$name <- "BipartiteCondDegreeDist"
  }
  proposal
}

InitMHP.ConstantEdges <- function(arguments, nw) {
  proposal <- list(name = "ConstantEdges", inputs=NULL)
  proposal
}

InitMHP.HammingConstantEdges <- function(arguments, nw) {
  proposal <- list(name = "HammingConstantEdges", inputs=NULL)
  if(is.bipartite(nw)){
    proposal$name <- "BipartiteHammingConstantEdges"
  }
  proposal
}

InitMHP.HammingTNT <- function(arguments, nw) {
  proposal <- list(name = "HammingTNT", inputs=NULL)
  if(is.bipartite(nw)){
    proposal$name <- "BipartiteHammingTNT"
  }
  proposal
}

InitMHP.randomtoggleNonObserved <- function(arguments, nw) {
  if(network.naedgecount(nw)==0){
   stop("The passed network does not have any non-observed dyads.\n Hence constraining to the observed will hold the network fixed at this network.\n Either the network or the constraint need to be altered.")
  }
  proposal <- list(name = "randomtoggleList", inputs=to_ergm_Cdouble(is.na(nw)))
  proposal
}

InitMHP.NonObservedTNT <- function(arguments, nw) {
  if(network.naedgecount(nw)==0){
   stop("The passed network does not have any non-observed dyads.\n Hence constraining to the observed will hold the network fixed at this network.\n Either the network or the constraint need to be altered.")
  }
  proposal <- list(name = "listTNT", inputs=to_ergm_Cdouble(is.na(nw)))
  proposal
}


InitMHP.fixedas <- function(arguments, nw){
	y0<-as.edgelist(arguments$constraints$fixedas$free_dyads, prototype=nw)
	## Given the list of toggleable dyads, no formation-specific proposal function is needed:
	proposal <- list(name = "randomtoggleList", inputs=to_ergm_Cdouble(y0), pkgname="ergm")
	
	proposal
	
}

InitMHP.fixedasTNT <- function(arguments, nw){
	y0<-as.edgelist(arguments$constraints$fixedas$free_dyads, prototype=nw)
	## Given the list of toggleable dyads, no formation-specific proposal function is needed:
	proposal <- list(name = "listTNT", inputs=to_ergm_Cdouble(y0), pkgname="ergm")
	
	proposal
	
}

InitMHP.fixallbut <- function(arguments, nw){
	y0<-as.edgelist(arguments$constraints$fixallbut$free_dyads, prototype=nw)
	## Given the list of toggleable dyads, no formation-specific proposal function is needed:
	proposal <- list(name = "randomtoggleList", inputs=to_ergm_Cdouble(y0), pkgname="ergm")
	
	proposal
	
}


InitMHP.fixallbutTNT <- function(arguments, nw){
	y0<-as.edgelist(arguments$constraints$fixallbut$free_dyads, prototype=nw)
	## Given the list of toggleable dyads, no formation-specific proposal function is needed:
	proposal <- list(name = "listTNT", inputs=to_ergm_Cdouble(y0), pkgname="ergm")
	
	proposal
	
}


InitMHP.RLE <- function(arguments, nw){
  proposal <- list(name = "RLE", inputs=to_ergm_Cdouble(as.rlebdm(arguments$constraints)), pkgname="ergm")
}

InitMHP.RLETNT <- function(arguments, nw){
  proposal <- list(name = "RLETNT", inputs=to_ergm_Cdouble(as.rlebdm(arguments$constraints)), pkgname="ergm")
}
