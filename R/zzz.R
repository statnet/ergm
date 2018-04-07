#  File R/zzz.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
.onAttach <- function(lib, pkg){
  #' @importFrom statnet.common statnetStartupMessage
  sm <- statnetStartupMessage("ergm", c("statnet","ergm.count","tergm"), TRUE)
  if(!is.null(sm)){
    packageStartupMessage(sm)
    packageStartupMessage(paste(c(strwrap(paste("NOTE: Versions before 3.6.1 had a bug in the implementation of the bd() constriant which distorted the sampled distribution somewhat. In addition, Sampson's Monks datasets had mislabeled vertices. See the NEWS and the documentation for more details.",sep="")),""),collapse="\n"))
  }
}

.onLoad <- function(lib, pkg){
  # . is used as a placeholder by stantet.common::NVL3().
  utils::globalVariables(".")

  .RegisterProposals()
  .RegisterConstraintImplications()
  .RegisterInitMethods()
}

.RegisterProposals <- function(){
  ergm_proposal_table("c", "Bernoulli", c("", "bd"),  0, "random", "randomtoggle")
  ergm_proposal_table("c", "Bernoulli", c("", "bd"),  1, "TNT", "TNT")
  ergm_proposal_table("c", "Bernoulli", c(".dyads",".dyads+bd"),  -2, "random", "RLE")
  ergm_proposal_table("c", "Bernoulli", c(".dyads",".dyads+bd"),  -1, "TNT", "RLETNT")
  ergm_proposal_table("c", "Bernoulli", "", -100, "TNT10", "TNT10")
  ergm_proposal_table("c", "Bernoulli", "degrees",  0, "random", "CondDegree")
  ergm_proposal_table("c", "Bernoulli", "degreesmix",  0, "random", "CondDegreeMix")
  ergm_proposal_table("c", "Bernoulli", c("idegrees+odegrees","b1degrees+b2degrees"),  0, "random", "CondDegree")
  ergm_proposal_table("c", "Bernoulli", "odegrees",  0, "random", "CondOutDegree")
  ergm_proposal_table("c", "Bernoulli", "idegrees",  0, "random", "CondInDegree")
  ergm_proposal_table("c", "Bernoulli", "b1degrees",  0, "random", "CondB1Degree")
  ergm_proposal_table("c", "Bernoulli", "b2degrees",  0, "random", "CondB2Degree")
  ergm_proposal_table("c", "Bernoulli", "degreedist",  0, "random", "CondDegreeDist")
  ergm_proposal_table("c", "Bernoulli", "idegreedist",  0, "random", "CondInDegreeDist")
  ergm_proposal_table("c", "Bernoulli", "odegreedist",  0, "random", "CondOutDegreeDist")
  ergm_proposal_table("c", "Bernoulli", c("bd+edges","edges"),  0, "random", "ConstantEdges")
  ergm_proposal_table("c", "Bernoulli", "edges+hamming",  0, "random", "HammingConstantEdges")
  ergm_proposal_table("c", "Bernoulli", "hamming",  0, "random", "HammingTNT")
  ergm_proposal_table("c", "Bernoulli", c("bd+observed","observed"),  0, "random", "randomtoggleNonObserved")
  ergm_proposal_table("c", "Bernoulli", c("bd+observed","observed"),  1, "TNT", "NonObservedTNT")
  ergm_proposal_table("c", "Bernoulli", c("blockdiag","bd+blockdiag"), 0, "random", "blockdiag")
  ergm_proposal_table("c", "Bernoulli", c("blockdiag","bd+blockdiag"), 1, "TNT", "blockdiagTNT")
  ergm_proposal_table("c", "Bernoulli", c("blockdiag+observed","bd+blockdiag+observed"),  0, "random", "blockdiagNonObserved")
  ergm_proposal_table("c", "Bernoulli", c("blockdiag+observed","bd+blockdiag+observed"),  1, "TNT", "blockdiagNonObservedTNT")
  ergm_proposal_table("c", "Bernoulli", "fixedas",  0, "random", "fixedas")
  ergm_proposal_table("c", "Bernoulli", "fixedas",  1, "TNT", "fixedasTNT")
  ergm_proposal_table("c", "Bernoulli", "fixallbut",  0, "random", "fixallbut")
  ergm_proposal_table("c", "Bernoulli", "fixallbut",  1, "TNT", "fixallbutTNT")

  
  ergm_proposal_table("c", "StdNormal", "",  0, "random", "StdNormal")
  ergm_proposal_table("c", "StdNormal", ".dyads",  0, "random", "DistRLE")

  ergm_proposal_table("c", "Unif", "",  0, "random", "Unif")
  ergm_proposal_table("c", "Unif", "observed",  0, "random", "UnifNonObserved")
  ergm_proposal_table("c", "Unif", ".dyads",  0, "random", "DistRLE")
  
  ergm_proposal_table("c", "DiscUnif", "",  0, "random", "DiscUnif")
  ergm_proposal_table("c", "DiscUnif", "observed",  0, "random", "DiscUnifNonObserved")  

  ergm_proposal_table("c", c("Unif","DiscUnif","StdNormal","Poisson","Binomial"), ".dyads",  -3, "random", "DistRLE")
}

.RegisterConstraintImplications <- function(){
  ergm.ConstraintImplications("edges", c())
  ergm.ConstraintImplications("degrees", c("edges", "idegrees", "odegrees", "idegreedist", "odegreedist", "degreedist", "bd"))
  ergm.ConstraintImplications("degreesmix", c("edges", "idegrees", "odegrees", "idegreedist", "odegreedist", "degreedist", "bd"))
  ergm.ConstraintImplications("odegrees", c("edges", "odegreedist"))
  ergm.ConstraintImplications("idegrees", c("edges", "idegreedist"))
  ergm.ConstraintImplications("b1degrees", c("edges"))
  ergm.ConstraintImplications("b2degrees", c("edges"))
  ergm.ConstraintImplications("degreedist", c("edges", "idegreedist", "odegreedist"))
  ergm.ConstraintImplications("idegreedist", c("edges"))
  ergm.ConstraintImplications("odegreedist", c("edges"))
  ergm.ConstraintImplications("bd", c())
  ergm.ConstraintImplications("hamming", c())
}

.RegisterInitMethods <- function(){
  ergm.init.methods("Bernoulli", c("MPLE", "CD", "zeros"))
  ergm.init.methods("StdNormal", c("CD","zeros"))
  ergm.init.methods("Unif", c("CD","zeros"))
  ergm.init.methods("DiscUnif", c("CD","zeros"))
}
