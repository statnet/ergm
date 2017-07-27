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
  sm <- statnetStartupMessage("ergm", c("statnet","ergm.count","tergm"), TRUE)
  if(!is.null(sm)){
    packageStartupMessage(sm)
    packageStartupMessage(paste(c(strwrap(paste("NOTE: Versions before 3.6.1 had a bug in the implementation of the bd() constriant which distorted the sampled distribution somewhat. In addition, Sampson's Monks datasets had mislabeled vertices. See the NEWS and the documentation for more details.",sep="")),""),collapse="\n"))
  }
}

.onLoad <- function(lib, pkg){
  .RegisterMHPs()
  .RegisterConstraintImplications()
  .RegisterInitMethods()
}

.RegisterMHPs <- function(){
  ergm.MHP.table("c", "Bernoulli", "",  0, "random", "randomtoggle")
  ergm.MHP.table("c", "Bernoulli", "bd",  0, "random", "randomtoggle")
  ergm.MHP.table("c", "Bernoulli", "",  1, "TNT", "TNT")
  ergm.MHP.table("c", "Bernoulli", "bd",  1, "TNT", "TNT")
  ergm.MHP.table("c", "Bernoulli", "", -1, "TNT10", "TNT10")
  ergm.MHP.table("c", "Bernoulli", "degrees",  0, "random", "CondDegree")
  ergm.MHP.table("c", "Bernoulli", "degreesmix",  0, "random", "CondDegreeMix")
  ergm.MHP.table("c", "Bernoulli", "idegrees+odegrees",  0, "random", "CondDegree")
  ergm.MHP.table("c", "Bernoulli", "b1degrees+b2degrees",  0, "random", "CondDegree")
  ergm.MHP.table("c", "Bernoulli", "odegrees",  0, "random", "CondOutDegree")
  ergm.MHP.table("c", "Bernoulli", "idegrees",  0, "random", "CondInDegree")
  ergm.MHP.table("c", "Bernoulli", "b1degrees",  0, "random", "CondB1Degree")
  ergm.MHP.table("c", "Bernoulli", "b2degrees",  0, "random", "CondB2Degree")
  ergm.MHP.table("c", "Bernoulli", "degreedist",  0, "random", "CondDegreeDist")
  ergm.MHP.table("c", "Bernoulli", "idegreedist",  0, "random", "CondInDegreeDist")
  ergm.MHP.table("c", "Bernoulli", "odegreedist",  0, "random", "CondOutDegreeDist")
  ergm.MHP.table("c", "Bernoulli", "bd+edges",  0, "random", "ConstantEdges")
  ergm.MHP.table("c", "Bernoulli", "edges",  0, "random", "ConstantEdges")
  ergm.MHP.table("c", "Bernoulli", "edges+hamming",  0, "random", "HammingConstantEdges")
  ergm.MHP.table("c", "Bernoulli", "hamming",  0, "random", "HammingTNT")
  ergm.MHP.table("c", "Bernoulli", "bd+observed",  0, "random", "randomtoggleNonObserved")
  ergm.MHP.table("c", "Bernoulli", "bd+observed",  1, "TNT", "NonObservedTNT")
  ergm.MHP.table("c", "Bernoulli", "observed",  0, "random", "randomtoggleNonObserved")
  ergm.MHP.table("c", "Bernoulli", "observed",  1, "TNT", "NonObservedTNT")
  ergm.MHP.table("c", "Bernoulli", "blockdiag", 0, "random", "blockdiag")
  ergm.MHP.table("c", "Bernoulli", "blockdiag", 1, "TNT", "blockdiagTNT")
  ergm.MHP.table("c", "Bernoulli", "bd+blockdiag", 0, "random", "blockdiag")
  ergm.MHP.table("c", "Bernoulli", "bd+blockdiag", 1, "TNT", "blockdiagTNT")
  ergm.MHP.table("c", "Bernoulli", "blockdiag+observed",  0, "random", "blockdiagNonObserved")
  ergm.MHP.table("c", "Bernoulli", "blockdiag+observed",  1, "TNT", "blockdiagNonObservedTNT")
  ergm.MHP.table("c", "Bernoulli", "bd+blockdiag+observed",  0, "random", "blockdiagNonObserved")
  ergm.MHP.table("c", "Bernoulli", "bd+blockdiag+observed",  1, "TNT", "blockdiagNonObservedTNT")
  ergm.MHP.table("c", "Bernoulli", "fixedas",  0, "random", "fixedas")
  ergm.MHP.table("c", "Bernoulli", "fixedas",  1, "TNT", "fixedasTNT")
  ergm.MHP.table("c", "Bernoulli", "fixallbut",  0, "random", "fixallbut")
  ergm.MHP.table("c", "Bernoulli", "fixallbut",  1, "TNT", "fixallbutTNT")
  
  
  ergm.MHP.table("c", "StdNormal", "",  0, "random", "StdNormal")

  ergm.MHP.table("c", "Unif", "",  0, "random", "Unif")
  ergm.MHP.table("c", "Unif", "observed",  0, "random", "UnifNonObserved")
  
  ergm.MHP.table("c", "DiscUnif", "",  0, "random", "DiscUnif")
  ergm.MHP.table("c", "DiscUnif", "observed",  0, "random", "DiscUnifNonObserved")  
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
  ergm.ConstraintImplications("blockdiag", c())
  ergm.ConstraintImplications("hamming", c())
  ergm.ConstraintImplications("observed", c())
}

.RegisterInitMethods <- function(){
  ergm.init.methods("Bernoulli", c("MPLE", "CD", "zeros"))
  ergm.init.methods("StdNormal", c("CD","zeros"))
  ergm.init.methods("Unif", c("CD","zeros"))
  ergm.init.methods("DiscUnif", c("CD","zeros"))
}

