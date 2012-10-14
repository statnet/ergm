.onAttach <- function(lib, pkg){
  packageStartupMessage(mkStartupMessage("ergm"))
  
  # If the following have already been defined in the latentnet package, don't duplicate. Otherwise, assign them.
  IFNOTEXISTS <- c("robust.inverse","mcmc.diagnostics","mcmc.diagnostics.default","gof","gof.default")
  for(fun in IFNOTEXISTS){
    if(!exists(fun, mode="function")){
      assign(fun, get(paste('.',fun,sep='')), pos="package:ergm")
    }
  }
  
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
  ergm.MHP.table("c", "Bernoulli", "observed",  0, "random", "randomtoggleNonObserved")
  ergm.MHP.table("c", "Bernoulli", "blockdiag", 0, "random", "blockdiag")
  ergm.MHP.table("c", "Bernoulli", "blockdiag", 1, "TNT", "blockdiagTNT")
  ergm.MHP.table("c", "Bernoulli", "bd+blockdiag", 0, "random", "blockdiag")
  ergm.MHP.table("c", "Bernoulli", "bd+blockdiag", 1, "TNT", "blockdiagTNT")
  ergm.MHP.table("c", "Bernoulli", "blockdiag+observed",  0, "random", "blockdiagNonObserved")
  ergm.MHP.table("c", "Bernoulli", "bd+blockdiag+observed",  0, "random", "blockdiagNonObserved")
}

.RegisterConstraintImplications <- function(){
  ergm.ConstraintImplications("edges", c())
  ergm.ConstraintImplications("degrees", c("edges", "idegrees", "odegrees", "idegreedist", "odegreedist", "degreedist", "bd"))
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
  ergm.init.methods("Bernoulli", c("MPLE", "zeros"))
}
