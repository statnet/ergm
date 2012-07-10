.onAttach <- function(lib, pkg){
  info <- packageDescription("ergm")

  # If the following have already been defined in the latentnet package, don't duplicate. Otherwise, assign them.
  IFNOTEXISTS <- c("robust.inverse","mcmc.diagnostics","mcmc.diagnostics.default","gof","gof.default","paste.and")
  for(fun in IFNOTEXISTS){
    if(!exists(fun, mode="function")){
      assign(fun, get(paste('.',fun,sep='')), pos="package:ergm")
    }
  }
  
  packageStartupMessage(
    paste('\nergm: version ', info$Version, ', created on ', info$Date, '\n',
          "Copyright (c) 2003, Mark S. Handcock, University of California-Los Angeles\n",
          "                    David R. Hunter, Penn State University\n",
          "                    Carter T. Butts, University of California-Irvine\n",
          "                    Steven M. Goodreau, University of Washington\n",
          "                    Pavel N. Krivitsky, Penn State University\n",
          "                    Martina Morris, University of Washington\n",
          'Based on "statnet" project software (statnet.org).\n',
          'For license and citation information see statnet.org/attribution\n',
          'or type citation("ergm").\n', sep="")
                        )
  
  .RegisterMHPs()
  .RegisterConstraintImplications()
}

.RegisterMHPs <- function(){
  ergm.MHP.table("f", "Bernoulli", "",  0, "random", "formation")
  ergm.MHP.table("f", "Bernoulli", "bd",  0, "random", "formation")
  ergm.MHP.table("f", "Bernoulli", "",  1, "TNT", "formationTNT")
  ergm.MHP.table("f", "Bernoulli", "bd",  1, "TNT", "formationTNT")
  ergm.MHP.table("d", "Bernoulli", "",  0, "random", "dissolution")
  ergm.MHP.table("d", "Bernoulli", "bd",  0, "random", "dissolution")
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
  ergm.MHP.table("c", "Bernoulli", "observed",  0, "random", "randomtoggleNonObserved")
  ergm.MHP.table("c", "Poisson", "",  0, "random", "Poisson")
  ergm.MHP.table("c", "Poisson", "",  0, "0inflated", "ZIPoisson")
  ergm.MHP.table("c", "Poisson", "observed",  0, "random", "PoissonNonObserved")
  ergm.MHP.table("c", "DescRank", "",  0, "random", "DescRank")
  ergm.MHP.table("c", "DescRank", "ranks",  0, "random", "DescRankEquivalent")
  ergm.MHP.table("c", "StdNormal", "",  0, "random", "StdNormal")
  ergm.MHP.table("c", "StdNormal", "ranks",  0, "random", "StdNormalRank")
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
  ergm.ConstraintImplications("observed", c())
}
