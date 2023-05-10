#  File R/zzz.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2023 Statnet Commons
################################################################################
#' @importFrom Rdpack reprompt
#' @importFrom MASS ginv
#' @importFrom Matrix nearPD
.onAttach <- function(libname, pkgname){
  #' @importFrom statnet.common statnetStartupMessage
  sm <- statnetStartupMessage("ergm", c("statnet","ergm.count","tergm"), TRUE)
  if(!is.null(sm)){
    packageStartupMessage(sm)
    packageStartupMessage(paste(c(strwrap(paste0(sQuote("ergm"), " 4 is a major update that introduces some backwards-incompatible changes. Please type ",sQuote("news(package=\"ergm\")"), " for a list of major changes.")),""),collapse="\n"))
  }
}

.onLoad <- function(libname, pkgname){
  # . is used as a placeholder by stantet.common::NVL3().
  utils::globalVariables(".")

  default_options(ergm.eval.loglik=TRUE,
                  ergm.loglik.warn_dyads=TRUE,
                  ergm.cluster.retries=5)

  eval(COLLATE_ALL_MY_CONTROLS_EXPR)

  .RegisterProposals()
  .RegisterKeywords()
}

## BEGIN boilerplate: should be kept in sync with statnet.common.
# TODO: Figure out some automatic way to keep this in sync with statnet.common.
#' @name snctrl
#'
#' @title Statnet Control
#'
#' @description A utility to facilitate argument completion of control lists, reexported from `statnet.common`.
#'
#' @section Currently recognised control parameters:
#' This list is updated as packages are loaded and unloaded.
#'
#' \Sexpr[results=rd,stage=render]{statnet.common::snctrl_names()}
#'
#' @seealso [statnet.common::snctrl()]
#' @docType import
NULL
#' @export
snctrl <- statnet.common::snctrl
eval(UPDATE_MY_SCTRL_EXPR)
## END boilerplate: should be kept in sync with statnet.common.

.RegisterProposals <- function(){
  ergm_proposal_table("c", "Bernoulli", "|.dyads|bd",  -2, "random", "randomtoggle")
  ergm_proposal_table("c", "Bernoulli", "|.dyads|bd&sparse",  -1, "TNT", "TNT")
  ergm_proposal_table("c", "Bernoulli", "|bdmax|blocks|strat&sparse",  -3, "BDStratTNT", "BDStratTNT")
  ergm_proposal_table("c", "Bernoulli", c("&bdmax|blocks|strat&sparse", "|bdmax&blocks|strat&sparse", "|bdmax|blocks&strat&sparse"),  0, "BDStratTNT", "BDStratTNT")
  ergm_proposal_table("c", "Bernoulli", "", -100, "TNT10", "TNT10")
  ergm_proposal_table("c", "Bernoulli", "&degrees",  0, "random", "CondDegree")
  ergm_proposal_table("c", "Bernoulli", "&degreesmix",  0, "random", "CondDegreeMix")
  ergm_proposal_table("c", "Bernoulli", c("&idegrees&odegrees","&b1degrees&b2degrees"),  0, "random", "CondDegree")
  ergm_proposal_table("c", "Bernoulli", "&odegrees",  0, "random", "CondOutDegree")
  ergm_proposal_table("c", "Bernoulli", "&idegrees",  0, "random", "CondInDegree")
  ergm_proposal_table("c", "Bernoulli", "&b1degrees",  0, "random", "CondB1Degree")
  ergm_proposal_table("c", "Bernoulli", "&b2degrees",  0, "random", "CondB2Degree")
  ergm_proposal_table("c", "Bernoulli", "&degreedist",  0, "random", "CondDegreeDist")
  ergm_proposal_table("c", "Bernoulli", "&idegreedist",  0, "random", "CondInDegreeDist")
  ergm_proposal_table("c", "Bernoulli", "&odegreedist",  0, "random", "CondOutDegreeDist")
  ergm_proposal_table("c", "Bernoulli", "|bd&edges",  0, "random", "ConstantEdges")
  ergm_proposal_table("c", "Bernoulli", "&edges&hamming",  0, "random", "HammingConstantEdges")
  ergm_proposal_table("c", "Bernoulli", "&hamming&sparse",  0, "random", "HammingTNT")
  ergm_proposal_table("c", "Bernoulli", "&dyadnoise&sparse",  1, "TNT", "dyadnoiseTNT")
  ergm_proposal_table("c", "Bernoulli", "&dyadnoise",  0, "random", "dyadnoise")

  ergm_proposal_table("c", "StdNormal", "",  0, "random", "StdNormal")
  ergm_proposal_table("c", "StdNormal", "|.dyads",  0, "random", "DistRLE")

  ergm_proposal_table("c", "Unif", "",  0, "random", "Unif")
  ergm_proposal_table("c", "Unif", "&observed",  0, "random", "UnifNonObserved")
  ergm_proposal_table("c", "Unif", "|.dyads",  0, "random", "DistRLE")

  ergm_proposal_table("c", "DiscUnif", "",  0, "random", "DiscUnif")
  ergm_proposal_table("c", "DiscUnif", "&observed",  0, "random", "DiscUnifNonObserved")

  ergm_proposal_table("c", "DiscUnif", "",  -1, "random2", "DiscUnif2")
  ergm_proposal_table("c", c("Unif","DiscUnif","StdNormal","Poisson","Binomial"), "|.dyads",  -3, "random", "DistRLE")
}


.RegisterKeywords <- function() {
  ergm_keyword(name="binary", short="bin", description="suitable for binary ERGMs", popular=TRUE, package="ergm")
  ergm_keyword(name="bipartite", short="bip", description="suitable for bipartite networks", popular=TRUE, package="ergm")
  ergm_keyword(name="categorical nodal attribute", short="cat nodal attr", description="involves a categorical nodal attribute", popular=FALSE, package="ergm")
  ergm_keyword(name="categorical dyadic attribute", short="cat dyad attr", description="involves a categorical dyadic attribute", popular=FALSE, package="ergm")
  ergm_keyword(name="categorical triadic attribute", short="cat triad attr", description="involves a categorical triadic attribute", popular=FALSE, package="ergm")
  ergm_keyword(name="continuous", short="cont", description="a continuous distribution for edge values", popular=FALSE, package="ergm")
  ergm_keyword(name="curved", short="curved", description="is a curved term", popular=FALSE, package="ergm")
  ergm_keyword(name="directed", short="dir", description="suitable for directed networks", popular=TRUE, package="ergm")
  ergm_keyword(name="discrete", short="discrete", description="a discrete distribution for edge values", popular=FALSE, package="ergm")
  ergm_keyword(name="dyad-independent", short="dyad-indep", description="does not induce dyadic dependence", popular=TRUE, package="ergm")
  ergm_keyword(name="finite", short="fin", description="finite edge values only", popular=FALSE, package="ergm")
  ergm_keyword(name="frequently-used", short="freq", description="is frequently used", popular=FALSE, package="ergm")
  ergm_keyword(name="nonnegative", short="nneg", description="only meaningful for nonnegative edge values", popular=FALSE, package="ergm")
  ergm_keyword(name="operator", short="op", description="a term operator", popular=TRUE, package="ergm")
  ergm_keyword(name="positive", short="pos", description="only meaningful for positive edge values", popular=FALSE, package="ergm")
  ergm_keyword(name="quantitative nodal attribute", short="quant nodal attr", description="involves a quantitative nodal attribute", popular=FALSE, package="ergm")
  ergm_keyword(name="quantitative dyadic attribute", short="quant dyad attr", description="involves a quantitative dyadic attribute", popular=FALSE, package="ergm")
  ergm_keyword(name="quantitative triadic attribute", short="quant triad attr", description="involves a quantitative triadic attribute", popular=FALSE, package="ergm")
  ergm_keyword(name="soft", short="soft", description="a constraint that does not necessarily forbid specific networks outright but reweights their probabilities", popular=FALSE, package="ergm")
  ergm_keyword(name="triad-related", short="triad rel", description="involves triangles, two-paths, and other triadic structures", popular=FALSE, package="ergm")
  ergm_keyword(name="valued", short="val", description="suitable for valued ERGMs", popular=TRUE, package="ergm")
  ergm_keyword(name="undirected", short="undir", description="suitable for undirected networks", popular=TRUE, package="ergm")
}

.onUnload <- function(libpath){
  library.dynam.unload("ergm",libpath)
}
