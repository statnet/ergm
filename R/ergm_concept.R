#  File R/ergm_concept.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################

ergm_concept <- local({
  cache <- data.frame(name=c(), short=c(), description=c(), popular=c(), package=c())

  function(name=NULL, short=NULL, description=NULL, popular=NULL, package=NULL) {
    if (all(is.null(name), is.null(short), is.null(description), is.null(popular), is.null(package))) {
      return(cache)
    } else if (any(is.null(name), is.null(short), is.null(description), is.null(popular), is.null(package))) {
      stop("All arguments are needed to register Ergm concept")
    } else if (!is.logical(popular)) {
      stop("Logical value expected for argument 'popular'")
    } else {
      cache <<- rbind.data.frame(cache, data.frame(name=name, short=short, description=description, popular=popular, package=package, stringsAsFactors=FALSE))
    }
  }
})

.RegisterConcepts <- function() {
  ergm_concept(name="binary", short="bin", description="", popular=FALSE, package="ergm")
  ergm_concept(name="bipartite", short="bip", description="", popular=FALSE, package="ergm")
  ergm_concept(name="categorical nodal attribute", short="cat nodal attr", description="", popular=FALSE, package="ergm")
  ergm_concept(name="continuous", short="", description="", popular=FALSE, package="ergm")
  ergm_concept(name="curved", short="curved", description="", popular=FALSE, package="ergm")
  ergm_concept(name="directed", short="dir", description="", popular=FALSE, package="ergm")
  ergm_concept(name="discrete", short="", description="", popular=FALSE, package="ergm")
  ergm_concept(name="dyad-independent", short="dyad-indep", description="", popular=FALSE, package="ergm")
  ergm_concept(name="finite", short="", description="", popular=FALSE, package="ergm")
  ergm_concept(name="frequently-used", short="freq", description="", popular=FALSE, package="ergm")
  ergm_concept(name="non-negative", short="non-neg", description="", popular=FALSE, package="ergm")
  ergm_concept(name="operator", short="op", description="", popular=FALSE, package="ergm")
  ergm_concept(name="positive", short="", description="", popular=FALSE, package="ergm")
  ergm_concept(name="quantative nodal attribute", short="quant nodal attr", description="", popular=FALSE, package="ergm")
  ergm_concept(name="soft", short="", description="", popular=FALSE, package="ergm")
  ergm_concept(name="triad-related", short="triad rel", description="", popular=FALSE, package="ergm")
  ergm_concept(name="valued", short="val", description="", popular=TRUE, package="ergm")
  ergm_concept(name="undirected", short="undir", description="", popular=FALSE, package="ergm")
}
