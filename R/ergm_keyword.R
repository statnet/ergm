#  File R/ergm_keyword.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################

#' Dynamic ERGM keyword registry
#'
#' A function to manage dynamic ERGM keywords. To register a keyword, call the function with all parameters
#' provided. To fetch all registered keywords, call the function with no parameters specified.
#'
#' @param name full name of the keyword
#' @param short abbreviation of the keyword name
#' @param description description of the keyword
#' @param popular logical to indicate if a keyword is popular
#' @param package package the keyword is first defined in
#' @return Returns a dataframe with the following columns:
#'   - name
#'   - short
#'   - description
#'   - popular
#'   - package
#' @keywords models internal
#' @export
ergm_keyword <- local({
  cache <- data.frame(name=c(), short=c(), description=c(), popular=c(), package=c())

  function(name=NULL, short=NULL, description=NULL, popular=NULL, package=NULL) {
    if (all(is.null(name), is.null(short), is.null(description), is.null(popular), is.null(package))) {
      return(cache)
    } else if (any(is.null(name), is.null(short), is.null(description), is.null(popular), is.null(package))) {
      stop("All arguments are needed to register Ergm keyword")
    } else if (!is.logical(popular)) {
      stop("Logical value expected for argument 'popular'")
    } else {
      cache <<- rbind.data.frame(cache, data.frame(name=name, short=short, description=description, popular=popular, package=package, stringsAsFactors=FALSE))
    }
  }
})

.RegisterKeywords <- function() {
  ergm_keyword(name="binary", short="bin", description="", popular=FALSE, package="ergm")
  ergm_keyword(name="bipartite", short="bip", description="", popular=FALSE, package="ergm")
  ergm_keyword(name="categorical nodal attribute", short="cat nodal attr", description="", popular=FALSE, package="ergm")
  ergm_keyword(name="continuous", short="cont", description="", popular=FALSE, package="ergm")
  ergm_keyword(name="curved", short="curved", description="", popular=FALSE, package="ergm")
  ergm_keyword(name="directed", short="dir", description="", popular=FALSE, package="ergm")
  ergm_keyword(name="discrete", short="discrete", description="", popular=FALSE, package="ergm")
  ergm_keyword(name="dyad-independent", short="dyad-indep", description="", popular=FALSE, package="ergm")
  ergm_keyword(name="finite", short="fin", description="", popular=FALSE, package="ergm")
  ergm_keyword(name="frequently-used", short="freq", description="", popular=FALSE, package="ergm")
  ergm_keyword(name="non-negative", short="non-neg", description="", popular=FALSE, package="ergm")
  ergm_keyword(name="operator", short="op", description="", popular=FALSE, package="ergm")
  ergm_keyword(name="positive", short="pos", description="", popular=FALSE, package="ergm")
  ergm_keyword(name="quantative nodal attribute", short="quant nodal attr", description="", popular=FALSE, package="ergm")
  ergm_keyword(name="soft", short="soft", description="", popular=FALSE, package="ergm")
  ergm_keyword(name="triad-related", short="triad rel", description="", popular=FALSE, package="ergm")
  ergm_keyword(name="valued", short="val", description="", popular=TRUE, package="ergm")
  ergm_keyword(name="undirected", short="undir", description="", popular=FALSE, package="ergm")
}
