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
      cache <<- rbind.data.frame(cache, data.frame(name=name, short=short, description=description, popular=popular, package=package))
    }
  }
})
