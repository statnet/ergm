#  File R/WtReferences.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
# For now, this file contains information about reference
# measures. Eventually, we should create an "InitErgmReference" or similar
# framework.



#' Set up the initial fitting methods for reference measure and query available
#' methods for that reference measure
#' 
#' This is a low-level function not intended to be called directly by end
#' users. This function sets up the available initial fitting methods for each
#' reference measure and queries them.
#' 
#' 
#' @param reference The reference measure used in the model.
#' @param new.methods If passed, prepends the new initial fitting methods to
#' the list for that reference measure.
#' @return A character vector listing initial methods for the reference measure
#' specified. (If `new.methods` is passed, does so invisibly.)
#' @export
ergm.init.methods <- local({
  init.methods <- list()
  function(reference, new.methods){
    if(!missing(new.methods)){
      init.methods[[reference]] <<- unique(c(new.methods, init.methods[[reference]]))
    }else{
      init.methods[[reference]]
    }
  }
})
