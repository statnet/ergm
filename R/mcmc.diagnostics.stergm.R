#  File ergm/R/mcmc.diagnostics.stergm.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
mcmc.diagnostics.stergm <- function(object, 
                                    center=TRUE,
                                    curved=TRUE,
                                    vars.per.page=3, ...){
  if(!is.null(object$formation.fit$sample)){
    cat("\n==========================\n")
    cat("Formation fit diagnostics\n")
    cat("==========================\n\n")
    mcmc.diagnostics.ergm(object$formation.fit, center=center, curved=curved, vars.per.page=vars.per.page, ...)
  }
  if(!is.null(object$dissolution.fit$sample)){
    cat("\n==========================\n")
    cat("Dissolution fit diagnostics\n")
    cat("==========================\n\n")
    mcmc.diagnostics.ergm(object$dissolution.fit, center=center, curved=curved, vars.per.page=vars.per.page, ...)
  }
  if(!is.null(object$sample)){
    cat("\n==========================\n")
    cat("EGMME diagnostics\n")
    cat("==========================\n\n")
    mcmc.diagnostics.ergm(object, center=center, curved=curved, vars.per.page=vars.per.page, ...)
  }
}
