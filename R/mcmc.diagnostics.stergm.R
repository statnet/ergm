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
