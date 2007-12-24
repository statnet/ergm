ergm.degenerate <- function(formula, ...)
{
    trms <- ergm.getterms(formula)
#   termnames <- ergm.gettermnames(trms)
    parent <- sys.parent()
    g <- try(as.network(eval(trms[[2]],parent)), silent=TRUE)
    while(inherits(g,"try-error") & parent > 1){
      parent <- parent - 1
      g <- try(as.network(eval(trms[[2]],parent)), silent=TRUE)
    }
    if(inherits(g,"try-error")){
     stop("Invalid graph. Is the left-hand-side of the formula correct?")
    }
#
#   Check for degenerate specifications
#
    m.all <- ergm.getmodel(trms, g, drop=FALSE)
    m <- ergm.getmodel(trms, g, drop=TRUE)
    degenerate <- is.na(match(m.all$coef.names, m$coef.names))
    names(degenerate) <- m.all$coef.names
    degenerate
}
