##########################################################################
# The <ergm.degenerate> function ensures that the given formula specifies
# a valid network on its left hand side (by halting if not) and determines 
# which of the model terms would be dropped because of their degeneracy
# status
#
# --PARAMETERS--
#   formula:  the ergm formula
#   ...    :  addtional parameters that may be passed from within;
#             these are ignored
#
# --RETURNED--
#   a logical vector of whether each model term would be dropped
#   because it is deemed degenerate
#
#########################################################################

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
#   Check for degenerate specifications
    m.all <- ergm.getmodel(trms, g)
    degenerate <- +(m.all$maxval==obs.stats)-(m.all$minval==obs.stats)
    names(degenerate) <- m.all$coef.names
    degenerate
}
