## A "dispatcher" function to distinguish plain ERGMs from separable
## temporal ERGMs.

simulate.formula <- function(object, nsim=1, seed=NULL, ...) {
  nw <- try(ergm.getnetwork(object))
  if(!inherits(nw,"try-error")){
    if("dissolution"%in%names(list(...))){
      simulate.formula.stergm(object, nsim=nsim, seed=seed, ...)
    }else{
      simulate.formula.ergm(object, nsim=nsim, seed=seed, ...)
    }
  }
}


