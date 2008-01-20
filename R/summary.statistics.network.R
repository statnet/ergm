summary.formula <- function(object, ...){
  current.warn <- options()$warn
  options(warn=0)
  trms<-terms(object)
  if (trms[[1]]!="~")
    stop ("Formula must be of form 'y ~ model'.")
  if(length(trms)<3){return(object)}
  parent <- sys.parent()
  rhs <- try(eval(trms[[2]],parent), silent = TRUE)
  while(inherits(rhs,"try-error") & parent > 1){                      
    parent <- parent - 1
    rhs <- try(eval(trms[[2]],parent), silent = TRUE)
  }
  options(warn=current.warn)
  UseMethod("summary.statistics",object=rhs)
}

summary.statistics <- function(object, ..., drop=FALSE, basis=NULL) {
  UseMethod("summary.statistics")
}

summary.statistics.formula <- function(object, ..., drop=FALSE, basis=NULL) {
  summary.statistics.network(object, ..., drop=drop, basis=basis)
}

summary.statistics.ergm <- function(object, ..., drop=FALSE, basis=NULL)
{
  summary.statistics.network(object$formula, ..., drop=drop, basis=basis)
}

summary.statistics.matrix <- 
summary.statistics.network <- function(object,...,drop=FALSE, basis=NULL) {
  current.warn <- options()$warn
  options(warn=0)
  if(is.network(basis)){
    nw <- basis
    formula <- as.formula(paste("basis",paste(as.character(object)[-2],collapse=" ")))
#   formula <- as.formula(paste(c("nw",as.character(formula)),collapse=" "))
  }else{
    formula <- object
    trms <- terms(formula)
    if(length(trms)>2){
      parent <- sys.parent()
      nw <- try(eval(trms[[2]],parent), silent = TRUE)
      while(inherits(nw,"try-error") & parent > 1){
        parent <- parent - 1
        nw <- try(eval(trms[[2]],parent), silent = TRUE)
      }
      if(class(nw) =="network.series")
        nw <- nw$networks[[1]]
      nw <- as.network(nw, ...)
    }else{
      stop("Must specify a network object")
    }
  }
    
  m <- ergm.getmodel(formula, nw, drop=drop)
  gs <- ergm.getglobalstats(nw, m)
  options(warn=current.warn)
  gs
}




