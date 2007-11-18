ergm.getnetwork <- function (formula) {
  current.warn <- options()$warn
# options(warn=0)
  if ((dc<-data.class(formula)) != "formula")
    stop (paste("Invalid formula of class ",dc))
  trms<-terms(formula)
  if (trms[[1]]!="~")
    stop ("Formula must be of the form 'network ~ model'.")

  nw.env<-environment(formula)
  if(!exists(x=paste(trms[[2]]),envir=nw.env)){
    stop(paste("The network in the formula '",capture.output(print(formula)),"' can not be found.",sep=""))
  }
  nw <- try(as.network(eval(trms[[2]],envir=nw.env), silent = TRUE))  
  if(inherits(nw,"try-error")){
      stop("Invalid network. Is the left-hand-side of the formula correct?")
  }
# options(warn=current.warn)
  nw
}
