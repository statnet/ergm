#################################################################################
# The <ergm.getnetwork> function ensures that the network in a given formula
# is valid; if so, the network is returned; if not, execution is halted with
# warnings
#
# --PARAMETERS--
#   formula     :  the formula as 'network ~ model.term(s)'
#   loopswarning:  whether warnings about loops should be printed (T or F);
#                  default=TRUE
#
# --RETURNED--
#   nw: the network from the formula IF (i) the formula was correctly structured
#       and (ii) the network is found within the formula's enviornment
#
###################################################################################

ergm.getnetwork <- function (form, loopswarning=TRUE) {
  current.warn <- options()$warn
# options(warn=0)
  if ((dc<-data.class(form)) != "formula")
    stop (paste("Invalid formula of class ",dc))
  trms<-terms(form)
  if (trms[[1]]!="~")
    stop ("Formula must be of the form 'network ~ model'.")

  nw.env<-environment(form)
  if(!exists(x=paste(trms[[2]]),envir=nw.env)){
    stop(paste("The network in the formula '",capture.output(print(form)),"' can not be found.",sep=""))
  }
  nw <- try(as.network(eval(trms[[2]],envir=nw.env), silent = TRUE))  
  if(inherits(nw,"try-error")){
      stop("Invalid network. Is the left-hand-side of the formula correct?")
  }
  # options(warn=current.warn)
  if (loopswarning) {
    e <- as.matrix.network.edgelist(nw)
    if(any(e[,1]==e[,2])) {
      print("Warning:  This network contains loops")
    } else if (has.loops(nw)) {
      print("Warning:  This network is allowed to contain loops")
    }
  }
  nw
}
