#  File R/ergm.getnetwork.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2015 Statnet Commons
#######################################################################
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
  if ((dc<-data.class(form)) != "formula")
    stop (paste("Invalid formula of class ",dc))
  trms<-terms(form)
  if (trms[[1]]!="~")
    stop ("Formula must be of the form 'network ~ model'.")

  nw.env<-environment(form)
  if(!exists(x=paste(trms[[2]]),envir=nw.env)){
    stop(paste("The network in the formula '",capture.output(print(form)),"' cannot be found.",sep=""))
  }
  nw <- try(
          {
            tmp <- eval(trms[[2]],envir=nw.env)
            if(is.network(tmp)) tmp else as.network(tmp)
          },
          silent = TRUE)
  if(inherits(nw,"try-error")){
      stop("Invalid network. Is the left-hand-side of the formula correct?")
  }
  if (loopswarning) {
    e <- as.edgelist(nw)
    if(any(e[,1]==e[,2])) {
      print("Warning:  This network contains loops")
    } else if (has.loops(nw)) {
      print("Warning:  This network is allowed to contain loops")
    }
  }
  nw
}
