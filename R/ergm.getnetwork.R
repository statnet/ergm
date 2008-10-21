#  File ergm/R/ergm.getnetwork.R
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
# Copyright 2003 Mark S. Handcock, University of Washington
#                David R. Hunter, Penn State University
#                Carter T. Butts, University of California - Irvine
#                Steven M. Goodreau, University of Washington
#                Martina Morris, University of Washington
# Copyright 2007 The statnet Development Team
######################################################################
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
