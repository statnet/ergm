#  File ergm/R/stergm.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
################################################################################
# stergm --- fit Separable Temporal ERGMs.
################################################################################

stergm <- function(nw, formation, dissolution, estimate, times=NULL, offset.coef.form=NULL, offset.coef.diss=NULL,
                   targets=NULL, target.stats=NULL,
                   eval.loglik=FALSE,
                 control=control.stergm(),
                 verbose=FALSE, ...) {
  current.warn <- options()$warn
  on.exit(options(warn=current.warn), add=TRUE)
  options(warn=0)
  if(!is.null(control$seed))  set.seed(as.integer(control$seed))

  estimate <- match.arg(estimate,c("CMLE","CMPLE","EGMME"))
  
  if(!inherits(formation,"formula") || !inherits(dissolution,"formula"))
    stop("Arguments formation and dissolution must be formulas.")

  if(length(formation)==3){
    warning("Formation formula has an LHS, which will be ignored in favor of nw.")
    formation <- formation[c(1,3)] # in a formula f<-y~x, f[1]=~, f[2]=y, and f[3]=x
  }

  if(length(dissolution)==3){
    warning("Dissolution formula has an LHS, which will be ignored in favor of nw.")
    dissolution <- dissolution[c(1,3)] # in a formula f<-y~x, f[1]=~, f[2]=y, and f[3]=x
  }

  out <- switch(estimate,
                CMLE=,
                CMPLE=stergm.CMLE(nw, formation, dissolution,
                  times, offset.coef.form, offset.coef.diss, eval.loglik,
                  estimate, control, verbose),
                EGMME=stergm.EGMME(nw, formation, dissolution, offset.coef.form, offset.coef.diss,
                  targets, target.stats, estimate, control, verbose)
                  )

  
  out$formation <- formation
  out$dissolution <- dissolution
  out$control <- control

  class(out)<-"stergm"
  out
}
