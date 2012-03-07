#  File ergm/R/logLik.stergm.CMLE.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################

logLik.stergm.CMLE<-function(object, add=FALSE, force.reeval=FALSE, eval.loglik=add || force.reeval, control=control.logLik.stergm(), ...){
  if(!is.null(control$seed))  set.seed(as.integer(control$seed))

  if(add){
    object$formation.fit <- logLik(object$formation.fit, add=TRUE, force.reeval=FALSE, eval.loglik = add || force.reeval, control=control$control.form)
    object$dissolution.fit <- logLik(object$dissolution.fit, add=TRUE, force.reeval=FALSE, eval.loglik = add || force.reeval, control=control$control.diss)
    
    object
  }else{
    llk.form <- logLik(object$formation.fit, add=FALSE, force.reeval=FALSE, eval.loglik = add || force.reeval, control=control$control.form)
    llk.diss <- logLik(object$dissolution.fit, add=FALSE, force.reeval=FALSE, eval.loglik = add || force.reeval, control=control$control.diss)

    llk <- llk.form + llk.diss
    attr(llk,"df") <- attr(llk.form,"df") + attr(llk.diss,"df")
    attr(llk,"nobs") <- attr(llk.form,"nobs") + attr(llk.diss,"nobs")
    
    llk
  }
}
