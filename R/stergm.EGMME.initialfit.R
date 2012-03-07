#  File ergm/R/stergm.EGMME.initialfit.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
stergm.EGMME.initialfit<-function(init.form, init.diss, nw, model.form, model.diss, model.mon, control, verbose=FALSE){
  if(!is.null(control$init.method) && control$init.method == "zeros"){
    init.form[is.na(init.form)]<-0
    init.diss[is.na(init.diss)]<-0
  }else if(!any(is.na(init.form)) && !any(is.na(init.diss))){
    # Don't need to do anything.
  }else if(all(model.form$formula[c(1,3)]==model.mon$formula[c(1,3)])
           && all(model.diss$etamap$offsettheta)
           && all(model.diss$coef.names %in% model.form$coef.names)
           && is.dyad.independent(model.diss$formula)){
    if(verbose) cat("Formation statistics are analogous to targeted statistics, dissolution is fixed, dissolution terms appear to have formation analogs, and dissolution process is dyad-independent, so using Carnegie-Krivitsky-Hunter-Goodreau approximation.\n")
    # Fit an ERGM to the formation terms:
    init.form<-coef(ergm(model.form$formula, eval.loglik=FALSE, control=control.ergm(init=init.form)))
    # Now, match up non-offset formation terms with dissolution terms.
    # In case it's not obvious (it's not to me) what the following
    # does, it takes non-offset elements of init.form, then, from
    # those, it takes those elements that correspond (by name) to the
    # dissolution coefficient names and decrements them by init.diss.
    #
    # Yes, I am also a little surprised that assigning to a
    # double-index works.
    init.form[!model.form$etamap$offsettheta][match(names(init.diss),names(init.form[!model.form$etamap$offsettheta]))] <-
      init.form[!model.form$etamap$offsettheta][match(names(init.diss),names(init.form[!model.form$etamap$offsettheta]))] - init.diss
  }else{
    stop("No initial parameter method for specified model and targets combination is implemented. Specify via control$init.form and control$init.diss .")
  }
  out <- list(formation = model.form$formula, dissolution = model.diss$formula, targets = model.mon$formula, target.stats=model.mon$target.stats, nw = nw, control = control, formation.fit = list(coef=init.form, etamap = model.form$etamap), dissolution.fit = list(coef=init.diss, etamap = model.diss$etamap))
  class(out)<-"stergm"
  out
}
