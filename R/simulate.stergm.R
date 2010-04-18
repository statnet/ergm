simulate.stergm<-function(object,
                          nsim=1, seed=NULL,
                          theta.form=object$coef.form,theta.diss=object$coef.diss,
                        time.burnin=0, time.interval=1, MH.burnin=object$MH.burnin,
                        constraints=object$constraints,
                        stergm.order=object$stergm.order,
                        control=control.simulate.stergm(),
                        toggles=TRUE,
                        verbose=FALSE, ...){
  simulate.formula.stergm(object$formation,object$dissolution,nsim=nsim,seed=seed,theta.form=theta.form, theta.diss=theta.diss,  time.burnin=time.burnin, time.interval=time.interval,MH.burnin=MH.burnin,constraints=constraints,stergm.order=stergm.order,control=control,verbose=verbose,...)
}

simulate.formula.stergm <- function(object, dissolution, nsim=1, seed=NULL,
                                theta.form,theta.diss,
                        time.burnin=0, time.interval=1, MH.burnin=1000,
                        constraints=~.,
                        stergm.order="DissAndForm",
                        control=control.simulate.stergm(),
                        toggles=TRUE,
                        verbose=FALSE, ...) {

  formation <- object
  
  if(!is.null(seed)) set.seed(as.integer(seed))
  nw <- ergm.getnetwork(formation)
  if(class(nw) =="network.series"){
    nw <- nw$networks[[1]]
  }
  nw <- as.network(nw)
  if(!is.network(nw)){
    stop("A network object on the LHS of the formation must be given")
  }

  dissolution<-ergm.update.formula(dissolution,nw~.)
  
  model.form <- ergm.getmodel(formation, nw, drop=control$drop)
  model.diss <- ergm.getmodel(dissolution, nw, stergm.order=stergm.order,drop=control$drop)

  verbose <- match(verbose,
                c("FALSE","TRUE", "very"), nomatch=1)-1
  if(missing(theta.form)) {
    theta.form <- rep(0,length(model.form$coef.names))
    warning("No parameter values given, using Bernouli network.\n\t")
  }

  if(missing(theta.diss)) {
    theta.diss <- rep(0,length(model.diss$coef.names))
    warning("No parameter values given, using Bernoulli dissolution.\nThis means that every time step, half the ties get dissolved!\n\t")
  }

  if((time.burnin!=0 || time.interval!=1) && toggles){
    warning("Burnin is present or interval isn't 1. Toggle list will not be returned.")
    toggles<-FALSE
  }
  
  if(!is.null(seed)) set.seed(as.integer(seed))
    
  MHproposal.form <- MHproposal(constraints,control$prop.args.form,nw,
                                                    model.form,weights=control$prop.weights.form,class="f")
  MHproposal.diss <- MHproposal(constraints,control$prop.args.diss,nw,
                                                    model.diss,weights=control$prop.weights.diss,class="d")
  MCMCparams <- c(control,list(samplesize=nsim, time.interval=time.interval,
                               time.burnin=time.burnin,
                               MH.burnin=MH.burnin,
                               parallel=0,
                               meanstats.form=summary(model.form$formula),
                               meanstats.diss=summary(model.diss$formula),
                               toggles=toggles))
  
  z <- stergm.getMCMCsample(nw, model.form, model.diss,
                             MHproposal.form, MHproposal.diss,
                             theta.form, theta.diss, MCMCparams, verbose)

  if(control$final){
   nw <- z$newnetwork
   return(nw)
  }else{
    changed<-z$changed
    attr(changed,"start")<-1
    attr(changed,"end")<-nsim
    out.list <- list(formation = formation,
                     dissolution = dissolution,
                     networks = nw,
                     constraints=constraints,
                     changed=changed, 
                     maxchanges=z$maxchanges,
                     stats.form = z$statsmatrix.form,stats.diss = z$statsmatrix.diss,
                     coef.form=theta.form,coef.diss=theta.diss)
    class(out.list) <- "network.series"
    return(out.list)
  }
}
