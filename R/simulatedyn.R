simulatedyn <- function(object, dissolve, nsteps=1, seed=NULL, theta,gamma,
                        burnin=0, interval=1, dyninterval=1000,
                        constraints=~.,
                        dissolve.order="DissThenForm",
                        control=control.simulatedyn(),
                        toggles=TRUE,
                        verbose=FALSE) {
  formula <- object

  control$dyninterval<-dyninterval
  
  if(!is.null(seed)) set.seed(as.integer(seed))
  nw <- ergm.getnetwork(formula)
  if(class(nw) =="network.series"){
    nw <- nw$networks[[1]]
  }
  nw <- as.network(nw)
  if(!is.network(nw)){
    stop("A network object on the LHS of the formula must be given")
  }

  dissolve<-safeupdate.formula(dissolve,nw~.)
  
  model.form <- ergm.getmodel(formula, nw, drop=control$drop)
  model.diss <- ergm.getmodel(dissolve, nw, dissolve.order=dissolve.order)

  verbose <- match(verbose,
                c("FALSE","TRUE", "very"), nomatch=1)-1
  if(missing(theta)) {
    theta <- rep(0,length(model.form$coef.names))
    warning("No parameter values given, using Bernouli network.\n\t")
  }

  if(missing(gamma)) {
    gamma <- rep(0,length(model.diss$coef.names))
    warning("No parameter values given, using Bernoulli dissolution.\nThis means that every time step, half the ties get dissolved!\n\t")
  }

  if((burnin!=0 || interval!=1) && toggles){
    warning("Burnin is present or interval isn't 1. Toggle list will not be returned.")
    toggles<-FALSE
  }
  
  if(!is.null(seed)) set.seed(as.integer(seed))
    
  MHproposal.form <- MHproposal(constraints,control$prop.args.form,nw,
                                                    model.form,weights=control$prop.weights.form,class="f")
  MHproposal.diss <- MHproposal(constraints,control$prop.args.diss,nw,
                                                    model.diss,weights=control$prop.weights.diss,class="d")
  MCMCparams <- c(control,list(samplesize=nsteps, interval=interval,
                               burnin=burnin,
                               parallel=0,
                               meanstats.form=summary(model.form$formula),
                               meanstats.diss=summary(model.diss$formula),
                               toggles=toggles))
  
  z <- ergm.getMCMCDynsample(nw, model.form, model.diss,
                             MHproposal.form, MHproposal.diss,
                             theta, gamma, MCMCparams, verbose)

  if(control$final){
   nw <- z$newnetwork
   return(nw)
  }else{
    changed<-z$changed
    attr(changed,"start")<-1
    attr(changed,"end")<-nsteps
    out.list <- list(formula = formula, networks = nw,
                     constraints=constraints,
                     changed=changed, 
                     maxchanges=z$maxchanges,
                     stats.form = z$statsmatrix.form,stats.diss = z$statsmatrix.diss,
                     coef.form=theta,coef.diss=gamma)
    class(out.list) <- "network.series"
    return(out.list)
  }
}
