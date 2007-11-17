simulatedyn <- function(object, dissolve=NULL, nsteps=1, seed=NULL, theta0,gamma0,
                        burnin=0, interval=1, dyninterval=1000,
                        constraints=~.,
                        dissolve.order="DissThenForm",
                        control=ergm.simulatedyn.control(),
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

  model.form <- ergm.getmodel(formula, nw, drop=control$drop)
  model.diss <- ergm.getmodel.dissolve(dissolve, nw, dissolve.order=dissolve.order)

  verbose <- match(verbose,
                c("FALSE","TRUE", "very"), nomatch=1)-1
  if(missing(theta0)) {
    theta0 <- rep(0,length(model.form$coef.names))
    warning("No parameter values given, using Bernouli network.\n\t")
  }

  if(missing(gamma0)) {
    gamma0 <- rep(0,length(model.diss$coef.names))
    warning("No parameter values given, using Bernoulli dissolution.\nThis means that every time step, half the ties get dissolved!\n\t")
  }

  if(burnin!=0 || interval!=1) warning("Burnin is present or interval isn't 1. Toggle list will not be returned.")
  
  if(!is.null(seed)) set.seed(as.integer(seed))
    
  MHproposal.form <- getMHproposal(constraints,control$prop.args.form,nw,
                                                    model.form,weights=control$prop.weights.form,class="f")
  MHproposal.diss <- getMHproposal(constraints,control$prop.args.diss,nw,
                                                    model.diss,weights=control$prop.weights.diss,class="d")
  MCMCparams <- c(control,list(nsteps=nsteps, interval=interval,
                           burnin=burnin,
                           parallel=0,
                           meanstats.form=theta0-theta0,
                           meanstats.diss=gamma0-gamma0))
  
  z <- ergm.getMCMCDynsample(nw, model.form, model.diss,
                             MHproposal.form, MHproposal.diss,
                             theta0, gamma0, MCMCparams, verbose)

  if(control$final){
   nw <- z$newnetwork
   return(nw)
  }else{
    out.list <- list(formula = formula, networks = nw,
                     constraints=constraints,
                     changed=z$changed, 
                     maxchanges=z$maxchanges,
                     stats.form = z$statsmatrix.form,stats.diss = z$statsmatrix.diss,
                     coef.form=theta0,coef.diss=gamma0)
    class(out.list) <- "network.series"
    return(out.list)
  }
}
