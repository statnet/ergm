simulatedyn <- function(object, dissolve=NULL, nsteps=1, seed=NULL, theta0,gamma0,
                        burnin=0, interval=1, dyninterval=1000,
                        constraint="none",
                        dissolve.order="DissThenForm",
                        control=ergm.simulatedyn.control(),
                        verbose=FALSE) {
  formula <- object

  control$dyninterval<-dyninterval
  
  if(is.null(seed)){seed <- sample(10000000, size=1)}
  set.seed(as.integer(seed))
  nw <- ergm.getnetwork(formula)
  if(class(nw) =="network.series"){
    nw <- nw$networks[[1]]
  }
  nw <- as.network(nw)
  if(!is.network(nw)){
    stop("A network object on the LHS of the formula must be given")
  }
  # Resolve conditioning
  # This is laborious to cover the various partial specifications
  #
  BD <- ergm.boundDeg(control$boundDeg, nnodes=network.size(nw))

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
  
  if(is.null(seed)){seed <- sample(10000000, size=1)}
  set.seed(as.integer(seed))
    
  MHproposal.form <- getMHproposal(lookupMHproposal("f",constraint,control$prop.weights.form),
                                   control$prop.args.form, nw, model.form)
  MHproposal.diss <- getMHproposal(lookupMHproposal("d",constraint,control$prop.weights.diss),
                                   control$prop.args.diss, nw, model.diss)
  MCMCparams <- c(control,list(nsteps=nsteps, interval=interval,
                           stats.form=matrix(summary(model.form$formula),ncol=length(model.form$coef.names),nrow=1),
                           stats.diss=matrix(summary(model.diss$formula),ncol=length(model.diss$coef.names),nrow=1),
                           burnin=burnin,
                           parallel=0,
                           meanstats.form=theta0-theta0,
                           meanstats.diss=gamma0-gamma0))
  
  z <- ergm.getMCMCDynsample(nw, model.form, model.diss,
                             MHproposal.form, MHproposal.diss,
                             theta0, gamma0, MCMCparams, verbose, BD)

  if(control$final){
   nw <- z$newnetwork
   return(nw)
  }else{
    out.list <- list(formula = formula, networks = nw,
                     changed=z$changed, 
                     maxchanges=z$maxchanges,
                     stats.form = z$statsmatrix.form,stats.diss = z$statsmatrix.diss,
                     coef.form=theta0,coef.diss=gamma0)
    class(out.list) <- "network.series"
    return(out.list)
  }
}
