ergm <- function(formula, theta0="MPLE", 
                 MPLEonly=FALSE, MLestimate=!MPLEonly, seed=NULL,
                 burnin=10000, MCMCsamplesize=10000, interval=100, maxit=3,
                 proposaltype="randomtoggle", proposalargs=NULL,
                 meanstats=NULL,
                 dissolve=NULL, gamma=-4.59512, dissolve.order="DissThenForm",
                 algorithm.control=list(),
                 verbose=FALSE, ...) {
  current.warn <- options()$warn
  options(warn=0)
  if(!is.null(seed))  set.seed(as.integer(seed))
  if (verbose) cat("Evaluating network in model\n")

  ## Defaults :
  con <- list(nr.maxit=100, calc.mcmc.se=TRUE, hessian=FALSE, 
              compress=FALSE,
              maxNumDyadTypes=10000, 
              maxedges=20000,
              maxchanges=1000000,
              MPLEsamplesize=50000, 
              trace=0,
              boundDeg=NULL,
              steplength=0.5,
              drop=TRUE,
              force.mcmc=FALSE,
              mcmc.precision=0.05,
              metric="Likelihood",
              method="BFGS",
              trustregion=20,
              style="Newton-Raphson",
              phase1_n=NULL, initial_gain=NULL, 
              nsubphases=maxit, niterations=NULL, phase3_n=NULL,
              dyninterval=1000,
              parallel=0,
              returnMCMCstats=TRUE
             )

  con[(namc <- names(algorithm.control))] <- algorithm.control

  nw <- ergm.getnetwork(formula)
  if(!is.null(meanstats)){ con$drop <- FALSE }
  if (verbose) cat("Fitting initial model.\n")
   if(!is.null(dissolve)){
       proposaltype <- "formationTNT"
   }
  model.initial <- ergm.getmodel(formula, nw, drop=con$drop, initialfit=TRUE)
  MHproposal <- getMHproposal(proposaltype, proposalargs, nw, model.initial)
#
  BD <- ergm.boundDeg(con$boundDeg, nnodes=network.size(nw))
  Clist.initial <- ergm.Cprepare(nw, model.initial)
  Clist.initial$meanstats=meanstats
  theta0copy <- theta0
  initialfit <- ergm.initialfit(theta0copy, MLestimate, Clist.initial,
                                model.initial, verbose=verbose, ...)
  if (MLestimate && 
      (   !ergm.independencemodel(model.initial)
       || !is.null(meanstats))
       || con$force.mcmc
      ) {
    theta0 <- initialfit$coef
    names(theta0) <- model.initial$coef.names
    theta0[is.na(theta0)] <- 0
  } else { # Just return initial (non-MLE) fit and exit.
    initialfit$network <- nw
    initialfit$newnetwork <- nw
    initialfit$formula <- formula
    return(initialfit)
  } 
  model <- ergm.getmodel(formula, nw, drop=con$drop, expanded=TRUE)
  theta0 <- ergm.revisetheta0(model, theta0)
  # revise theta0 to reflect additional parameters

  Clist <- ergm.Cprepare(nw, model)
  Clist$meanstats <- meanstats
  Clist$obs <- summary(model$formula)

  if (verbose) cat("ergm.mainfitloop\n")
  MCMCparams=c(con,list(samplesize=MCMCsamplesize, burnin=burnin, interval=interval,maxit=maxit))
  MHproposal=list(package=con$proposalpackage, type=proposaltype)
  styles <- c("Newton-Raphson","Robbins-Monro","Stochastic-Approximation")
  con$style <- styles[pmatch(con$style,styles,nomatch=1)]
  if(!is.null(dissolve)){
   if (verbose) cat("Fitting Dynamic ERGM.\n")
   model.dissolve <- ergm.getmodel.dissolve(dissolve, nw, dissolve.order)
   v <- switch(con$style,
    "Robbins-Monro" = ergm.robmon.dyn(theta0, nw, model, model.dissolve,
                    Clist, BD, gamma, 
                    MCMCparams=MCMCparams, MHproposal=MHproposal,
                    verbose),
                      ergm.mainfitloop.dyn(theta0, nw,
                          model, model.dissolve, Clist,
                          BD, gamma, initialfit,
                          MCMCparams=MCMCparams, MHproposal=MHproposal,
                          verbose=verbose, 
                          ...)
              )
  }else{
   if (verbose) cat("Fitting ERGM.\n")
   v <- switch(con$style,
    "Robbins-Monro" = ergm.robmon(theta0, nw, model, Clist, BD, burnin, interval,
                      proposaltype, verbose, con),
    "Stochastic-Approximation" = ergm.stocapprox(theta0, nw, model, 
                                 Clist, BD, 
                                 MCMCparams=MCMCparams, MHproposal=MHproposal,
                                 verbose),
                      ergm.mainfitloop(theta0, nw,
                          model, Clist,
                          BD, initialfit,
                          MCMCparams=MCMCparams, MHproposal=MHproposal,
                          verbose=verbose, 
                          ...)
              )
  }

  if (!con$returnMCMCstats)
    v$sample <- NULL
  options(warn=current.warn)
  v
}
