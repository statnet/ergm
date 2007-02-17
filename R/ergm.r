ergm <- function(formula, theta0="MPLE", 
                 MPLEonly=FALSE, MLestimate=!MPLEonly, seed=NULL,
                 burnin=10000, MCMCsamplesize=10000, interval=100, maxit=3,
                 proposaltype="randomtoggle", 
                 meanstats=NULL,
                 dissolve=NULL, gamma=0.01,
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
              MPLEsamplesize=50000, 
              trace=0,
              boundDeg=NULL,
              steplength=0.5,
              drop=TRUE,
              proposalpackage="statnet",
              mcmc.precision=0.05,
              metric="Likelihood",
              method="BFGS",
              trustregion=20,
              style="Newton-Raphson",
              phase1_n=NULL, initial_gain=NULL, 
              nsubphases=maxit, niterations=NULL, phase3_n=NULL,
              dyninterval=1000,
              returnMCMCstats=TRUE
             )

  con[(namc <- names(algorithm.control))] <- algorithm.control

  nw <- ergm.getnetwork(formula)
  if(!is.null(meanstats)){ con$drop <- FALSE }
  if (verbose) cat("Fitting initial model.\n")
  if(is.bipartite(nw)){
   if(proposaltype=="randomtoggle"){proposaltype <- "Bipartiterandomtoggle"}
   if(proposaltype=="TNT"){proposaltype <- "BipartiteTNT"}
   if(proposaltype=="ConstantEdges"){proposaltype <- "BipartiteConstantEdges"}
   if(proposaltype=="HammingConstantEdges"){proposaltype <- "BipartiteHammingConstantEdges"}
   if(proposaltype=="Hamming"){proposaltype <- "BipartiteHamming"}
   if(proposaltype=="formation"){proposaltype <- "BipartiteFormation"}
   if(proposaltype=="formationTNT"){proposaltype <- "BipartiteFormationTNT"}
   if(!is.null(dissolve)){
     proposaltype <- "BipartiteFormationTNT"
   }
  }
  model.initial <- ergm.getmodel(formula, nw, drop=con$drop, initialfit=TRUE)
#
# Check for redundant terms
#
  if(proposaltype=="ConstantEdges" && "edges" %in% model.initial$coef.names){
    stop("The model contains an 'edges' term and a proposal that holds the edges fixed. One of them is redundant. Please restate the model.")
  }
  BD <- ergm.boundDeg(con$boundDeg, nnodes=network.size(nw))
  Clist.initial <- ergm.Cprepare(nw, model.initial)
  Clist.initial$meanstats=meanstats
  theta0copy <- theta0
  initialfit <- ergm.initialfit(theta0copy, MLestimate, Clist.initial,
                                model.initial, verbose=verbose, ...)
  if (MLestimate && 
      (!ergm.independencemodel(model.initial) || !is.null(meanstats))) {
    theta0 <- initialfit$coef
    names(theta0) <- model.initial$coef.names
    theta0[is.na(theta0)] <- 0
  } else { # Just return initial (non-MLE) fit and exit.
    initialfit$network <- nw
    initialfit$newnetwork <- nw
    initialfit$formula <- formula
    return(initialfit)
  } 
  if (verbose) cat("Fitting ERGM.\n")
  model <- ergm.getmodel(formula, nw, drop=con$drop, expanded=TRUE)
  theta0 <- ergm.revisetheta0(model, theta0)
  # revise theta0 to reflect additional parameters

  Clist <- ergm.Cprepare(nw, model)
  Clist$meanstats <- meanstats
  Clist$obs <- summary(model$formula)

  if (verbose) cat("ergm.mainfitloop\n")
  styles <- c("Newton-Raphson","Robbins-Monro")
  con$style <- styles[pmatch(con$style,styles,nomatch=1)]
  if(!is.null(dissolve)){
    model.dissolve <- ergm.getmodel.dissolve(dissolve, nw)
    v <- ergm.robmon.dyn(theta0, nw, model, model.dissolve,
                    Clist, BD, gamma, burnin, interval,
                    proposaltype, verbose, con)
  }else{
   if(con$style == "Robbins-Monro"){
    v <- ergm.robmon(theta0, nw, model, Clist, BD, burnin, interval,
                     proposaltype, verbose, con)
   }else{
    v <- ergm.mainfitloop(theta0, nw,
                          model, Clist,
                          BD, initialfit, burnin, MCMCsamplesize,
                          interval, maxit, proposaltype, con$proposalpackage,
                          compress=con$compress, verbose=verbose, 
                          mcmc.precision=con$mcmc.precision,
                          nr.maxit=con$nr.maxit, calc.mcmc.se=con$calc.mcmc.se,
                          hessian=con$hessian, trustregion=con$trustregion,
                          steplength=con$steplength,
                          ...)
   }
  }

  if (!con$returnMCMCstats)
    v$sample <- NULL
  options(warn=current.warn)
  v
}
