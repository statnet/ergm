#  File R/ergm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################
###############################################################################
# The <ergm> function fits ergms from a specified formula returning either
# MPLEs or approximate MLE's based on MCMC estimation.
#
#
# --PARAMETERS--
#   formula       :  a formula of the form 'nw ~ model term(s)'
#   init        :  a vector of starting values for estimation or offset values, or optionally
#                    if these are to be estimated, NULL (the default);
#   constraints   :  a one-sided formula of the constraint terms; options are
#                         bd        degrees        nodegrees
#                         edges     degreedist     idegreedist
#                         observed  odegreedist
#                    default="~ ."
#   target.stats     :  a vector of the mean value parameters;
#                    default=the observed statistic from the 'nw' in formula
#   control       :  a list of control parameters returned from <control.ergm>;
#                    default=control.ergm()
#   verbose       :  whether ergm should be verbose (T or F); default=FALSE
#
#
# --RETURNED--
#   because a stergm object is the return type of several functions, and
#   because this is a rather lengthy list, and because the returned items
#   of this function borrow from the other stergm.* functions, this list
#   provides the returned items for all funtions returning a stergm.
#   The symbol preceding each component indicates which function returns it,
#   but remember that, <stergm> will additionally return the items from
#   one of the other stergm functions as well:                               
#   because an ergm object is the return type of several functions, and
#   because this is a rather lengthy list, this list represents the return
#   type for all funtions returning an ergm. The symbol preceding each
#   component indicates which function returns it:
#       <ergm>             = $
#       <ergm.mainfitloop> = *
#       <ergm.mple>        = !
#       <ergm.stepping>    = @
#       <ergm.stocapprox>  = %
#       <ergm.estimate>    = ^
#       <ergm.robmon>      = &
#       <ergm.mapl>        = #
#       <ergm.maple>       = ~
#
#   the components include:
#
#    $*!@%^&#~+  coef            :  the vector of estimated model coefficients
#    $* @%^&  +  sample          :  the row-binded matrix of network statistics from
#                                   each sample; 'sample' will also have 2 attributes:
#                   mcpar        : the following vector taken from control:
#                                        c(burnin+1, endrun, interval)
#                                  where 'endrun' is defined as
#                                     burnin+interval*(samplesize-1)
#                   class        : "mcmc"
#    $*!@%^&#~+  iterations      :  the number of Newton-Raphson iterations required
#                                   before convergence
#    $*!@%^&#~+  MCMCtheta       :  the vector of natural parameters used to produce
#                                   the MCMC sample
#    $* @%^&  +  loglikelihood   :  the estimated change in log-likelihood in the last
#                                   iteration
#    $*!@%^&#~+  gradient        :  the value of the gradient of the approximated log-
#                                   likelihood function at the maximizing value
#    $*!@%^&#~+  covar           :  the approximated covariance matrix for the MLE
#    $*!@%^&#~+  samplesize      :  the size of the MCMC sample
#    $*!@%^&#~+  failure         :  whether estimation failed (T or F)
#    $*!@%^&#~+  mc.se           :  the standard error estimates
#    $* @% &# +  newnetwork      :  the final network sampled; in the ergm returned from
#                                   <ergm.robmom>, this='network'
#    $* @% &# +  network         :  the 'nw' inputted to <ergm> via the 'formula'
#    $* @% &  +  theta.original  :  the theta values at the start of the MCMC sampling
#    $* @        mplefit         :  the MPLE fit as a glm object, and returned by
#                                   <ergm.mple>
#    $*!@% &#~+  mle.lik         :  the approximate log-likelihood for the MLE, if computed
#    $* @        etamap          :  the set of function mapping theta -> eta;
#                                   see <etamap>? for the components of this list
#    $           degeneracy.value:  the degeneracy value assigned by <ergm.degeneracy>
#    $           degeneracy.type :  a vector of length 2, as returned by
#                                   <ergm.compute.degeneracy> (found in the
#                                   <ergm.degeracy> file)
#    $      #    formula         :  the 'formula' value inputted to <ergm>
#    $      #    constraints     :  the 'constraints' value inputted to <ergm>
#           #    prop.args       :  the list of arguments that were passed onto the
#                                   <InitMHP> routines
#    $      #    prop.weights    :  the MCMC proposal weights inputted to <ergm> via
#                                  'control'
#    $      #    offset          :  a vector of whether each model parameter was set at
#                                  a fixed value (not estimated)
#    $      #    drop            :  list of dropped terms
#     * @%^&  +  sample.obs      :  the matrix of sample network statistics for observed
#                                   data
#     * @        parallel        :  the number of additional threads used when sampling
#      !    #~   glm             :  the fit established by MPL estimation and returned
#                                   by <ergm.logitreg>, <ergm.pen.glm> or <glm>
#                                   depending on the 'MPLEtype';
#      !    #~   glm.null        :  the null fit established by MPL estimation and
#                                   returned by <ergm.logitreg>, <ergm.pen.glm> or <glm>
#                                   depending on the 'MPLEtype';
#      !   #~    theta1          :  the vector of ??
#         &      rm.coef         :  the robmon coefficients used as 'init' in the final
#                                   estimation
#      !   #~   loglikelihoodratio: the log-likelihood corresponding to
#                                   'coef'
#
#####################################################################################    

ergm <- function(formula, response=NULL,
                 reference=~Bernoulli,
                 constraints=~.,
                 offset.coef=NULL,
                 target.stats=NULL,
                 eval.loglik=TRUE,
                 estimate=c("MLE", "MPLE"),
                 control=control.ergm(),
                 verbose=FALSE,...) {
  check.control.class()

  estimate <- match.arg(estimate)
  # Backwards-compatibility:
  control<-control.ergm.toplevel(control,...)
  if(!is.null(list(...)$MPLEonly) && list(...)$MPLEonly){
    warning("Argument MPLEonly is deprecated. Use ``estimate=\"MPLE\"'' instead." )
    estimate <- "MPLE"
  }

  if(!is.null(control$seed))  set.seed(as.integer(control$seed))
  if (verbose) cat("Evaluating network in model\n")
  
  nw <- ergm.getnetwork(formula)
  proposalclass <- "c"


  # Missing data handling only needs to happen if the sufficient
  # statistics are not specified. If the sufficient statistics are
  # specified, the nw's dyad states are irrelevant.
  if(network.naedgecount(nw) && !is.null(target.stats)){
    warning("Target statistics specified in a network with missing dyads. Missingness will be overridden.")
    nw[as.matrix(is.na(nw),matrix.type="edgelist")] <- 0
  }
 
  MHproposal.obs <- if(network.naedgecount(nw)==0) NULL else append.rhs.formula(constraints, list(as.name("observed")), TRUE)

  ## Construct approximate response network if target.stats are given.
  
  if(!is.null(target.stats)){
    nw.stats<-summary(remove.offset.formula(formula),response=response)
    target.stats <- vector.namesmatch(target.stats, names(nw.stats))
    target.stats <- na.omit(target.stats)
    if(length(nw.stats)!=length(target.stats)){
      stop("Incorrect length of the target.stats vector: should be ", length(nw.stats), " but is ",length(target.stats),". Note that offset() terms should *not* get target statistics.")
    }
    
    # no need to pass the offset term's init to SAN
    offset.terms <- offset.info.formula(formula)$term
    san.control <- control$SAN.control
    san.control$coef <- san.control$coef[!offset.terms]
    
    if(verbose) cat("Constructing an approximate response network.\n")
    ## If target.stats are given, overwrite the given network and formula
    ## with SAN-ed network and formula.
    if(control$SAN.maxit > 0){
     for(srun in 1:control$SAN.maxit){
      nw<-san(remove.offset.formula(formula), target.stats=target.stats,
              response=response,
              reference=reference,
              constraints=constraints,
              control=san.control,
              verbose=verbose)
      formula<-ergm.update.formula(formula,nw~., from.new="nw")
      nw.stats <- summary(remove.offset.formula(formula),response=response)
      srun <- srun + 1
      if(verbose){
        cat(paste("Finished SAN run",srun,"\n"))
      }
      if(verbose){
        cat("SAN summary statistics:\n")
        print(nw.stats)
        cat("Meanstats Goal:\n")
        print(target.stats)
        cat("Difference: SAN target.stats - Goal target.stats =\n")
        print(round(nw.stats-target.stats,0))
      }
      if(sum((nw.stats-target.stats)^2) <= 5) break
     }
    }

    offinfo <- offset.info.formula(formula,response=response)
    tmp <- rep(NA, length(offinfo$eta))
    tmp[!offinfo$eta] <- target.stats
    names(tmp)[!offinfo$eta] <- names(target.stats)
    s <- summary(formula,response=response)[offinfo$eta]
    tmp[offinfo$eta] <- s
    names(tmp)[offinfo$eta] <- names(s)
    
    # From this point on, target.stats has NAs corresponding to the
    # offset terms.
    #
    # TODO: Only have target.stats contain non-offset terms'
    # statistics, and have the rest of the code handle it
    # intelligently.
    target.stats <- tmp
  } else {
    if (network.edgecount(nw) == 0) warning("Network is empty and no target stats are specified.")
  }
  
  if (verbose) cat("Initializing Metropolis-Hastings proposal(s):") 
  
  MHproposal <- MHproposal(constraints, weights=control$MCMC.prop.weights, control$MCMC.prop.args, nw, class=proposalclass,reference=reference,response=response)
  if (verbose) cat(" ",MHproposal$pkgname,":MH_",MHproposal$name,sep="")

  
  # Note:  MHproposal function in CRAN version does not use the "class" argument for now
  if(!is.null(MHproposal.obs)){
      MHproposal.obs <- MHproposal(MHproposal.obs, weights=control$MCMC.prop.weights, control$MCMC.prop.args, nw, class=proposalclass, reference=reference, response=response)
      if (verbose) cat(" ",MHproposal.obs$pkgname,":MH_",MHproposal.obs$name,sep="")
  }

  if(verbose) cat("\n")

  # conddeg MPLE only handles tetrad toggles, so it must be restricted:
  conddeg <- switch(!is.directed(nw) && ("degrees" %in% names(MHproposal$arguments$constraints) ||
                                         all(c("b1degrees","b2degrees") %in% names(MHproposal$arguments$constraints))),
                    control$drop,
                    NULL)
  
  if (verbose) cat("Initializing model.\n")
  
  # Construct the initial model.
  control$init.method <- match.arg(control$init.method, ergm.init.methods(MHproposal$reference$name))
  model.initial <- ergm.getmodel(formula, nw, response=response, initialfit=control$init.method=="MPLE")
   
  # If some control$init is specified...
  if(!is.null(control$init)){
    # Check length of control$init.
    if (length(control$init)!=length(model.initial$etamap$offsettheta)) {
      if(verbose) cat("control$init is", control$init, "\n", "number of statistics is",length(model.initial$coef.names), "\n")
      stop(paste("Invalid starting parameter vector control$init:",
                 "wrong number of parameters.",
                 "If you are passing output from another ergm run as control$init,",
                 "in a model with curved terms, see help(enformulate.curved)."))
    }
  }else control$init <- rep(NA, length(model.initial$etamap$offsettheta)) # Set the default value of control$init.

  if(!is.null(offset.coef)){
      # TODO: Names matching here?
      if(length(control$init[model.initial$etamap$offsettheta])!=length(offset.coef))
          stop("Invalid offset parameter vector offset.coef: ",
               "wrong number of parameters: expected ",
               length(control$init[model.initial$etamap$offsettheta]),
               " got ",length(offset.coef),".")
      control$init[model.initial$etamap$offsettheta]<-offset.coef
  }
  
  # Make sure any offset elements are given in control$init.
  if(any(is.na(control$init) & model.initial$etamap$offsettheta)) stop("The model contains offset terms whose parameter values have not been specified:", paste.and(model.initial$coef.names[is.na(control$init)|model.initial$offsettheta]), ".", sep="")

  # Check if any terms are constrained to a constant and issue a warning.
  constrcheck <- ergm.checkconstraints.model(model.initial, MHproposal, control$init)
  model.initial <- constrcheck$model; control$init <- constrcheck$init

  # Check if any terms are at their extremes and handle them depending on control$drop.
  extremecheck <- ergm.checkextreme.model(model=model.initial, nw=nw, init=control$init, response=response, target.stats=target.stats, drop=control$drop)
  model.initial <- extremecheck$model; control$init <- extremecheck$init

  
  # Construct the curved model, and check if it's different from the initial model. If so, we know that it's curved.
  model <- ergm.getmodel(formula, nw, response=response, expanded=TRUE, silent=TRUE)
  # MPLE is not supported for curved ERGMs.
  if(length(model$etamap$offsetmap)!=length(model.initial$etamap$offsetmap) && estimate=="MPLE") stop("Maximum Pseudo-Likelihood (MPLE) estimation for curved ERGMs is not implemented at this time. You may want to pass fixed=TRUE parameter in curved terms to specify the curved parameters as fixed.")
  
  if (verbose) { cat("Fitting initial model.\n") }

  MPLE.is.MLE <- (MHproposal$reference$name=="Bernoulli"
                  && is.dyad.independent(model.initial)
                  && !control$force.main
                  && is.dyad.independent(MHproposal$arguments$constraints,
                                         MHproposal.obs$arguments$constraints))

  # If all other criteria for MPLE=MLE are met, _and_ SAN network matches target.stats directly, we can get away with MPLE.
  MCMCflag <- (estimate=="MLE" && (!MPLE.is.MLE
                               || (!is.null(target.stats) && !isTRUE(all.equal(target.stats,nw.stats)))
                              )
               || control$force.main)

  # Short-circuit the optimization if all terms are either offsets or dropped.
  if(MCMCflag && all(model.initial$etamap$offsettheta)){
    # Note that this cannot be overridden with control$force.main.
    MCMCflag <- FALSE
    warning("All terms are either offsets or extreme values. Skipping MCMC.")
  }

  model.initial$nw.stats <- summary(model.initial$formula, response=response, initialfit=control$init.method=="MPLE")
  model.initial$target.stats <- if(!is.null(target.stats)) target.stats else model.initial$nw.stats
  
  initialfit <- ergm.initialfit(init=control$init, initial.is.final=!MCMCflag,
                                formula=formula, nw=nw, reference=reference, 
                                m=model.initial, method=control$init.method,
                                MPLEtype=control$MPLE.type, 
                                conddeg=conddeg, control=control,
                                MHproposal=MHproposal,
                                MHproposal.obs=MHproposal.obs,
                                verbose=verbose, response=response,
                                maxNumDyadTypes=control$MPLE.max.dyad.types,
                                ...)
  
  if (!MCMCflag){ # Just return initial (non-MLE) fit and exit.
    initialfit$offset <- model.initial$etamap$offsettheta
    initialfit$drop <- if(control$drop) extremecheck$extremeval.theta
    initialfit$estimable <- constrcheck$estimable
    initialfit$network <- nw
    initialfit$reference <- reference
    initialfit$response <- response
    initialfit$newnetwork <- nw
    initialfit$formula <- formula
    initialfit$constrained <- MHproposal$arguments$constraints
    initialfit$constrained.obs <- MHproposal.obs$arguments$constraints
    initialfit$constraints <- constraints
    initialfit$target.stats <- model.initial$target.stats
    initialfit$target.esteq <- if(!is.null(model.initial$target.stats)){
      tmp <- .ergm.esteq(initialfit$coef, model.initial, rbind(model.initial$target.stats))
      structure(c(tmp), names=colnames(tmp))
    }
    initialfit$estimate <- estimate

    initialfit$control<-control

    if(eval.loglik) initialfit$null.lik <- logLikNull.ergm(initialfit)
    if(any(!model.initial$etamap$offsettheta) && eval.loglik){
      if(verbose) cat("Evaluating log-likelihood at the estimate.\n")
      initialfit<-logLik.ergm(initialfit, add=TRUE, control=control$loglik.control)
    }
    return(initialfit)
  }

  # Otherwise, set up the main phase of estimation:
  # Revise the initial value, if necessary:
  init <- initialfit$coef
  init[is.na(init)] <- 0
  if(control$init.method=="MPLE"){ # Only MPLE requires these kludges.
      names(init) <- model.initial$coef.names
      # revise init to reflect additional parameters
      init <- ergm.reviseinit(model, init)
  }
  names(init) <- .coef.names.model(model, FALSE)

  # Check if any terms are constrained to a constant and issue a warning.
  constrcheck <- ergm.checkconstraints.model(model, MHproposal, init=init, silent=TRUE)
  model <- constrcheck$model; control$init <- constrcheck$init
  
  # Check if any terms are at their extremes and handle them depending on control$drop.
  extremecheck <- ergm.checkextreme.model(model=model, nw=nw, init=init, response=response, target.stats=target.stats, drop=control$drop, silent=TRUE)
  model <- extremecheck$model; init <- extremecheck$init

  model$nw.stats <- summary(model$formula, response=response)
  model$target.stats <- if(!is.null(target.stats)) target.stats else model$nw.stats

  if (verbose) cat("Fitting ERGM.\n")
  mainfit <- switch(control$main.method,
    "Robbins-Monro" = ergm.robmon(init, nw, model, 
      MHproposal=MHproposal, verbose=verbose, control=control),
    "Stochastic-Approximation" = ergm.stocapprox(init, nw, model, 
                                 control=control, MHproposal=MHproposal,
                                 verbose),
    "Stepping" = ergm.stepping(init, nw, model, initialfit, constraints,
				#nstats=nstats, 
				#approx=lognormapprox, filename.prefix=NULL, 
				#control=control.ergm(nsim1=100, nsim2=1000, gridsize=100),  # simulation parameters
				#plots=FALSE,  # currently useless, but plots can be reimplemented
				control=control, 
				MHproposal=MHproposal, MHproposal.obs=MHproposal.obs, 
				verbose=verbose,...),
    "MCMLE" = ergm.MCMLE(init, nw,
                          model, 
                          # no need to pass initialfit to MCMLE
                          initialfit=(initialfit<-NULL),
                          control=control, MHproposal=MHproposal,
                          MHproposal.obs=MHproposal.obs,
                          verbose=verbose,
                      response=response,
                          ...),


              stop("Method ", control$main.method, " is not implemented.")
              )
  
  initialfit <- NULL

  if(!is.null(control$MCMLE.check.degeneracy) && control$MCMLE.check.degeneracy && (is.null(mainfit$theta1$independent) || !all(mainfit$theta1$independent))){
    if(verbose) {
      cat("Checking for degeneracy.\n")
    }
    degeneracy <- ergm.degeneracy(mainfit, test.only=TRUE)
  } else {
    degeneracy <- list(degeneracy.value=NULL, degeneracy.type=NULL)
  }
  mainfit$degeneracy.value <- degeneracy$degeneracy.value
  mainfit$degeneracy.type <- degeneracy$degeneracy.type

  mainfit$formula <- formula
  mainfit$target.stats <- model$target.stats
  mainfit$target.esteq <- if(!is.null(model$target.stats)){
    tmp <- .ergm.esteq(mainfit$coef, model, rbind(model$target.stats))
    structure(c(tmp), names=colnames(tmp))
  }

  mainfit$constrained <- MHproposal$arguments$constraints
  mainfit$constrained.obs <- MHproposal.obs$arguments$constraints
  mainfit$constraints <- constraints

  # unless the main fitting algorithm passes back a modified control
  if (is.null(mainfit$control)) mainfit$control<-control

  mainfit$response<-response
  mainfit$reference<-reference
  mainfit$estimate <- estimate

  mainfit$offset <- model$etamap$offsettheta
  mainfit$drop <- if(control$drop) extremecheck$extremeval.theta
  mainfit$estimable <- constrcheck$estimable
  mainfit$etamap <- model$etamap

  mainfit$null.lik<-logLikNull.ergm(mainfit)
  
  if (!control$MCMC.return.stats)
    mainfit$sample <- NULL

  if(eval.loglik){
    if(verbose) cat("Evaluating log-likelihood at the estimate.\n")
    mainfit<-logLik.ergm(mainfit, add=TRUE, control=control$loglik.control)
  }

  if (MCMCflag) {
    cat("\nThis model was fit using MCMC.  To examine model diagnostics", 
        "and check for degeneracy, use the mcmc.diagnostics() function.\n")
  }

  mainfit
}
