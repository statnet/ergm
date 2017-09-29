#  File R/ergm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
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
                 obs.constraints=~-observed,
                 offset.coef=NULL,
                 target.stats=NULL,
                 eval.loglik=TRUE,
                 estimate=c("MLE", "MPLE", "CD"),
                 control=control.ergm(),
                 verbose=FALSE,...) {
  check.control.class("ergm", "ergm")
  control.toplevel(control,...)
  
  estimate <- match.arg(estimate)

  if(estimate=="CD"){
    control$init.method <- "CD"
    eval.loglik <- FALSE
  }

  if(estimate=="MPLE"){
    control$init.method <- "MPLE"
  }
  
  if(!is.null(control$seed))  set.seed(as.integer(control$seed))
  if (verbose) message("Evaluating network in model.")
  
  nw <- ergm.getnetwork(formula)
  proposalclass <- "c"
  
  # Handle the observation process constraints.
  tmp <- .handle.obs.constraints(nw, constraints, obs.constraints, target.stats)
  nw <- tmp$nw
  constraints.obs <- tmp$constraints.obs

  model <- ergm.getmodel(formula, nw, response=response)

  ## Construct approximate response network if target.stats are given.
  if(!is.null(target.stats)){
    nw.stats <- ergm.getglobalstats(nw, model, response=response)[!model$etamap$offsetmap]
    target.stats <- vector.namesmatch(target.stats, names(nw.stats))
    target.stats <- na.omit(target.stats)
    if(length(nw.stats)!=length(target.stats)){
      stop("Incorrect length of the target.stats vector: should be ", length(nw.stats), " but is ",length(target.stats),". Note that offset() terms should *not* get target statistics.")
    }
    
    # no need to pass the offset term's init to SAN
    san.control <- control$SAN.control
    san.control$coef <- san.control$coef[!model$etamap$offsettheta]
    
    if(verbose) message("Constructing an approximate response network.")
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
          message(paste("Finished SAN run",srun,""))
        }
        if(verbose){
          message("SAN summary statistics:")
          .message_print(nw.stats)
          message("Meanstats Goal:")
          .message_print(target.stats)
          message("Difference: SAN target.stats - Goal target.stats =")
          .message_print(round(nw.stats-target.stats,0))
        }
        if(sum((nw.stats-target.stats)^2) <= 5) break
      }
    }
    
    # From this point on, target.stats has NAs corresponding to the
    # offset terms.
    #
    # TODO: Only have target.stats contain non-offset terms'
    # statistics, and have the rest of the code handle it
    # intelligently.
    target.stats <- .align.target.stats.offset(model, target.stats)   

  } else {
    if (network.edgecount(nw) == 0) warning("Network is empty and no target stats are specified.")
  }
  
  if (verbose) message("Initializing Metropolis-Hastings proposal(s):",appendLF=FALSE) 
  
  MHproposal <- MHproposal(constraints, weights=control$MCMC.prop.weights, control$MCMC.prop.args, nw, class=proposalclass,reference=reference,response=response)
  if (verbose) message(" Unconstrained ",MHproposal$pkgname,":MH_",MHproposal$name, " ", appendLF=FALSE)

  if(!is.null(MHproposal$auxiliaries)){
    if(verbose) message("(requests auxiliaries: reinitializing model)")
    model <- ergm.getmodel(formula, nw, response=response, extra.aux=list(MHproposal$auxiliaries))
  }else message()
    
  if(!is.null(constraints.obs)){
    MHproposal.obs <- MHproposal(constraints.obs, weights=control$obs.MCMC.prop.weights, control$obs.MCMC.prop.args, nw, class=proposalclass, reference=reference, response=response)
    if (verbose) message("Constrained ",MHproposal.obs$pkgname,":MH_",MHproposal.obs$name, " ", appendLF=FALSE)

    if(!is.null(MHproposal.obs$auxiliaries)){
      if(verbose) message("(requests auxiliaries: reinitializing model)")
      model$obs.model <- ergm.getmodel(formula, nw, response=response, extra.aux=list(MHproposal.obs$auxiliaries))
    }else message()
  }else MHproposal.obs <- NULL

  
  
  
  # conddeg MPLE has been superceded, but let the user know:
  if(!is.directed(nw) && ("degrees" %in% names(MHproposal$arguments$constraints) ||
                                           all(c("b1degrees","b2degrees") %in% names(MHproposal$arguments$constraints)))) message("Note that degree-conditional MPLE has been removed in version 4.0, having been superceded by Contrastive Divergence.")  
  
  if (verbose) message("Initializing model.")
  
  # Construct the initial model.
  
  # The following kludge knocks out MPLE if the sample space
  # constraints are not dyad-independent. For example, ~observed
  # constraint is dyad-independent, while ~edges is not.
  #
  # TODO: Create a flexible and general framework to manage methods
  # for obtaining initial values.
  init.candidates <- ergm.init.methods(MHproposal$reference$name)
  if("MPLE" %in% init.candidates && !is.dyad.independent(MHproposal$arguments$constraints,
                                                         MHproposal.obs$arguments$constraints)){
    init.candidates <- init.candidates[init.candidates!="MPLE"]
    if(verbose) message("MPLE cannot be used for this constraint structure.")
  }

  control$init.method <- match.arg(control$init.method, init.candidates)
  if(verbose) message(paste0("Using initial method '",control$init.method,"'."))

  if(length(model$etamap$curved)){ # Curved model: use ergm.logitreg() rather than glm().
      control$MPLE.type <- "logitreg"
  }
  
  # If some control$init is specified...
  if(!is.null(control$init)){
    # Check length of control$init.
    if (length(control$init)!=length(model$etamap$offsettheta)) {
      if(verbose){
        message("control$init =")
        .message_print(control$init)
        message("number of statistics is ",length(model$coef.names), "")
      }
      stop(paste("Invalid starting parameter vector control$init:",
                 "wrong number of parameters.",
                 "If you are passing output from another ergm run as control$init,",
                 "in a model with curved terms, see help(enformulate.curved)."))
    }
  }else control$init <- rep(NA, length(model$etamap$offsettheta)) # Set the default value of control$init.
  
  if(!is.null(offset.coef)){
      # TODO: Names matching here?
      if(length(control$init[model$etamap$offsettheta])!=length(offset.coef))
          stop("Invalid offset parameter vector offset.coef: ",
               "wrong number of parameters: expected ",
               length(control$init[model$etamap$offsettheta]),
               " got ",length(offset.coef),".")
      control$init[model$etamap$offsettheta]<-offset.coef
  }
  
  # Make sure any offset elements are given in control$init.
  if(any(is.na(control$init) & model$etamap$offsettheta)) stop("The model contains offset terms whose parameter values have not been specified:", paste.and(model$coef.names[is.na(control$init)|model$offsettheta]), ".", sep="")
  
  # Check if any terms are constrained to a constant and issue a warning.
  constrcheck <- ergm.checkconstraints.model(model, MHproposal, control$init)
  model <- constrcheck$model; control$init <- constrcheck$init
  
  # Check if any terms are at their extremes and handle them depending on control$drop.
  extremecheck <- ergm.checkextreme.model(model=model, nw=nw, init=control$init, response=response, target.stats=target.stats, drop=control$drop)
  model <- extremecheck$model; control$init <- extremecheck$init
  
  # MPLE is not supported for valued ERGMs.
  if(estimate=="MPLE"){
    if(!is.null(response)) stop("Maximum Pseudo-Likelihood (MPLE) estimation for valued ERGM terms is not implemented at this time. You may want to use Contrastive Divergence by passing estimate='CD' instead.")
    if(!is.dyad.independent(MHproposal$arguments$constraints,
                            MHproposal.obs$arguments$constraints))
      stop("Maximum Pseudo-Likelihood (MPLE) estimation for ERGMs with dyad-dependent constraints is only implemented for certain degree constraints at this time. You may want to use Contrastive Divergence by passing estimate='CD' instead.")
  }
  
  if (verbose) { message("Fitting initial model.") }
  
  MPLE.is.MLE <- (MHproposal$reference$name=="Bernoulli"
                  && is.dyad.independent(model)
                  && !is.curved(formula, response=response)
                  && !control$force.main
                  && is.dyad.independent(MHproposal$arguments$constraints,
                                         MHproposal.obs$arguments$constraints))
  
  # If all other criteria for MPLE=MLE are met, _and_ SAN network matches target.stats directly, we can get away with MPLE.
  if (!is.null(target.stats) && !isTRUE(all.equal(target.stats,nw.stats))) message("Unable to match target stats. Using MCMLE estimation.")
  MCMCflag <- (estimate=="MLE" && (!MPLE.is.MLE
                                   || (!is.null(target.stats) && !isTRUE(all.equal(target.stats,nw.stats)))
  )
  || control$force.main)
  
  # Short-circuit the optimization if all terms are either offsets or dropped.
  if(all(model$etamap$offsettheta)){
    # Note that this cannot be overridden with control$force.main.
    message("All terms are either offsets or extreme values. No optimization is performed.")
    return(structure(list(coef=control$init,
                          iterations=0,
                          loglikelihood=NA,
                          mle.lik=NULL,
                          gradient=rep(NA,length=length(control$init)),
                          failure=TRUE,
                          offset=model$etamap$offsettheta,
                          drop=if(control$drop) extremecheck$extremeval.theta,
                          estimable=constrcheck$estimable,
                          network=nw,
                          reference=reference,
                          response=response,
                          newnetwork=nw,
                          formula=formula,
                          constrained=MHproposal$arguments$constraints,
                          constrained.obs=MHproposal.obs$arguments$constraints,
                          constraints=constraints,
                          target.stats=target.stats,
                          target.esteq=if(!is.null(target.stats)){
                            tmp <- .ergm.esteq(initialfit$coef, model, rbind(target.stats))
                            structure(c(tmp), names=colnames(tmp))
                          },
                          estimate=estimate,
                          control=control
    ),
    class="ergm"))
    
  }
  
  model$nw.stats <- ergm.getglobalstats(nw, model, response=response)
  model$target.stats <- target.stats
  
  if(control$init.method=="CD") if(is.null(names(control$init)))
      names(control$init) <- coef.names.model(model, FALSE)
  
  initialfit <- ergm.initialfit(init=control$init, initial.is.final=!MCMCflag,
                                formula=formula, nw=nw, reference=reference, 
                                m=model, method=control$init.method,
                                MPLEtype=control$MPLE.type, 
                                control=control,
                                MHproposal=MHproposal,
                                MHproposal.obs=MHproposal.obs,
                                verbose=verbose, response=response,
                                maxNumDyadTypes=control$MPLE.max.dyad.types,
                                ...)
  
  if (!MCMCflag){ # Just return initial (non-MLE) fit and exit.
    initialfit$offset <- model$etamap$offsettheta
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
    initialfit$obs.constraints <- obs.constraints 
    initialfit$target.stats <- suppressWarnings(na.omit(model$target.stats))
    initialfit$nw.stats <- model$nw.stats
      initialfit$etamap <- model$etamap
    initialfit$target.esteq <- suppressWarnings(na.omit(if(!is.null(model$target.stats)){
      tmp <- .ergm.esteq(initialfit$coef, model, rbind(model$target.stats))
      structure(c(tmp), names=colnames(tmp))
    }))
    initialfit$estimate <- estimate
    
    initialfit$control<-control
    
    if(eval.loglik) initialfit$null.lik <- logLikNull.ergm(initialfit, verbose=verbose)
    if(any(!model$etamap$offsettheta) && eval.loglik){
      message("Evaluating log-likelihood at the estimate. ",appendLF=FALSE)
      initialfit<-logLik.ergm(initialfit, add=TRUE, control=control$loglik.control, verbose=verbose)
      message("")
    }
    return(initialfit)
  }
  
  # Otherwise, set up the main phase of estimation:
  
  parallel.toplevel <- NULL     # top level reminder to stop cluster
  if (inherits(control$parallel,"cluster")) {
    clus <- ergm.getCluster(control, verbose)
  } else if(is.numeric(control$parallel) && control$parallel!=0){
    clus <- ergm.getCluster(control, verbose)
    ergm.cluster.started(FALSE)
    parallel.toplevel <- control$parallel
    control$parallel <- clus
  } else {
    clus <- NULL
    ergm.cluster.started(FALSE)
    if (!is.numeric(control$parallel))
      warning("Unrecognized value passed to parallel control parameter.")
  }
  
  # Revise the initial value, if necessary:
  init <- initialfit$coef
  init[is.na(init)] <- 0
  names(init) <- coef.names.model(model, FALSE)
  
  if (verbose) message("Fitting ERGM.")
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
  
  # done with main fit
  
  if(!is.null(control$MCMLE.check.degeneracy) && control$MCMLE.check.degeneracy && (is.null(mainfit$theta1$independent) || !all(mainfit$theta1$independent))){
    if(verbose) {
      message("Checking for degeneracy.")
    }
    degeneracy <- ergm.degeneracy(mainfit, test.only=TRUE)
  } else {
    degeneracy <- list(degeneracy.value=NULL, degeneracy.type=NULL)
  }
  mainfit$degeneracy.value <- degeneracy$degeneracy.value
  mainfit$degeneracy.type <- degeneracy$degeneracy.type
  
  mainfit$formula <- formula
  mainfit$target.stats <- suppressWarnings(na.omit(model$target.stats))
  mainfit$nw.stats <- model$nw.stats
  mainfit$target.esteq <- suppressWarnings(na.omit(if(!is.null(model$target.stats)){
    tmp <- .ergm.esteq(mainfit$coef, model, rbind(model$target.stats))
    structure(c(tmp), names=colnames(tmp))
  }))
  
  mainfit$constrained <- MHproposal$arguments$constraints
  mainfit$constrained.obs <- MHproposal.obs$arguments$constraints
  mainfit$constraints <- constraints
  mainfit$obs.constraints <- obs.constraints
  
  # unless the main fitting algorithm passes back a modified control
  if (is.null(mainfit$control)) mainfit$control<-control
  
  mainfit$response<-response
  mainfit$reference<-reference
  mainfit$estimate <- estimate
  
  mainfit$offset <- model$etamap$offsettheta
  mainfit$drop <- if(control$drop) extremecheck$extremeval.theta
  mainfit$estimable <- constrcheck$estimable
  mainfit$etamap <- model$etamap
  
  mainfit$null.lik<-logLikNull.ergm(mainfit, verbose=verbose)
  
  if (!control$MCMC.return.stats)
    mainfit$sample <- NULL
  
  if(eval.loglik){
    message("Evaluating log-likelihood at the estimate. ", appendLF=FALSE)
    mainfit<-logLik.ergm(mainfit, add=TRUE, control=control$loglik.control, verbose=verbose)
  }
  
  # done with parallel cluster
  if (!is.null(parallel.toplevel)) {
    mainfit$control$parallel <- parallel.toplevel
    ergm.cluster.started(TRUE)
  }
  ergm.stopCluster(clus)
  
  if (MCMCflag) {
    message("This model was fit using MCMC.  To examine model diagnostics ", 
        "and check for degeneracy, use the mcmc.diagnostics() function.")
  }
  
  mainfit
}
