###############################################################################
# The <ergm> function fits ergms from a specified formula returning either
# MPLEs or approximate MLE's based on MCMC estimation.
#
#
# --PARAMETERS--
#   formula       :  a formula of the form 'nw ~ model term(s)'
#   theta0        :  a vector of starting values for estimation, or optionally
#                    if these are to be estimated, the string "MPLE";
#                    default="MPLE"
#   MPLEonly      :  whether MPL estimation should be used (T or F); this is
#                    ignored if 'MLestimate' is set; default=FALSE
#   MLestimate    :  this is a logical indicating whether ML
#                    estimation should be used (TRUE or FALSE)
#                    default='!MPLEonly'
#   seed          :  an integer starting value for the random number generator;
#                    default=NULL
#   burnin        :  the number of proposals to ignore before MCMC sampling
#                    begins; default=10,000
#   MCMCsamplesize:  the number of network statistics to sample;
#                    default=10,000
#   interval      :  the number of proposals between sampled statistics;
#                    default=100
#   maxit         :  the number of MCMC parameter updates to the value
#                    maximizing the MCMC likelihood; default=3
#   constraints   :  a one-sided formula of the constraint terms; options are
#                         bd        degrees        nodegrees
#                         edges     degreedist     indegreedist
#                         observed  outdegreedist
#                    default="~ ."
#   meanstats     :  a vector of the mean value parameters;
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
#                   mcpar        : the following vector taken from MCMCparams:
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
#    $* @% &  +  burnin          :  the 'burnin' value inputted to <ergm>
#    $* @% &  +  interval        :  the 'interval' value inputted to <ergm>
#    $* @% &  +  theta.original  :  the theta values at the start of the MCMC sampling
#    $* @        mplefit         :  the MPLE fit as a glm object, and returned by
#                                   <ergm.mple>
#    $*!@% &#~+  null.deviance   :  the deviance of the null model
#    $*!@% &#~+  mle.lik         :  the approximate log-likelihood for the MLE
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
#         &      rm.coef         :  the robmon coefficients used as 'theta0' in the final
#                                   estimation
#          #~    aic             :  the Akaike information criterion value
#      !   #~   loglikelihoodratio: the log-likelihood corresponding to
#                                   'coef'
#
#####################################################################################    

ergm <- function(formula, response=NULL, theta0="MPLE",
                 MPLEonly=FALSE, MLestimate=!MPLEonly, seed=NULL,
                 burnin=10000, MCMCsamplesize=10000, interval=100,
                 maxit=3,
                 reference="Bernoulli",obs=NA,
                 constraints=~.,
                 meanstats=NULL,
                 control=control.ergm(),
                 eval.loglik=FALSE,
                 verbose=FALSE,...) {
  current.warn <- options()$warn
  options(warn=0)
  if(!is.null(seed))  set.seed(as.integer(seed))
  if (verbose) cat("Evaluating network in model\n")

  nw <- ergm.getnetwork(formula)
  proposalclass <- "c"

  # Construct the constraint for the observation process.
  # There may be a better way to specify this in the future.
  
  MHproposal.obs<-constraints
  
  MHproposal.obs<-switch(tolower(obs),
                         detrank = ergm.update.formula(MHproposal.obs,~.+ranks),
                         MHproposal.obs)
  
  if(network.naedgecount(nw)) MHproposal.obs<-ergm.update.formula(MHproposal.obs,~.+observed)
  
  if(constraints==MHproposal.obs) MHproposal.obs<-NULL

  if(!is.null(meanstats)){
   if(!(!is.null(control$SAN.burnin) && is.na(control$SAN.burnin))){
    netsumm<-summary(formula,response=response)
    if(length(netsumm)!=length(meanstats))
      stop("Incorrect length of the meanstats vector: should be ", length(netsumm), " but is ",length(meanstats),".")
    
    control$drop <- FALSE

    if(verbose) cat("Constructing an approximate response network.\n")
    ## If meanstats are given, overwrite the given network and formula
    ## with SAN-ed network and formula.
    srun <- 0
    obs <- meanstats-meanstats
    while(sum((obs-meanstats)^2) > 5){
     nw<-san(formula, meanstats=meanstats,
             theta0=if(is.numeric(theta0)) theta0,
             response=response,
             reference=reference,
             constraints=constraints,
             verbose=verbose,
             burnin=
             if(is.null(control$SAN.burnin)) burnin
             else control$SAN.burnin,
             interval=interval)
     formula<-ergm.update.formula(formula,nw~.)
     obs <- summary(formula,response=response, basis=nw)
     srun <- srun + 1
     if(verbose){
      cat(paste("Finished SAN run",srun,"\n"))
     }
    if(verbose){
      cat("SAN summary statistics:\n")
      print(obs)
      cat("Meanstats Goal:\n")
      print(meanstats)
      cat("Difference: SAN meanstats - Goal meanstats =\n")
      print(round(obs-meanstats,0))
    }
    }
   }
  }
  if(control$nsubphases=="maxit") control$nsubphases<-maxit
  
  if (verbose) cat("Initializing model.\n")

  if(control$drop){
   model.initial <- ergm.getmodel(formula, nw, response=response, drop=FALSE, initialfit=TRUE)
   model.initial.drop <- ergm.getmodel(formula, nw, response=response, drop=TRUE, initialfit=TRUE)
   namesmatch <- match(model.initial$coef.names, model.initial.drop$coef.names)
   droppedterms <- rep(FALSE, length=length(model.initial$etamap$offsettheta))
   droppedterms[is.na(namesmatch)] <- TRUE
   model.initial$etamap$offsettheta[is.na(namesmatch)] <- TRUE
   model.initial$etamap$offsetmap[model.initial$etamap$canonical[droppedterms & (model.initial$etamap$canonical>0)]] <- TRUE
  }else{
   model.initial <- ergm.getmodel(formula, nw, response=response, drop=control$drop, initialfit=TRUE)
   droppedterms <- rep(FALSE, length=length(model.initial$etamap$offsettheta))
  }
  if (verbose) { cat("Initializing Metropolis-Hastings proposal.\n") }

  MHproposal <- MHproposal(constraints, weights=control$prop.weights, control$prop.args, nw, model.initial,class=proposalclass,reference=reference,response=response)
  # Note:  MHproposal function in CRAN version does not use the "class" argument for now
  if(!is.null(MHproposal.obs)) MHproposal.obs <- MHproposal(MHproposal.obs, weights=control$prop.weights, control$prop.args, nw, model.initial, class=proposalclass, reference=reference, response=response)

  conddeg <- switch(MHproposal$name %in% c("CondDegree","CondDegreeSimpleTetrad","BipartiteCondDegHexadToggles","BipartiteCondDegTetradToggles"),control$drop,NULL)
  MCMCparams=c(control,
   list(samplesize=MCMCsamplesize, burnin=burnin, interval=interval,
        maxit=maxit,
	mcmc.precision=control$mcmc.precision))


  if (verbose) { cat("Fitting initial model.\n") }
  if(reference!="Bernoulli" && theta0=="MPLE") stop("MPLE initial values are not implemented for weithed network ERGMs. Please specify theta0 manually.")
  theta0copy <- theta0
  initialfit <- ergm.initialfit(theta0=theta0copy, MLestimate=MLestimate, 
                                formula=formula, nw=nw, meanstats=meanstats,
                                m=model.initial,
                                MPLEtype=control$MPLEtype, 
                                initial.loglik=control$initial.loglik,
                                conddeg=conddeg, MCMCparams=MCMCparams, MHproposal=MHproposal,
                                force.MPLE=(reference=="Bernoulli" && ergm.independencemodel(model.initial)
                                            && !control$force.mcmc
                                            && constraints==(~.)),
                                verbose=verbose, 
                                compressflag = control$compress, 
                                maxNumDyadTypes=control$maxNumDyadTypes,
                                ...)
  MCMCflag <- ((MLestimate && (!ergm.independencemodel(model.initial)
                               || !is.null(meanstats)
                               || constraints!=(~.)))
                || control$force.mcmc || reference!="Bernoulli")
  if (MCMCflag) {
    theta0 <- initialfit$coef
    names(theta0) <- model.initial$coef.names
    theta0[is.na(theta0)] <- 0
  } else { # Just return initial (non-MLE) fit and exit.
    initialfit$offset <- model.initial$etamap$offsettheta
    initialfit$drop <- droppedterms
    initialfit$network <- nw
    initialfit$reference <- reference
    initialfit$newnetwork <- nw
    initialfit$formula <- formula
    initialfit$constraints <- constraints
    initialfit$prop.args <- control$prop.args
    initialfit$prop.weights <- control$prop.weights
    initialfit<-logLik.ergm(initialfit, nsteps=loglik.nsteps, add=TRUE)
    return(initialfit)
  } 
  if(control$drop){
   model <- ergm.getmodel(formula, nw, response=response, drop=FALSE, expanded=TRUE, silent=TRUE)
   # revise theta0 to reflect additional parameters
   theta0 <- ergm.revisetheta0(model, theta0)
   model.drop <- ergm.getmodel(formula, nw, response=response, drop=TRUE, expanded=TRUE, silent=TRUE)
#            silent="MPLE" %in% theta0copy)
   eta0 <- ergm.eta(theta0, model$etamap)
   namesdrop <- model$coef.names[is.na(match(model$coef.names, model.drop$coef.names))]
   names(model$etamap$offsettheta) <- names(theta0)
   droppedterms <- rep(FALSE, length=length(model$etamap$offsettheta))
   droppedterms[is.na(namesmatch)] <- TRUE
   theta0[droppedterms] <- -Inf
   model$etamap$offsettheta[names(model$etamap$offsettheta) %in% namesdrop] <- TRUE
   model$etamap$offsetmap[model$etamap$canonical[droppedterms & (model$etamap$canonical>0)]] <- TRUE
  }else{
   model <- ergm.getmodel(formula, nw, response=response, drop=control$drop, expanded=TRUE)
   # revise theta0 to reflect additional parameters
   theta0 <- ergm.revisetheta0(model, theta0)
   names(model$etamap$offsettheta) <- names(theta0)
   droppedterms <- rep(FALSE, length=length(model$etamap$offsettheta))
  }

  Clist <- ergm.Cprepare(nw, model, response=response)
  Clist$obs <- summary(model$formula, drop=FALSE, response=response)

  Clist$meanstats <- Clist$obs
  if(!is.null(meanstats)){
   if (is.null(names(meanstats))){
    if(length(meanstats) == length(Clist$obs)){
     names(meanstats) <- names(Clist$obs)
     Clist$meanstats <- meanstats
    }else{
     namesmatch <- names(summary(model$formula, drop=FALSE, response=response))
     if(length(meanstats) == length(namesmatch)){
       namesmatch <- match(names(meanstats), namesmatch)
       Clist$meanstats <- meanstats[namesmatch]
     }
    }
   }else{
    namesmatch <- match(names(Clist$obs), names(meanstats))
    Clist$meanstats[!is.na(namesmatch)] <- meanstats[namesmatch[!is.na(namesmatch)]]
   }
  }

  MCMCparams <- c(control,
                  list(samplesize=MCMCsamplesize, burnin=burnin,
                       interval=interval,
                       maxit=maxit,
                       mcmc.precision=control$mcmc.precision))


  if (verbose) cat("Fitting ERGM.\n")
  v <- switch(control$style,
    "Robbins-Monro" = ergm.robmon(theta0, nw, model, Clist, burnin, interval,
                      MHproposal(constraints,weights=control$prop.weights, control$prop.args, nw, model, response=response), verbose, control),
    "Stochastic-Approximation" = ergm.stocapprox(theta0, nw, model, 
                                 Clist, 
                                 MCMCparams=MCMCparams, MHproposal=MHproposal,
                                 verbose),
    "Stepping" = ergm.stepping(theta0, nw, model, Clist, initialfit, 
				#nstats=nstats, 
				#approx=lognormapprox, filename.prefix=NULL, 
				#control=control.ergm(nsim1=100, nsim2=1000, gridsize=100),  # simulation parameters
				#plots=FALSE,  # currently useless, but plots can be reimplemented
				MCMCparams=MCMCparams, 
				MHproposal=MHproposal, MHproposal.obs=MHproposal.obs, 
				verbose=verbose,...),
     ergm.mainfitloop(theta0, nw,
                          model, Clist, 
                          initialfit,
                          MCMCparams=MCMCparams, MHproposal=MHproposal,
                          MHproposal.obs=MHproposal.obs,
                          verbose=verbose,
                      response=response,
                          ...)
              )

  if(!is.null(MCMCparams$check.degeneracy) && MCMCparams$check.degeneracy && (is.null(v$theta1$independent) || !all(v$theta1$independent))){
    if(verbose) {
      cat("Checking for degeneracy.\n")
    }
    degeneracy <- ergm.degeneracy(v, test.only=TRUE)
  } else {
    degeneracy <- list(degeneracy.value=NULL, degeneracy.type=NULL)
  }
  v$degeneracy.value <- degeneracy$degeneracy.value
  v$degeneracy.type <- degeneracy$degeneracy.type
  if(exists("formula.passed")){
    v$formula <- formula.passed
  }else{
    v$formula <- formula
  }
  v$constraints <- constraints
  v$prop.args <- control$prop.args
  v$prop.weights <- control$prop.weights

  v$response<-response
  v$reference<-reference

  v$offset <- model$etamap$offsettheta
  v$drop <- droppedterms
  v$etamap <- model$etamap
  if (!control$returnMCMCstats)
    v$sample <- NULL
  options(warn=current.warn)
  if (MCMCflag) {
    cat("\nThis model was fit using MCMC.  To examine model diagnostics", 
        "and check for degeneracy, use the mcmc.diagnostics() function.\n")
  }
  if(eval.loglik)
    v<-logLik.ergm(v, nsteps=control$loglik.nsteps, add=TRUE)
  v
}
