#  File R/ergm.bridge.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################

## This is a helper function that constructs and returns the network
## object to be used and the model object.
ergm.bridge.preproc<-function(object, basis, response, ...){

  # If basis is not null, replace network in formula by basis.
  # In either case, let nw be network object from formula.
  if(is.null(nw <- basis)) {
    nw <- ergm.getnetwork(object)
  }

  nw <- ensure_network(nw)
  
  # New formula (no longer use 'object'):
  form <- nonsimp_update.formula(object, nw ~ ., from.new="nw")
  
  list(nw=nw, form=form, model=ergm_model(form, nw, response=response, ...))
}


#' Bridge sampling to evaluate ERGM log-likelihoods and log-likelihood ratios
#' 
#' \code{ergm.bridge.llr} uses bridge sampling with geometric spacing to
#' estimate the difference between the log-likelihoods of two parameter vectors
#' for an ERGM via repeated calls to \code{\link{simulate.formula.ergm}}.
#' 
#' 
#' 
#' @param object A model formula. See \code{\link{ergm}} for details.
#' @template response
#' @param constraints A one-sided formula specifying one or more
#'   constraints on the support of the distribution of the networks
#'   being simulated. See the documentation for a similar argument for
#'   \code{\link{ergm}} for more information.
#' @param from,to The initial and final parameter vectors.
#' @param basis An optional \code{\link[network]{network}} object to
#'   start the Markov chain.  If omitted, the default is the
#'   left-hand-side of the \code{object}.
#' @param verbose Logical: If TRUE, print detailed information.
#' @param \dots Further arguments to \code{ergm.bridge.llr} and
#'   \code{\link{simulate.formula.ergm}}.
#' @param llronly Logical: If TRUE, only the estiamted log-ratio will
#'   be returned by `ergm.bridge.llr`.
#' @param control Control arguments.  See
#'   \code{\link{control.ergm.bridge}} for details.
#' @param coef A vector of coefficients for the configuration of
#'   interest.
#' @param llkonly Whether only the estiamted log-likelihood should be
#'   returned by the `ergm.bridge.0.llk` and
#'   `ergm.bridge.dindstart.llk`.  (Defaults to TRUE.)
#' @return If `llronly=TRUE` or `llkonly=TRUE`, these functions return
#'   the scalar log-likelihood-ratio or the log-likelihood.
#'   Otherwise, they return a list with the following components:
#'   \item{llr}{The estimated log-ratio.}  \item{llrs}{The estimated
#'   log-ratios for each of the \code{nsteps} bridges.}  \item{path}{A
#'   numeric matrix with nsteps rows, with each row being the
#'   respective bridge's parameter configuration.}  \item{stats}{A
#'   numeric matrix with nsteps rows, with each row being the
#'   respective bridge's vector of simulated statistics.}
#'   \item{Dtheta.Du}{The gradient vector of the parameter values with
#'   respect to position of the bridge.}
#' @seealso \code{\link{simulate.formula.ergm}}
#' @references Hunter, D. R. and Handcock, M. S. (2006)
#'   \emph{Inference in curved exponential family models for
#'   networks}, Journal of Computational and Graphical Statistics.
#' @keywords model
#' @export
ergm.bridge.llr<-function(object, response=NULL, constraints=~., from, to, basis=NULL, verbose=FALSE, ..., llronly=FALSE, control=control.ergm.bridge()){
  check.control.class("ergm.bridge", "ergm.bridge.llr")
  control.toplevel(..., myname="ergm.bridge")

  if(!is.null(control$seed)) {set.seed(as.integer(control$seed))}
  
  ## Here, we need to get the model object to get the likelihood and gradient functions.
  tmp<-ergm.bridge.preproc(object,basis,response,term.options=control$term.options)
  nw<-tmp$nw; m<-tmp$model; form<-tmp$form; rm(tmp)


  ## Generate the path.
  path<-t(rbind(sapply(seq(from=0+1/2/(control$nsteps+1),to=1-1/2/(control$nsteps+1),length.out=control$nsteps),function(u) cbind(to*u + from*(1-u)))))

  stats<-matrix(NA,control$nsteps,m$etamap$etalength)
  
  if(network.naedgecount(nw)){
    constraints.obs<-nonsimp_update.formula(constraints,~.+observed)
    form.obs<-form
    stats.obs <- matrix(NA,control$nsteps,m$etamap$etalength)  
  }else stats.obs<-matrix(summary(form,response=response, term.options=control$term.options),control$nsteps,m$etamap$etalength,byrow=TRUE)  

  message("Using ", control$nsteps, " bridges: ", appendLF=FALSE)
  for(i in seq_len(control$nsteps)){
    theta<-path[i,]
    if(verbose==0) message(i," ",appendLF=FALSE)
    if(verbose>0) message("Running theta=[",paste(format(theta),collapse=","),"].")
    if(verbose>1) message("Burning in...")
    ## First burn-in has to be longer, but those thereafter should be shorter if the bridges are closer together.
    nw.state<-simulate(form, coef=theta, nsim=1, response=response, constraints=constraints, output="pending_update_network", verbose=max(verbose-1,0),
                       control=control.simulate.formula(MCMC.burnin=if(i==1) control$MCMC.burnin else ceiling(control$MCMC.burnin/sqrt(control$nsteps)),
                                                        term.options=control$term.options,
                         MCMC.interval=1,
                         MCMC.prop.args=control$MCMC.prop.args,
                         MCMC.prop.weights=control$MCMC.prop.weights,
                         MCMC.packagenames=control$MCMC.packagenames,
                         parallel=control$parallel,
                         parallel.type=control$parallel.type,
                         parallel.version.check=control$parallel.version.check
                                                        ), ...)
    nonsimp_update.formula(form,nw.state~., from.new="nw.state")
    stats[i,]<-colMeans(as.matrix(simulate(form, coef=theta, response=response, constraints=constraints, output="stats", verbose=max(verbose-1,0),
                              control=control.simulate.formula(MCMC.burnin=0,
                                                               MCMC.interval=control$MCMC.interval,
                                                               term.options=control$term.options,
                                                               MCMC.prop.args=control$MCMC.prop.args,
                                                               MCMC.prop.weights=control$MCMC.prop.weights,
                                                               MCMC.packagenames=control$MCMC.packagenames,
                                                               parallel=control$parallel,
                                                               parallel.type=control$parallel.type,
                                                               parallel.version.check=control$parallel.version.check
                                                               ),
                              nsim=ceiling(control$MCMC.samplesize/control$nsteps), ...)))
    
    if(network.naedgecount(nw)){
      nw.state.obs<-simulate(form.obs, coef=theta, nsim=1, response=response, constraints=constraints.obs, output="pending_update_network", verbose=max(verbose-1,0),
                             control=control.simulate.formula(MCMC.burnin=if(i==1) control$obs.MCMC.burnin else ceiling(control$obs.MCMC.burnin/sqrt(control$nsteps)),
                                                              term.options=control$term.options,
                               MCMC.interval=1,
                               MCMC.prop.args=control$MCMC.prop.args,
                               MCMC.prop.weights=control$MCMC.prop.weights,
                               MCMC.packagenames=control$MCMC.packagenames,
                               parallel=control$parallel,
                               parallel.type=control$parallel.type,
                               parallel.version.check=control$parallel.version.check), ...)
      nonsimp_update.formula(form.obs,nw.state.obs~., from.new="nw.state.obs")
      stats.obs[i,]<-colMeans(as.matrix(simulate(form.obs, coef=theta, response=response, constraints=constraints.obs, output="stats", verbose=max(verbose-1,0),
                                control=control.simulate.formula(MCMC.burnin=0,
                                  MCMC.interval=control$obs.MCMC.interval,
                                  parallel=control$parallel,
                                  parallel.type=control$parallel.type,
                                  parallel.version.check=control$parallel.version.check,
                                  term.options=control$term.options),
                                nsim=ceiling(control$obs.MCMC.samplesize/control$nsteps), ...)))
    }
  }
  message(".")
    
  Dtheta.Du<-to-from

  llrs<--sapply(seq_len(control$nsteps), function(i) crossprod(Dtheta.Du,ergm.etagradmult(path[i,],stats[i,]-stats.obs[i,],m$etamap)))/control$nsteps
  llr<-sum(llrs)
  if(llronly) llr
  else list(llr=llr,from=from,to=to,llrs=llrs,path=path,stats=stats,stats.obs=stats.obs,Dtheta.Du=Dtheta.Du)
}

#' @rdname ergm.bridge.llr
#'
#' @description \code{ergm.bridge.0.llk} is a convenience wrapper that
#'   returns the log-likelihood of configuration \eqn{\theta}
#'   \emph{relative to the reference measure}. That is, the
#'   configuration with \eqn{\theta=0} is defined as having log-likelihood of
#'   0.
#'
#' @return \code{ergm.bridge.0.llk} result list also includes an `llk`
#'   element, with the log-likelihood itself (with the reference
#'   distribution assumed to have likelihood 0).
#' 
#' 
#' @export
ergm.bridge.0.llk<-function(object, response=response, constraints=~., coef, ..., llkonly=TRUE, control=control.ergm.bridge()){
  check.control.class("ergm.bridge", "ergm.bridge.0.llk")
  control.toplevel(...,myname="ergm.bridge")

  br<-ergm.bridge.llr(object, constraints=constraints, from=rep(0,length(coef)), to=coef, response=response, control=control, ...)
  if(llkonly) br$llr
  else c(br,llk=br$llr)
}

#' @rdname ergm.bridge.llr
#'
#' @description `ergm.bridge.dindstart.llk` is a wrapper that uses a
#'   dyad-independent ERGM as a starting point for bridge sampling to
#'   estimate the log-likelihood for a given dyad-dependent model and
#'   parameter configuration.  Note that it only handles binary ERGMs
#'   (`response=NULL`) and with constraints (`constraints=`) that that
#'   do not induce dyadic dependence.
#'
#' @param dind A one-sided formula with the dyad-independent model to use as a
#' starting point. Defaults to the dyad-independent terms found in the formula
#' \code{object} with an overal density term (\code{edges}) added if not
#' redundant.
#' @param coef.dind Parameter configuration for the dyad-independent starting
#' point. Defaults to the MLE of \code{dind}.
#'
#' @return \code{ergm.bridge.dindstart.llk} result list also includes
#'   an `llk` element, with the log-likelihood itself and an
#'   `llk.dind` element, with the log-likelihood of the nearest
#'   dyad-independent model.
#' 
#' @export
ergm.bridge.dindstart.llk<-function(object, response=NULL, constraints=~., coef, dind=NULL, coef.dind=NULL,  basis=NULL, ..., llkonly=TRUE, control=control.ergm.bridge()){
  check.control.class("ergm.bridge", "ergm.bridge.dindstart.llk")
  control.toplevel(...,myname="ergm.bridge")

  if(!is.null(response)) stop("Only binary ERGMs are supported at this time.")

  ## Here, we need to get the model object to get the list of
  ## dyad-independent terms.
  tmp<-ergm.bridge.preproc(object,basis,response,term.options=control$term.options)
  nw<-tmp$nw; m<-tmp$model; form<-tmp$form; rm(tmp)

  p.pos.full <- c(0,cumsum(nparam(m, byterm=TRUE)))
  
  constraints.obs <- if(network.naedgecount(nw))
    nonsimp_update.formula(constraints,~.+observed)
  else
    NULL

  if(!is.dyad.independent(ergm_conlist(constraints,nw), ergm_conlist(constraints.obs,nw))) stop("Bridge sampling with dyad-independent start does not work with dyad-dependent constraints.")

  ## By default, take dyad-independent terms in the formula, fit a
  ## model with these terms and "edges". Terms that are redundant (NA)
  ## get their coefficients zeroed out below.
  offset.dind <- c()
  if(is.null(dind)){
    dind<-~.
    terms.full<-list_rhs.formula(form)[!m$term.skipped] # Ensure that terms to be added to the dyad-independent formula are aligned with terms that had actually made it into the model.
    for(i in seq_along(terms.full))
      if(NVL(m$terms[[i]]$dependence, TRUE) == FALSE){
        dind<-append_rhs.formula(dind,list(terms.full[[i]]))
        if(m$offset[i]) offset.dind <- c(offset.dind, coef[(p.pos.full[i]+1):p.pos.full[i+1]])
      }
    dind<-append_rhs.formula(dind,list(as.name("edges")))
    environment(dind) <- environment(object)
  }
  
  dind<-nonsimp_update.formula(dind,nw~., from.new="nw")

  if(!is.dyad.independent(dind, term.options=control$term.options))
    stop("Reference model `dind' must be dyad-independent.")

  ergm.dind<-suppressMessages(suppressWarnings(ergm(dind,estimate="MPLE",constraints=constraints,eval.loglik=FALSE,control=control.ergm(drop=FALSE, term.options=control$term.options, MPLE.max.dyad.types=control$MPLE.max.dyad.types), offset.coef = offset.dind)))
  
  if(is.null(coef.dind)){
    coef.dind<-ifelse(is.na(coef(ergm.dind)),0,coef(ergm.dind))
    llk.dind<--ergm.dind$glm$deviance/2 - -ergm.dind$glm$null.deviance/2
  }else{
    lin.pred <- model.matrix(ergm.dind$glm) %*% coef.dind
    llk.dind<-
      crossprod(lin.pred,ergm.dind$glm$y*ergm.dind$glm$prior.weights)-sum(log1p(exp(lin.pred))*ergm.dind$glm$prior.weights) -
      sum(as.rlebdm(ergm.dind$constrained, ergm.dind$constrained.obs, which="informative"))*log(1/2)
  }

  ## Construct the augmented formula.
  form.aug<-append_rhs.formula(object, list_rhs.formula(dind))

  from<-c(rep(0,length(coef)),coef.dind)
  to<-c(coef,rep(0,length(coef.dind)))

  br<-ergm.bridge.llr(form.aug, response=response, constraints=constraints, from=from, to=to, basis=basis, control=control)
  
  if(llkonly) llk.dind + br$llr
  else c(br,llk.dind=llk.dind, llk=llk.dind + br$llr)
}


