#  File R/ergm.bridge.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################

## This is a helper function that constructs and returns the network
## object to be used and the model object.
ergm.bridge.preproc<-function(object, basis, response){

  # If basis is not null, replace network in formula by basis.
  # In either case, let nw be network object from formula.
  if(is.null(nw <- basis)) {
    nw <- ergm.getnetwork(object)
  }

  nw <- as.network(nw)
  if(!is.network(nw)){
    stop("A network object on the LHS of the formula or via",
         " the 'basis' argument must be given")
  }
  
  # New formula (no longer use 'object'):
  form <- ergm.update.formula(object, nw ~ ., from.new="nw")
  
  list(nw=nw, form=form)
}


## The workhorse function: Uses bridge sampling to estimate the
## log-likelihood-ratio between two configurations `to' and `from' for
## a formula `object', using `nsteps' MCMC samples. If llronly==TRUE,
## returns only the estimate. Otherwise, returns a list with more
## details. Other parameters are same as simulate.ergm.
ergm.bridge.llr<-function(object, response=NULL, reference=~Bernoulli, constraints=~., from, to, obs.constraints=~-observed, target.stats=NULL, basis=NULL, verbose=FALSE, ..., llronly=FALSE, control=control.ergm.bridge()){
  check.control.class("ergm.bridge", "ergm.bridge.llr")
  control.toplevel(..., myname="ergm.bridge")

  if(!is.null(control$seed)) {set.seed(as.integer(control$seed))}
  
  ## Here, we need to get the model object to get the likelihood and gradient functions.
  tmp<-ergm.bridge.preproc(object,basis,response)
  nw<-tmp$nw; form<-tmp$form; rm(tmp)

  ## Generate the path.
  path<-t(rbind(sapply(seq(from=0+1/2/(control$nsteps+1),to=1-1/2/(control$nsteps+1),length.out=control$nsteps),function(u) cbind(to*u + from*(1-u)))))

  tmp <- .handle.obs.constraints(nw, constraints, obs.constraints, target.stats)
  nw <- tmp$nw
  constraints.obs <- tmp$constraints.obs
  

  ## Preinitialize MHproposals and set "observed" statistics:
  MHproposal <- MHproposal(constraints,arguments=control$MCMC.prop.args,
                           nw=nw, weights=control$MCMC.prop.weights, class="c",reference=reference,response=response)  
  m<-ergm.getmodel(object, nw, response=response, extra.aux=list(MHproposal$auxiliaries))

  if(!is.null(constraints.obs)){
    MHproposal.obs <- MHproposal(constraints.obs,arguments=control$obs.MCMC.prop.args,
                                 nw=nw, weights=control$obs.MCMC.prop.weights, class="c",reference=reference,response=response)
    m.obs<-ergm.getmodel(object, nw, response=response, extra.aux=list(MHproposal.obs$auxiliaries))

    stats.obs <- matrix(NA,control$nsteps,m$etamap$etalength)
  }else
    stats.obs<-matrix(NVL(target.stats,ergm.getglobalstats(nw, m, response=response)),control$nsteps,m$etamap$etalength,byrow=TRUE)

  stats<-matrix(NA,control$nsteps,m$etamap$etalength)
  
  message("Using ", control$nsteps, " bridges: ", appendLF=FALSE)
  
  for(i in seq_len(control$nsteps)){
    theta<-path[i,]
    if(verbose==0) message(i," ",appendLF=FALSE)
    if(verbose>0) message("Running theta=[",paste(format(theta),collapse=","),"].")
    if(verbose>1) message("Burning in...")
    ## First burn-in has to be longer, but those thereafter should be shorter if the bridges are closer together.
    nw.state<-simulate(m, coef=theta, nsim=1, response=response, reference=reference, constraints=MHproposal, basis=nw, statsonly=FALSE, verbose=max(verbose-1,0),
                       control=control.simulate.formula(MCMC.burnin=if(i==1) control$MCMC.burnin else ceiling(control$MCMC.burnin/sqrt(control$nsteps)),
                         MCMC.interval=1,
                         MCMC.packagenames=control$MCMC.packagenames,
                         parallel=control$parallel,
                         parallel.type=control$parallel.type,
                         parallel.version.check=control$parallel.version.check), ...)

    stats[i,]<-colMeans(simulate(m, coef=theta, response=response, reference=reference, constraints=MHproposal, basis=nw.state, statsonly=TRUE, verbose=max(verbose-1,0),
                              control=control.simulate.formula(MCMC.burnin=0,
                                                               MCMC.interval=control$MCMC.interval,
                                                               MCMC.packagenames=control$MCMC.packagenames,
                                                               parallel=control$parallel,
                                                               parallel.type=control$parallel.type,
                                                               parallel.version.check=control$parallel.version.check),
                              nsim=ceiling(control$MCMC.samplesize/control$nsteps), ...))
    
    if(!is.null(constraints.obs)){
      nw.state.obs<-simulate(m.obs, coef=theta, nsim=1, response=response, reference=reference, constraints=MHproposal.obs, basis=nw, statsonly=FALSE, verbose=max(verbose-1,0),
                             control=control.simulate.formula(MCMC.burnin=if(i==1) control$obs.MCMC.burnin else ceiling(control$obs.MCMC.burnin/sqrt(control$nsteps)),
                               MCMC.interval=1,
                               MCMC.packagenames=control$MCMC.packagenames,
                               parallel=control$parallel,
                               parallel.type=control$parallel.type,
                               parallel.version.check=control$parallel.version.check), ...)

      stats.obs[i,]<-colMeans(simulate(m.obs, coef=theta, response=response, reference=reference, constraints=MHproposal.obs, basis=nw.state.obs, statsonly=TRUE, verbose=max(verbose-1,0),
                                control=control.simulate.formula(MCMC.burnin=0,
                                  MCMC.interval=control$obs.MCMC.interval,
                                  MCMC.packagenames=control$MCMC.packagenames,
                                  parallel=control$parallel,
                                  parallel.type=control$parallel.type,
                                  parallel.version.check=control$parallel.version.check),
                                nsim=ceiling(control$obs.MCMC.samplesize/control$nsteps), ...))
    }
  }
  message(".")
    
  Dtheta.Du<-to-from

  llrs<--sapply(seq_len(control$nsteps), function(i) crossprod(Dtheta.Du,ergm.etagradmult(path[i,],stats[i,]-stats.obs[i,],m$etamap)))/control$nsteps
  llr<-sum(llrs)
  if(llronly) llr
  else list(llr=llr,from=from,to=to,llrs=llrs,path=path,stats=stats,stats.obs=stats.obs,Dtheta.Du=Dtheta.Du)
}

## A convenience wrapper around ergm.bridge.llr: returns the
## log-likelihood of configuration `theta' *relative to the reference
## measure*. That is, the configuration with theta=0 is defined as
## having log-likelihood of 0.
ergm.bridge.0.llk<-function(object, response=NULL, reference=~Bernoulli, coef, ..., llkonly=TRUE, control=control.ergm.bridge()){
  check.control.class("ergm.bridge", "ergm.bridge.0.llk")
  control.toplevel(...,myname="ergm.bridge")
  br<-ergm.bridge.llr(object, from=rep(0,length(coef)), to=coef, response=response, reference=reference, control=control, ...)
  if(llkonly) br$llr
  else c(br,llk=br$llr)
}

## A wrapper around ergm.bridge.llr that uses a specified
## dyad-independence model `dind' (specified as RHS-only formula),
## either at the its MLE (the default) or at a value specified by
## coef.dind, as a starting point for the bridge sampling. The terms
## in the dyad-independent model may overlap with the terms in the
## model whose likelihood is being evaluated, but don't have to.
## `dind' defaults to the dyad-independent terms of the `object'
## formula with an edges term added unless redundant.
ergm.bridge.dindstart.llk<-function(object, response=NULL, constraints=~., coef, obs.constraints=~-observed, target.stats=NULL, dind=NULL, coef.dind=NULL,  basis=NULL, ..., llkonly=TRUE, control=control.ergm.bridge()){
  check.control.class("ergm.bridge", "ergm.bridge.dindstart.llk")
  control.toplevel(...,myname="ergm.bridge")

  if(!is.null(response)) stop("Only binary ERGMs are supported at this time.")

  ## Here, we need to get the model object to get the list of
  ## dyad-independent terms.
  tmp<-ergm.bridge.preproc(object,basis,response)
  nw<-tmp$nw; form<-tmp$form; rm(tmp)

  m<-ergm.getmodel(object, nw, response=response)
  
  q.pos.full <- c(0,cumsum(coef.sublength.model(m)))
  p.pos.full <- c(0,cumsum(eta.sublength.model(m)))
  
  tmp <- .handle.obs.constraints(nw, constraints, obs.constraints, target.stats)
  nw <- tmp$nw
  constraints.obs <- tmp$constraints.obs
 
  if(!is.dyad.independent(ergm_conlist(constraints,nw), ergm_conlist(constraints.obs,nw))) stop("Bridge sampling with dyad-independent start does not work with dyad-dependent constraints.")

  # If target.stats are given, then we need between passed network and
  # target stats, if any. It also means that the dyad-independent
  # submodel cannot contain any statistics that target.stats does not,
  # so edges are not added on.
  if(!is.null(target.stats) && !is.null(dind)) stop("Non-default dyad-independent model is not supported when target.stats is passed and passed network's statistics do not match it.")

  ## From this point on, target.stats has NAs corresponding to offset
  ## terms.
  if(!is.null(target.stats)) target.stats <- .align.target.stats.offset(m, target.stats)
  
  ## By default, take dyad-independent terms in the formula, fit a
  ## model with these terms and "edges". Terms that are redundant (NA)
  ## get their coefficients zeroed out below.
  ## FIXME: What to do about dyad-independent curved terms?
  offset.dind <- c()
  if(!is.null(target.stats)) ts.dind <- c()
  if(is.null(dind)){
    dind<-~.
    terms.full<-term.list.formula(form[[3]])[!m$term.skipped] # Ensure that terms to be added to the dyad-independent formula are aligned with terms that had actually made it into the model.
    for(i in seq_along(terms.full))
      if(NVL(m$terms[[i]]$dependence, TRUE) == FALSE){
        dind<-append.rhs.formula(dind,list(terms.full[[i]]))
        if(!is.null(target.stats) && !m$offset[i]) ts.dind <- c(ts.dind, target.stats[(p.pos.full[i]+1):p.pos.full[i+1]])
        if(m$offset[i]) offset.dind <- c(offset.dind, coef[(q.pos.full[i]+1):q.pos.full[i+1]])
      }
    if(is.null(target.stats)) dind<-append.rhs.formula(dind,list(as.name("edges")))
    environment(dind) <- environment(object)
  }
  
  dind<-ergm.update.formula(dind,nw~., from.new="nw")

  if(!is.dyad.independent(dind))
    stop("Reference model `dind' must be dyad-independent.")

  ergm.dind<-suppressMessages(suppressWarnings(ergm(dind,estimate="MPLE",constraints=constraints,eval.loglik=FALSE,control=control.ergm(drop=FALSE, MPLE.max.dyad.types=control$MPLE.max.dyad.types), offset.coef = offset.dind)))
  
  if(is.null(coef.dind)){
    coef.dind <- coef(ergm.dind)[!ergm.dind$etamap$offsettheta]
    coef.dind <- ifelse(is.na(coef.dind),0,coef.dind)
    llk.dind<--ergm.dind$glm$deviance/2 - -ergm.dind$glm$null.deviance/2
  }else{
    lin.pred <- model.matrix(ergm.dind$glm) %*% coef.dind
    llk.dind <- 
      crossprod(lin.pred,ergm.dind$glm$y*ergm.dind$glm$prior.weights)-sum(log1p(exp(lin.pred))*ergm.dind$glm$prior.weights) -
        (network.dyadcount(ergm.dind$network,FALSE) - network.edgecount(NVL(as.rlebdm(ergm.dind$constrained, ergm.dind$constrained.obs,which="missing"),network.initialize(1))))*log(1/2)
  }
  
  # If there are target.stats we need to adjust the log-likelihood in
  # case they are different from those to which the dyad-independent
  # submodel was actually fit:
  # l(theta,ts)-l(theta,ns)=sum(theta*(ts-ns)).
  if(!is.null(target.stats)) llk.dind <- llk.dind + c(crossprod(coef.dind, NVL(c(ts.dind), ergm.dind$nw.stats[!ergm.dind$etamap$offsetmap]) - ergm.dind$nw.stats[!ergm.dind$etamap$offsetmap]))

  ## Construct the augmented formula.
  form.aug<-append.rhs.formula(object, term.list.formula(dind[[3]])[!ergm.dind$etamap$offset])

  from<-c(rep(0,length(coef)),coef.dind)
  to<-c(coef,rep(0,length(coef.dind)))

  if(!is.null(target.stats) && any(is.na(target.stats))){
    warning("Using target.stats for a model with offset terms may produce an inaccurate estimate of the log-likelihood and derived quantities (deviance, AIC, BIC, etc.), because some of the target stats must be imputed.")
    target.stats[m$etamap$offsetmap] <- ergm.getglobalstats(nw, m)[m$etamap$offsetmap]
  }

  br<-ergm.bridge.llr(form.aug, response=response, constraints=constraints, from=from, to=to, basis=basis, target.stats=c(target.stats, if(!is.null(target.stats)) ts.dind), control=control)
  
  if(llkonly) llk.dind + br$llr
  else c(br,llk.dind=llk.dind, llk=llk.dind + br$llr)
}


