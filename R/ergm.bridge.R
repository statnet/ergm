#  File R/ergm.bridge.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
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
  
  list(nw=nw, form=form, model=ergm.getmodel(form, nw, response=response))
}


## The workhorse function: Uses bridge sampling to estimate the
## log-likelihood-ratio between two configurations `to' and `from' for
## a model `object', using `nsteps' MCMC samples. If llronly==TRUE,
## returns only the estimate. Otherwise, returns a list with more
## details. Other parameters are same as simulate.ergm.
ergm.bridge.llr<-function(object, response=NULL, constraints=~., from, to, basis=NULL, verbose=FALSE, ..., llronly=FALSE, control=control.ergm.bridge()){
  check.control.class("ergm.bridge")

  if(!is.null(control$seed)) {set.seed(as.integer(control$seed))}
  if(!is.null(basis)) ergm.update.formula(form,basis~., from.new="basis")
  
  ## Here, we need to get the model object to get the likelihood and gradient functions.
  tmp<-ergm.bridge.preproc(object,basis,response)
  nw<-tmp$nw; m<-tmp$model; form<-tmp$form; rm(tmp)


  ## Generate the path.
  path<-t(rbind(sapply(seq(from=0+1/2/(control$nsteps+1),to=1-1/2/(control$nsteps+1),length.out=control$nsteps),function(u) cbind(to*u + from*(1-u)))))

  stats<-matrix(NA,control$nsteps,m$etamap$etalength)
  
  if(network.naedgecount(nw)){
    constraints.obs<-ergm.update.formula(constraints,~.+observed)
    form.obs<-form
    stats.obs <- matrix(NA,control$nsteps,m$etamap$etalength)  
  }else stats.obs<-matrix(summary(form,response=response),control$nsteps,m$etamap$etalength,byrow=TRUE)  

  for(i in seq_len(control$nsteps)){
    theta<-path[i,]
    if(verbose) cat("Running theta=[",paste(format(theta),collapse=","),"].\n",sep="")
    if(verbose>1) cat("Burning in...\n",sep="")
    ## First burn-in has to be longer, but those thereafter should be shorter if the bridges are closer together.
    nw.state<-simulate(form, coef=theta, nsim=1, response=response, constraints=constraints, statsonly=FALSE, verbose=max(verbose-1,0),
                       control=control.simulate.formula(MCMC.burnin=if(i==1) control$MCMC.burnin else ceiling(control$MCMC.burnin/sqrt(control$nsteps)),
                         MCMC.interval=1,
                         MCMC.prop.args=control$MCMC.prop.args,
                         MCMC.prop.weights=control$MCMC.prop.weights,
                         MCMC.packagenames=control$MCMC.packagenames), ...)
    ergm.update.formula(form,nw.state~., from.new="nw.state")
    stats[i,]<-apply(simulate(form, coef=theta, response=response, constraints=constraints, statsonly=TRUE, verbose=max(verbose-1,0),
                              control=control.simulate.formula(MCMC.burnin=0,
                                MCMC.interval=control$MCMC.interval),
                              nsim=ceiling(control$MCMC.samplesize/control$nsteps), ...),2,mean)
    
    if(network.naedgecount(nw)){
      nw.state.obs<-simulate(form.obs, coef=theta, nsim=1, response=response, constraints=constraints.obs, statsonly=FALSE, verbose=max(verbose-1,0),
                             control=control.simulate.formula(MCMC.burnin=if(i==1) control$obs.MCMC.burnin else ceiling(control$obs.MCMC.burnin/sqrt(control$nsteps)),
                               MCMC.interval=1,
                               MCMC.prop.args=control$MCMC.prop.args,
                               MCMC.prop.weights=control$MCMC.prop.weights,
                               MCMC.packagenames=control$MCMC.packagenames), ...)
      ergm.update.formula(form.obs,nw.state.obs~., from.new="nw.state.obs")
      stats.obs[i,]<-apply(simulate(form.obs, coef=theta, response=response, constraints=constraints.obs, statsonly=TRUE, verbose=max(verbose-1,0),
                                control=control.simulate.formula(MCMC.burnin=0,
                                  MCMC.interval=control$obs.MCMC.interval),
                                nsim=ceiling(control$obs.MCMC.samplesize/control$nsteps), ...),2,mean)
    }
  }
    
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
ergm.bridge.0.llk<-function(object, response=response, coef, ..., llkonly=TRUE, control=control.ergm.bridge()){
  check.control.class("ergm.bridge")
  br<-ergm.bridge.llr(object, from=rep(0,length(coef)), to=coef, response=response, control=control, ...)
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
ergm.bridge.dindstart.llk<-function(object, response=NULL, constraints=~., coef, dind=NULL, coef.dind=NULL,  basis=NULL, ..., llkonly=TRUE, control=control.ergm.bridge()){
  check.control.class("ergm.bridge")
  if(!is.null(response)) stop("Only binary ERGMs are supported at this time.")

  ## Here, we need to get the model object to get the list of
  ## dyad-independent terms.
  tmp<-ergm.bridge.preproc(object,basis,response)
  nw<-tmp$nw; m<-tmp$model; form<-tmp$form; rm(tmp)

  p.pos.full <- c(0,cumsum(coef.sublength.model(m)))
  
  constraints.obs <- if(network.naedgecount(nw))
    ergm.update.formula(constraints,~.+observed)
  else
    NULL

  if(!is.dyad.independent(mk.conlist(constraints,nw), mk.conlist(constraints.obs,nw))) stop("Bridge sampling with dyad-independent start does not work with dyad-dependent constraints.")

  ## By default, take dyad-independent terms in the formula, fit a
  ## model with these terms and "edges". Terms that are redundant (NA)
  ## get their coefficients zeroed out below.
  offset.dind <- c()
  if(is.null(dind)){
    dind<-~.
    terms.full<-term.list.formula(form[[3]])
    for(i in seq_along(terms.full))
      if(!is.null(m$terms[[i]]$dependence) && m$terms[[i]]$dependence==FALSE){
        dind<-append.rhs.formula(dind,list(terms.full[[i]]))
        if(m$offset[i]) offset.dind <- c(offset.dind, coef[(p.pos.full[i]+1):p.pos.full[i+1]])
      }
    dind<-append.rhs.formula(dind,list(as.name("edges")))
    environment(dind) <- environment(object)
  }
  
  dind<-ergm.update.formula(dind,nw~., from.new="nw")

  if(!is.dyad.independent(dind))
    stop("Reference model `dind' must be dyad-independent.")

  ergm.dind<-ergm(dind,estimate="MPLE",constraints=constraints,eval.loglik=FALSE,control=control.ergm(drop=FALSE, MPLE.max.dyad.types=control$MPLE.max.dyad.types), offset.coef = offset.dind)
  
  if(is.null(coef.dind)){
    coef.dind<-ifelse(is.na(coef(ergm.dind)),0,coef(ergm.dind))
    llk.dind<--ergm.dind$glm$deviance/2 - -ergm.dind$glm$null.deviance/2
  }else{
    lin.pred <- model.matrix(ergm.dind$glm) %*% coef.dind
    llk.dind<-
      crossprod(lin.pred,ergm.dind$glm$y*ergm.dind$glm$prior.weights)-sum(log1p(exp(lin.pred))*ergm.dind$glm$prior.weights) -
        (network.dyadcount(ergm.dind$network,FALSE) - network.edgecount(NVL(get.miss.dyads(ergm.dind$constrained, ergm.dind$constrained.obs),network.initialize(1))))*log(1/2)
  }  

  ## Construct the augmented formula.
  form.aug<-append.rhs.formula(object, term.list.formula(dind[[3]]))

  from<-c(rep(0,length(coef)),coef.dind)
  to<-c(coef,rep(0,length(coef.dind)))

  br<-ergm.bridge.llr(form.aug, response=response, constraints=constraints, from=from, to=to, basis=basis, control=control)
  
  if(llkonly) llk.dind + br$llr
  else c(br,llk.dind=llk.dind, llk=llk.dind + br$llr)
}

## ## A wrapper around ergm.bridge.llr that uses a model with a Hamming
## ## distance to the LHS network itself as a starting point, either with
## ## a specified coefficient `hamming.start' or with a coefficient such
## ## that the log-likelihood for it is llk.guess.
## ##
## ## The idea is to use the Hamming term as "scaffolding", which is
## ## slowly removed as the real model terms approach their objective
## ## values.
## ergm.bridge.hammingstart.llk<-function(object, response=NULL, coef, hamming.start=NULL, llk.guess=NULL, basis=NULL, ..., llkonly=TRUE, control=control.ergm.bridge()){
##   check.control.class("ergm.bridge")
##   if(!is.null(response)) stop("Only binary ERGMs are supported at this time.")
##   # If basis is not null, replace network in formula by basis.
##   # In either case, let nw be network object from formula.
##   if(is.null(nw <- basis)) {
##     nw <- ergm.getnetwork(object)
##   }
  
##   nw <- as.network(nw)
##   if(!is.network(nw)){
##     stop("A network object on the LHS of the formula or via",
##          " the 'basis' argument must be given")
##   }

##   if(is.null(hamming.start)){
##     if(is.null(llk.guess))  llk.guess<-logLik(ergm(nw~edges)$mle.lik

##     hamming.start<-log(expm1(-llk.guess/network.dyadcount(nw)))
##   }

##   form.aug<-ergm.update.formula(object, . ~ . + hamming(nw), from.new="nw")
##   from<-c(rep(0,length(coef)), hamming.start)
##   to<-c(coef,0)
  
##   llk.hamming<--network.dyadcount(nw)*log1p(exp(hamming.start))
##   br<-ergm.bridge.llr(form.aug, response=response, from=from, to=to, basis=basis, control=control)

##   if(llkonly) llk.hamming + br$llr
##   else c(br,llk.hamming=llk.hamming, llk=llk.hamming + br$llr) 
## }

