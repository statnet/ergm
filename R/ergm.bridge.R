#  File ergm/R/ergm.bridge.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################

## This is a helper function that constructs and returns the network
## object to be used and the model object.
ergm.bridge.preproc<-function(object, basis){

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
  form <- ergm.update.formula(object, nw ~ .)
  
  list(nw=nw, form=form, model=ergm.getmodel(form, nw))
}


## The workhorse function: Uses bridge sampling to estimate the
## log-likelihood-ratio between two configurations `to' and `from' for
## a model `object', using `nsteps' MCMC samples. If llronly==TRUE,
## returns only the estimate. Otherwise, returns a list with more
## details. Other parameters are same as simulate.ergm.
ergm.bridge.llr<-function(object, constraints=~., from, to, basis=NULL, verbose=FALSE, ..., llronly=FALSE, control=control.ergm.bridge()){

  if(!is.null(control$seed)) {set.seed(as.integer(control$seed))}
  if(!is.null(basis)) ergm.update.formula(form,basis~.)
  
  ## Here, we need to get the model object to get the likelihood and gradient functions.
  tmp<-ergm.bridge.preproc(object,basis)
  nw<-tmp$nw; m<-tmp$model; form<-tmp$form; rm(tmp)


  ## Generate the path.
  path<-t(rbind(sapply(seq(from=0+1/2/(control$nsteps+1),to=1-1/2/(control$nsteps+1),length.out=control$nsteps),function(u) cbind(to*u + from*(1-u)))))

  stats<-matrix(NA,control$nsteps,m$etamap$etalength)
  
  if(network.naedgecount(nw)){
    constraints.obs<-ergm.update.formula(constraints,~.+observed)
    form.obs<-form
    stats.obs <- matrix(NA,control$nsteps,m$etamap$etalength)  
  }else stats.obs<-matrix(summary(form),control$nsteps,m$etamap$etalength,byrow=TRUE)  

  for(i in seq_len(control$nsteps)){
    theta<-path[i,]
    if(verbose) cat("Running theta=[",paste(format(theta),collapse=","),"].\n",sep="")
    if(verbose>1) cat("Burning in...\n",sep="")
    ## First burn-in has to be longer, but those thereafter should be shorter if the bridges are closer together.
    nw.state<-simulate(form, coef=theta, nsim=1, constraints=constraints, statsonly=FALSE, verbose=max(verbose-1,0),
                       control=control.simulate.ergm(MCMC.burnin=if(i==1) control$MCMC.burnin else ceiling(control$MCMC.burnin/sqrt(control$nsteps)),
                         MCMC.interval=1,
                         MCMC.prop.args=control$MCMC.prop.args,
                         MCMC.prop.weights=control$MCMC.prop.weights,
                         MCMC.packagenames=control$MCMC.packagenames))
    ergm.update.formula(form,nw.state~.)
    stats[i,]<-apply(simulate(form, coef=theta, constraints=constraints, statsonly=TRUE, verbose=max(verbose-1,0),
                              control=control.simulate.ergm(MCMC.burnin=0,
                                MCMC.interval=control$MCMC.interval),
                              nsim=ceiling(control$MCMC.samplesize/control$nsteps)),2,mean)
    
    if(network.naedgecount(nw)){
      nw.state.obs<-simulate(form.obs, coef=theta, nsim=1, constraints=constraints.obs, statsonly=FALSE, verbose=max(verbose-1,0),
                             control=control.simulate.ergm(MCMC.burnin=if(i==1) control$obs.MCMC.burnin else ceiling(control$obs.MCMC.burnin/sqrt(control$nsteps)),
                               MCMC.interval=1,
                               MCMC.prop.args=control$MCMC.prop.args,
                               MCMC.prop.weights=control$MCMC.prop.weights,
                               MCMC.packagenames=control$MCMC.packagenames))
      ergm.update.formula(form.obs,nw.state.obs~.)
      stats.obs[i,]<-apply(simulate(form.obs, coef=theta, constraints=constraints.obs, statsonly=TRUE, verbose=max(verbose-1,0),
                                control=control.simulate.ergm(MCMC.burnin=0,
                                  MCMC.interval=control$obs.MCMC.interval),
                                nsim=ceiling(control$obs.MCMC.samplesize/control$nsteps)),2,mean)
    }
  }
    
  Dtheta.Du<-to-from

  llrs<--sapply(seq_len(control$nsteps), function(i) crossprod(Dtheta.Du,ergm.etagradmult(path[i,],stats[i,]-stats.obs[i,],m$etamap)))/control$nsteps
  llr<-sum(llrs)
  if(llronly) llr
  else list(llr=llr,llrs=llrs,path=path,stats=stats,stats.obs=stats.obs,Dtheta.Du=Dtheta.Du)
}

## A convenience wrapper around ergm.bridge.llr: returns the
## log-likelihood of configuration `theta' *relative to the reference
## measure*. That is, the configuration with theta=0 is defined as
## having log-likelihood of 0.
ergm.bridge.0.llk<-function(object, coef, ..., llkonly=TRUE, control=control.ergm.bridge()){
  br<-ergm.bridge.llr(object, from=rep(0,length(coef)), to=coef, control=control)
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
ergm.bridge.dindstart.llk<-function(object, coef, dind=NULL, coef.dind=NULL,  basis=NULL, ..., llkonly=TRUE, control=control.ergm.bridge()){
  ## Here, we need to get the model object to get the list of
  ## dyad-independent terms.
  tmp<-ergm.bridge.preproc(object,basis)
  nw<-tmp$nw; m<-tmp$model; form<-tmp$form; rm(tmp)

  ## By default, take dyad-independent terms in the formula, fit a
  ## model with these terms and "edges", then drop the terms that get
  ## NAs (i.e. are redundant).
   if(is.null(dind)){
    dind<-nw~edges
    terms.full<-term.list.formula(form[[3]])
    for(i in seq_along(terms.full))
      if(!is.null(m$terms[[i]]$dependence) && m$terms[[i]]$dependence==FALSE)
        dind<-append.rhs.formula(dind,list(terms.full[[i]]))
    ergm.dind<-ergm(dind,estimate="MPLE")
    ## If any terms are redundant...
    if(any(is.na(coef(ergm.dind)))){
      dind<-~nw
      terms.dind<-term.list.formula(ergm.dind$formula[[3]])
      for(i in seq_along(terms.dind))
        if(!is.na(coef(ergm.dind)[i]))
          dind<-append.rhs.formula(dind,list(terms.dind[[i]]))
      ergm.dind<-ergm(dind,estimate="MPLE")
    }
  }else{
    dind<-ergm.update.formula(dind,nw~.)  
    ergm.dind<-ergm(dind,estimate="MPLE")
  }  
  
  if(!is.dyad.independent(dind))
    stop("Reference model `dind' must be dyad-independent.")

  if(is.null(coef.dind)){
    coef.dind<-coef(ergm.dind)
    llk.dind<--ergm.dind$glm$deviance/2
  }else{
    lin.pred <- model.matrix(ergm.dind$glm) %*% coef.dind
    llk.dind<- crossprod(lin.pred,ergm.dind$glm$y*ergm.dind$glm$prior.weights)-sum(log1p(exp(lin.pred))*ergm.dind$glm$prior.weights)
  }
  

  ## Construct the augmented formula.
  form.aug<-append.rhs.formula(object, term.list.formula(dind[[3]]))

  from<-c(rep(0,length(coef)),coef.dind)
  to<-c(coef,rep(0,length(coef.dind)))

  br<-ergm.bridge.llr(form.aug, from=from, to=to, basis=basis, control=control)
  
  if(llkonly) llk.dind + br$llr
  else c(br,llk.dind=llk.dind, llk=llk.dind + br$llr)
}

## A wrapper around ergm.bridge.llr that uses a model with a Hamming
## distance to the LHS network itself as a starting point, either with
## a specified coefficient `hamming.start' or with a coefficient such
## that the log-likelihood for it is llk.guess.
##
## The idea is to use the Hamming term as "scaffolding", which is
## slowly removed as the real model terms approach their objective
## values.
ergm.bridge.hammingstart.llk<-function(object, coef, hamming.start=NULL, llk.guess=NULL, basis=NULL, ..., llkonly=TRUE, control=control.ergm.bridge()){
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

  if(is.null(hamming.start)){
    if(is.null(llk.guess))  llk.guess<-ergm(nw~edges)$mle.lik

    hamming.start<-log(expm1(-llk.guess/network.dyadcount(nw)))
  }

  form.aug<-ergm.update.formula(object, . ~ . + hamming(nw))
  from<-c(rep(0,length(coef)), hamming.start)
  to<-c(coef,0)
  
  llk.hamming<--network.dyadcount(nw)*log1p(exp(hamming.start))
  br<-ergm.bridge.llr(form.aug, from=from, to=to, basis=basis, control=control)

  if(llkonly) llk.hamming + br$llr
  else c(br,llk.hamming=llk.hamming, llk=llk.hamming + br$llr) 
}

