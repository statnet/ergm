
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
  form <- ergm.update.formula(object, nw ~ .)
  
  list(nw=nw, form=form, model=ergm.getmodel(form, nw, response=response))
}


## The workhorse function: Uses bridge sampling to estimate the
## log-likelihood-ratio between two configurations `to' and `from' for
## a model `object', using `nsteps' MCMC samples. If llronly==TRUE,
## returns only the estimate. Otherwise, returns a list with more
## details. Other parameters are same as simulate.ergm.
ergm.bridge.llr<-function(object, response=NULL, from, to, nsteps=20, sample.size=10000, burnin=10000, basis=NULL, verbose=FALSE, llronly=FALSE, ...){

  ## Here, we need to get the model object to get the likelihood and gradient functions.
  tmp<-ergm.bridge.preproc(object,basis,response)
  nw<-tmp$nw; m<-tmp$model; form<-tmp$form; rm(tmp)

  ## Generate the path.

  path<-t(rbind(sapply(seq(from=0+1/2/(nsteps+1),to=1-1/2/(nsteps+1),length.out=nsteps),function(u) cbind(to*u + from*(1-u)))))

  obs<-summary(form,response=response)

  stats<-matrix(NA,nsteps,length(obs))


  for(i in seq_len(nsteps)){
    theta<-path[i,]
    if(verbose) cat("Running theta=[",paste(format(theta),collapse=","),"].\n",sep="")
    if(verbose>1) cat("Burning in...\n",sep="")
    ## First burn-in has to be longer, but those thereafter should be shorter if the bridges are closer together.
    nw.state<-simulate(form, theta0=theta, nsim=1, response=response, basis=basis, statsonly=FALSE, verbose=max(verbose-1,0), burnin=if(i==1) burnin else ceiling(burnin/sqrt(nsteps)), interval=1, ...)
    ergm.update.formula(form,nw.state~.)
    stats[i,]<-apply(simulate(form, theta0=theta, response=response, basis=basis, statsonly=TRUE, verbose=max(verbose-1,0), burnin=0, nsim=ceiling(sample.size/nsteps), ...),2,mean)-obs
  }
    
  Dtheta.Du<-to-from

  llrs<--sapply(seq_len(nsteps), function(i) crossprod(Dtheta.Du,ergm.etagradmult(path[i,],stats[i,],m$etamap)))/nsteps
  llr<-sum(llrs)
  if(llronly) llr
  else list(llr=llr,llrs=llrs,path=path,stats=stats,Dtheta.Du=Dtheta.Du)
}

## A convenience wrapper around ergm.bridge.llr: returns the
## log-likelihood of configuration `theta' *relative to the reference
## measure*. That is, the configuration with theta=0 is defined as
## having log-likelihood of 0.
ergm.bridge.0.llk<-function(object, response=response, theta, nsteps=20, llkonly=TRUE, ...){
  br<-ergm.bridge.llr(object, from=rep(0,length(theta)), to=theta, nsteps=nsteps, response=response, ...)
  if(llkonly) br$llr
  else c(br,llk=br$llr)
}

## A wrapper around ergm.bridge.llr that uses a specified
## dyad-independence model `dind' (specified as RHS-only formula),
## either at the its MLE (the default) or at a value specified by
## theta.dind, as a starting point for the bridge sampling. The terms
## in the dyad-independent model may overlap with the terms in the
## model whose likelihood is being evaluated, but don't have to.
## `dind' defaults to the dyad-independent terms of the `object'
## formula with an edges term added unless redundant.
ergm.bridge.dindstart.llk<-function(object, response=NULL, theta, nsteps=20, dind=NULL, theta.dind=NULL,  basis=NULL, llkonly=TRUE, ...){
  if(!is.null(response)) stop("Only binary ERGMs are supported at this time.")

  ## Here, we need to get the model object to get the list of
  ## dyad-independent terms.
  tmp<-ergm.bridge.preproc(object,basis,response)
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
    ergm.dind<-ergm(dind,MPLEonly=TRUE)
    ## If any terms are redundant...
    if(any(is.na(coef(ergm.dind)))){
      dind<-~nw
      terms.dind<-term.list.formula(ergm.dind$formula[[3]])
      for(i in seq_along(terms.dind))
        if(!is.na(coef(ergm.dind)[i]))
          dind<-append.rhs.formula(dind,list(terms.dind[[i]]))
      ergm.dind<-ergm(dind,MPLEonly=TRUE)
    }
  }else{
    dind<-ergm.update.formula(dind,nw~.)  
    ergm.dind<-ergm(dind,MPLEonly=TRUE)
  }  
  
  if(!is.dyad.independent(dind))
    stop("Reference model `dind' must be dyad-independent.")

  if(is.null(theta.dind)){
    coef.dind<-coef(ergm.dind)
    llk.dind<--ergm.dind$glm$deviance/2
  }else{
    coef.dind<-theta.dind
    lin.pred <- model.matrix(ergm.dind$glm) %*% theta.dind
    llk.dind<- crossprod(lin.pred,ergm.dind$glm$y*ergm.dind$glm$prior.weights)-sum(log1p(exp(lin.pred))*ergm.dind$glm$prior.weights)
  }
  

  ## Construct the augmented formula.
  form.aug<-append.rhs.formula(object, term.list.formula(dind[[3]]))

  from<-c(rep(0,length(theta)),coef.dind)
  to<-c(theta,rep(0,length(coef.dind)))

  br<-ergm.bridge.llr(form.aug, response=response, from=from, to=to, basis=basis, nsteps=nsteps, ...)
  
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
ergm.bridge.hammingstart.llk<-function(object, response=NULL, theta, nsteps, hamming.start=NULL, llk.guess=NULL, basis=NULL, llkonly=TRUE, ...){
  if(!is.null(response)) stop("Only binary ERGMs are supported at this time.")
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
  from<-c(rep(0,length(theta)), hamming.start)
  to<-c(theta,0)
  
  llk.hamming<--network.dyadcount(nw)*log1p(exp(hamming.start))
  br<-ergm.bridge.llr(form.aug, response=response, from=from, to=to, basis=basis, nsteps=nsteps, ...)

  if(llkonly) llk.hamming + br$llr
  else c(br,llk.hamming=llk.hamming, llk=llk.hamming + br$llr) 
}

