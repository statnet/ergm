ergm.mapl <- function(formula, theta0="MPLE", 
                 MPLEonly=TRUE, MLestimate=!MPLEonly, nsim=0,
                 burnin=10000000,
                 maxit=3,
                 constraints=~.,
                 proposaltype="TNT10",
                 meanstats=NULL,
                 control=ergm.control(MPLEtype="penalized"),
                 tau=1, invcov=NULL,
                 verbose=FALSE, ...) {
  current.warn <- options()$warn
  options(warn=0)
  if (verbose) cat("Evaluating network in model\n")

  nw <- ergm.getnetwork(formula)
  if(!is.null(meanstats)){ control$drop <- FALSE }
  
  if (verbose) cat("Fitting initial model.\n")

  proposalclass <- "c"
    
  if(control$drop){
   model.initial <- ergm.getmodel(formula, nw, drop=FALSE, initialfit=TRUE)
   model.initial.drop <- ergm.getmodel(formula, nw, drop=TRUE, initialfit=TRUE)
   namesmatch <- match(model.initial$coef.names, model.initial.drop$coef.names)
   droppedterms <- rep(FALSE, length=length(model.initial$etamap$offsettheta))
   droppedterms[is.na(namesmatch)] <- TRUE
   model.initial$etamap$offsettheta[is.na(namesmatch)] <- TRUE
  }else{
   model.initial <- ergm.getmodel(formula, nw, drop=control$drop, initialfit=TRUE)
   droppedterms <- rep(FALSE, length=length(model.initial$etamap$offsettheta))
  }
  MHproposal <- MHproposal(constraints, weights=control$prop.weights, control$prop.args, nw, model.initial,class=proposalclass)
  MHproposal.miss <- MHproposal("randomtoggleNonObserved", control$prop.args, nw, model.initial)

  # MPLE & Meanstats -> need fake network
  if(!missing(meanstats)){
    nw<-san(formula, meanstats=meanstats, 
            constraints=~., #constraints=constraints,
            proposaltype=proposaltype,
            tau=tau, invcov=invcov, burnin=burnin, verbose=verbose)
    if(verbose){print(summary(formula, basis=nw)-meanstats)}
  }
  
  Clist.initial <- ergm.Cprepare(nw, model.initial)
  Clist.miss.initial <- ergm.design(nw, model.initial, initialfit=TRUE,
                                    verbose=verbose)
  Clist.initial$meanstats=meanstats
  theta0copy <- theta0
  pl <- ergm.pl.ihs(Clist=Clist.initial, Clist.miss=Clist.miss.initial,
                    m=model.initial,
                    verbose=verbose)
  initialfit <- ergm.maple(pl=pl, model.initial,
                          MPLEtype=control$MPLEtype, 
                          verbose=verbose, ...)
  if("MPLE" %in% theta0){theta0 <- initialfit$coef}

  if(!missing(nsim)){
   nw.initial<-simulate(formula, constraints=constraints,
                       theta0=theta0,
                       verbose=verbose)

   meanstats <- summary(formula, basis=nw)
   sim <- nw.initial
   for(i in 1:nsim){
    sim <- san(formula, meanstats=meanstats, verbose=verbose,
               proposaltype=proposaltype,
               tau=tau, invcov=invcov, burnin=burnin, constraints=constraints, basis=sim)
    if(verbose){print(summary(formula, basis=sim)-meanstats)}
    if(verbose){print(sum(sim[,] != nw[,]))}
    Clist.initial <- ergm.Cprepare(sim, model.initial)
    Clist.miss.initial <- ergm.design(sim, model.initial, initialfit=TRUE,
                                verbose=verbose)
    Clist.initial$meanstats=meanstats
    sim.pl <- ergm.pl.ihs(Clist=Clist.initial, Clist.miss=Clist.miss.initial,
                     m=model.initial,
                     verbose=verbose)
    pl$zy <- c(pl$zy,sim.pl$zy)
    pl$foffset <- c(pl$foffset,sim.pl$foffset)
    pl$xmat <- rbind(pl$xmat,sim.pl$xmat)
    pl$wend <- c(pl$wend,sim.pl$wend)
    pl$zy.full <- c(pl$zy.full,sim.pl$zy.full)
    pl$foffset.full <- c(pl$foffset.full,sim.pl$foffset.full)
    pl$xmat.full <- rbind(pl$xmat.full,sim.pl$xmat.full)
    pl$wend.full <- c(pl$wend.full,sim.pl$wend.full)
   }
   pl$wend <- pl$wend / nsim
   pl$wend.full <- pl$wend.full / nsim
   initialfit <- ergm.maple(pl=pl, model.initial,
                            MPLEtype=control$MPLEtype, 
                            verbose=verbose, ...)
  }

    initialfit$offset <- model.initial$etamap$offsettheta
    initialfit$drop <- droppedterms
    initialfit$network <- nw
    initialfit$newnetwork <- nw
    initialfit$formula <- formula
    initialfit$constraints <- constraints
    initialfit$prop.args <- control$prop.args
    initialfit$prop.weights <- control$prop.weights
    initialfit
}
