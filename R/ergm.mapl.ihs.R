ergm.mapl <- function(formula, theta0="MPLE", 
                 MPLEonly=TRUE, MLestimate=!MPLEonly, nsim=0,
                 burnin=100000,
                 maxit=3,
                 constraints=~.,
                 meanstats=NULL,
                 control=ergm.control(MPLEtype="penalized"),
                 verbose=FALSE, ...) {
  current.warn <- options()$warn
  options(warn=0)
  if (verbose) cat("Evaluating network in model\n")

  nw <- ergm.getnetwork(formula)
  if(!is.null(meanstats)){ control$drop <- FALSE }
  if(control$nsubphases=="maxit") control$nsubphases<-maxit

  
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
  if("MPLE" %in% theta0 && !is.null(meanstats)){
  # if IHS 
    nw.initial<-san(formula, meanstats=meanstats, verbose=verbose)
  #else
  # nw.initial<-mk.pseudonet(meanstats, formula, nw, verbose=verbose)
  # IHS end
  }
  else nw.initial<-nw
  
  Clist.initial <- ergm.Cprepare(nw.initial, model.initial)
  Clist.miss.initial <- ergm.design(nw.initial, model.initial, initialfit=TRUE,
                                verbose=verbose)
  Clist.initial$meanstats=meanstats
  theta0copy <- theta0
  pl <- ergm.pl.ihs(Clist=Clist.initial, Clist.miss=Clist.miss.initial,
                    m=model.initial,
                    verbose=verbose)
  if(!missing(nsim)){
   meanstats <- summary(formula, basis=nw)
   for(i in seq(along=sim)){
    sim <- san(formula, meanstats=meanstats, verbose=verbose,
                    burnin=burnin)
    Clist.initial <- ergm.Cprepare(sim, model.initial)
    Clist.miss.initial <- ergm.design(sim, model.initial, initialfit=TRUE,
                                verbose=verbose)
    Clist.initial$meanstats=meanstats
    aaa <- ergm.pl.ihs(Clist=Clist.initial, Clist.miss=Clist.miss.initial,
                     m=model.initial,
                     verbose=verbose)
    pl$zy <- c(pl$zy,aaa$zy)
    pl$foffset <- c(pl$foffset,aaa$foffset)
    pl$xmat <- rbind(pl$xmat,aaa$xmat)
    pl$wend <- c(pl$wend,aaa$wend)
    pl$zy.full <- c(pl$zy.full,aaa$zy.full)
    pl$foffset.full <- c(pl$foffset.full,aaa$foffset.full)
    pl$xmat.full <- rbind(pl$xmat.full,aaa$xmat.full)
    pl$wend.full <- c(pl$wend.full,aaa$wend.full)
   }
  }
  initialfit <- ergm.maple(pl=pl, model.initial,
                          MPLEtype=control$MPLEtype, 
                          verbose=verbose, ...)

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
