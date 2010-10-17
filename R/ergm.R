ergm <- function(formula, response=NULL, theta0="MPLE",
                 MPLEonly=FALSE, MLestimate=!MPLEonly, seed=NULL,
                 burnin=10000, MCMCsamplesize=10000, interval=100,
                 maxit=3,
                 family="PseudoBernoulli.logit",
                 constraints=~.,
                 meanstats=NULL,
                 control=control.ergm(),
                 verbose=FALSE, ...) {
  current.warn <- options()$warn
  options(warn=0)
  if(!is.null(seed))  set.seed(as.integer(seed))
  if (verbose) cat("Evaluating network in model\n")

  nw <- ergm.getnetwork(formula)
  proposalclass <- "c"
  # Next for conditional MLE in dynamic model 
  if(is.character(MLestimate) && 
     (MLestimate=="formation") | (MLestimate=="dissolution")){
   lhs <- terms(formula)[[2]]
   if( is.call(lhs) && (lhs[[1]]=="|")){
   if(MLestimate=="formation"){
    proposalclass <- "fmle"
    y0 <- ergm.getnetwork(as.formula(paste("~",lhs[[3]])))
    y1 <- ergm.getnetwork(as.formula(paste("~",lhs[[2]])))
    #
    # nw contents:
    # y0 y1  nw  nw[2] initial
    #  0  0   0    0    0
    #  0  1   1    0    1
    #  0 NA  NA    0    NA
    #  1  0   1    1    NA
    #  1  1   1    1    NA
    #  1 NA   1    1    NA
    # NA  0  NA    0    NA
    # NA  1   1    0    NA
    # NA NA  NA    0    NA
    #
    nw <- network.copy(y1)
    y0edges <- as.sociomatrix(y0)
    y0edges[is.na(y0edges)] <- 0
    y1edges <- as.sociomatrix(y1)
    y1edges[is.na(y1edges)] <- 0
    ydiscordantedges <- as.matrix(network((!y1edges)&y0edges,
      directed=is.directed(y0)),matrix.type="edgelist")
    add.edges(nw,ydiscordantedges[,1],ydiscordantedges[,2])
    nwm <- network.copy(y1)
    y0edges <- as.sociomatrix(y0)
    y0edges[is.na(y0edges)] <- 1
    y0edges <- as.matrix(network(y0edges,directed=is.directed(y0)),
                         matrix.type="edgelist")
    if(nrow(y0edges)>0){
     for(i in 1:nrow(y0edges)){
      nwm[y0edges[i,1],y0edges[i,2]] <- NA
     }
    }
    MHproposal.miss <- "formationNonObservedMLE"
   }
   if(MLestimate=="dissolution"){
    proposalclass <- "dmle"
    y0 <- ergm.getnetwork(as.formula(paste("~",lhs[[3]])))
    y1 <- ergm.getnetwork(as.formula(paste("~",lhs[[2]])))
    # nw contents:
    # y0 y1  nw  nw[2] initial
    #  0  0   0    0    NA
    #  0  1   0    0    NA
    #  0 NA   0    0    NA
    #  1  0   0    1     0
    #  1  1   1    1     1
    #  1 NA  NA    1    NA
    # NA  0   0    0    NA
    # NA  1  NA    0    NA
    # NA  NA NA    0    NA
    nw <- network.copy(y1)
    y0edges <- as.sociomatrix(y0)
    y0edges[is.na(y0edges)] <- 1
    y1edges <- as.sociomatrix(y1)
    y1edges[is.na(y1edges)] <- 1
    yminusedges <- as.matrix(network((!y0edges)&y1edges,directed=is.directed(y0)),matrix.type="edgelist")
    y1edges <- cbind(unlist(lapply(y1$mel, "[[", "outl")),unlist(lapply(y1$mel, "[[", "inl")))
    y1delete <- match(yminusedges[,1]+yminusedges[,2]*network.size(y1),
                      y1edges[,1]+y1edges[,2]*network.size(y1))
    delete.edges(nw,y1delete)
    #
    y0edges <- as.sociomatrix(is.na(y0))
    y1edges <- as.sociomatrix(y1)
    y1edges[is.na(y1edges)] <- 1
    yminusedges <- as.matrix(network(y0edges&y1edges,directed=is.directed(y0)),matrix.type="edgelist")
    if(nrow(yminusedges)>0){
     for(i in 1:nrow(yminusedges)){
      nw[yminusedges[i,1],yminusedges[i,2]] <- NA
     }
    }
    #
    y0edges <- as.sociomatrix(is.na(y0))
    y1edges <- as.sociomatrix(y1)
    y1edges[is.na(y1edges)] <- 0
#   yminusedges <- as.matrix(network(y0edges&(!y1edges),directed=is.directed(y0)),matrix.type="edgelist")
    yminusedges <- as.matrix(network(y0edges,directed=is.directed(y0)),matrix.type="edgelist")
    set.edge.attribute(y0,attrname="na",value=FALSE)
    nwm <- network.copy(y1)
    y0edges <- as.sociomatrix(y0)
    y0edges[is.na(y0edges)] <- 0
    y0edges <- as.matrix(network(!y0edges,directed=is.directed(y0)),
                         matrix.type="edgelist")
    if(nrow(y0edges)>0){
     for(i in 1:nrow(y0edges)){
      nwm[y0edges[i,1],y0edges[i,2]] <- NA
     }
    }
    if(nrow(yminusedges)>0){
     for(i in 1:nrow(yminusedges)){
      y0[yminusedges[i,1],yminusedges[i,2]] <- 0
     }
    }
    MHproposal.miss <- "dissolutionNonObservedMLE"
   }
   formula.passed<-formula
   formula<-ergm.update.formula(formula,nw~.)
   MLestimate=!MPLEonly
  }}else{
   nwm <- network.copy(nw)
   MHproposal.miss <- "randomtoggleNonObserved"
  }
  # End conditional MLE in dynamic model
  if(!is.null(meanstats)){
   control$drop <- FALSE
   if(!(!is.null(control$SAN.burnin) && is.na(control$SAN.burnin))){
    netsumm<-summary(formula,response=response)
    if(length(netsumm)!=length(meanstats))
      stop("Incorrect length of the meanstats vector: should be ", length(netsumm), " but is ",length(meanstats),".")
    
    if(verbose) cat("Constructing an approximate response network.\n")
    ## If meanstats are given, overwrite the given network and formula
    ## with SAN-ed network and formula.
    nw<-san(formula, meanstats=meanstats,
            theta0=if(is.numeric(theta0)) theta0,
            response=response,
            family=family,
            constraints=constraints,
            verbose=verbose,
            burnin=
            if(is.null(control$SAN.burnin)) burnin
            else control$SAN.burnin,
            interval=interval)
    formula<-ergm.update.formula(formula,nw~.)
    if (verbose) {
     cat("Original meanstats:\n")
     print(meanstats)
     cat("SAN meanstats - Original meanstats:\n")
     print(summary(formula,response=response, basis=nw)-meanstats)
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
  if (verbose) cat("Initializing Metropolis-Hastings proposal.\n")

  MHproposal <- MHproposal(constraints, weights=control$prop.weights, control$prop.args, nw, model.initial,class=proposalclass,family=family)
  # Note:  MHproposal function in CRAN version does not use the "class" argument for now
  MHproposal.miss <- MHproposal(MHproposal.miss, control$prop.args, nw, model.initial)

  conddeg <- switch(MHproposal$name=="CondDegree",control$drop,NULL)
  MCMCparams=c(control,
   list(samplesize=MCMCsamplesize, burnin=burnin, interval=interval,
        maxit=maxit, Clist.miss=NULL, Clist.dt=NULL,
	mcmc.precision=control$mcmc.precision))


  if (verbose) cat("Fitting initial model.\n")
  theta0copy <- theta0
  initialfit <- ergm.initialfit(theta0=theta0copy, MLestimate=MLestimate, 
                                formula=formula, nw=nwm, meanstats=meanstats,
                                m=model.initial,
                                MPLEtype=control$MPLEtype, 
                                initial.loglik=control$initial.loglik,
                                conddeg=conddeg, MCMCparams=MCMCparams, MHproposal=MHproposal,
                                force.MPLE=(!control$force.mcmc && ergm.independencemodel(model.initial)
                                            && constraints==(~.)),
                                verbose=verbose, 
                                compressflag = control$compress, 
                                maxNumDyadTypes=control$maxNumDyadTypes,
                                ...)
  MCMCflag <- ((MLestimate && (!ergm.independencemodel(model.initial)
                               || !is.null(meanstats)
                               || constraints!=(~.)))
                || control$force.mcmc)
  if (MCMCflag) {
    theta0 <- initialfit$coef
    names(theta0) <- model.initial$coef.names
    theta0[is.na(theta0)] <- 0
  } else { # Just return initial (non-MLE) fit and exit.
    initialfit$offset <- model.initial$etamap$offsettheta
    initialfit$drop <- droppedterms
    initialfit$network <- nw
    initialfit$newnetwork <- nw
    initialfit$formula <- formula
    initialfit$constraints <- constraints
    initialfit$prop.args <- control$prop.args
    initialfit$prop.weights <- control$prop.weights
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
  Clist.miss <- ergm.design(nw, model, verbose=FALSE)
  Clist.miss <- ergm.design(nw, model, verbose=FALSE)
  if((MHproposal$name!="FormationMLE")&(MHproposal$name!="DissolutionMLE")){
    Clist.dt <- list(heads=NULL, tails=NULL, nedges=0, dir=is.directed(nw))
  }else{
    Clist.dt <- ergm.Cprepare(y0, model)
  }
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

  MCMCparams=c(control,
   list(samplesize=MCMCsamplesize, burnin=burnin, interval=interval,
        maxit=maxit,Clist.miss=Clist.miss,Clist.dt=Clist.dt, mcmc.precision=control$mcmc.precision))


  if (verbose) cat("Fitting ERGM.\n")
  v <- switch(control$style,
    "Robbins-Monro" = ergm.robmon(theta0, nw, model, Clist, burnin, interval,
                      MHproposal(constraints,weights=control$prop.weights, control$prop.args, nw, model), verbose, control),
    "PILA" = ergm.PILA(theta0, nw, model, Clist,
      MHproposal(constraints,weights=control$prop.weights, control$prop.args, nw, model),
      MCMCparams, control, verbose),
    #  PILA stuff:  Not in CRAN version
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
				MHproposal=MHproposal, MHproposal.miss=MHproposal.miss, 
				verbose=verbose,...),
     ergm.mainfitloop(theta0, nw,
                          model, Clist, 
                          initialfit,
                          MCMCparams=MCMCparams, MHproposal=MHproposal,
                          MHproposal.miss=MHproposal.miss,
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

  v$response=response
  v$family=family

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
  v
}
