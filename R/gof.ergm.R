#  File R/gof.ergm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################

#' Conduct Goodness-of-Fit Diagnostics on a Exponential Family Random Graph
#' Model
#' 
#' \code{\link{gof}} calculates \eqn{p}-values for geodesic distance, degree,
#' and reachability summaries to diagnose the goodness-of-fit of exponential
#' family random graph models.  See \code{\link{ergm}} for more information on
#' these models.
#' 
#' A sample of graphs is randomly drawn from the specified model.  The first
#' argument is typically the output of a call to \code{\link{ergm}} and the
#' model used for that call is the one fit.
#' 
#' For \code{GOF = ~model}, the model's observed sufficient statistics are
#' plotted as quantiles of the simulated sample. In a good fit, the observed
#' statistics should be near the sample median (0.5).
#'
#' By default, the sample consists of 100 simulated networks, but this sample
#' size (and many other settings) can be changed using the \code{control} 
#' argument described above.
#'
#' @aliases gof.default
#' @param object Either a formula or an \code{\link{ergm}} object.
#' See documentation for \code{\link{ergm}}.
#' @param \dots Additional arguments, to be passed to lower-level functions.
#' @param coef When given either a formula or an object of class ergm,
#' \code{coef} are the parameters from which the sample is drawn. By default
#' set to a vector of 0.
#' @param GOF formula; an formula object, of the form \code{~ <model terms>}
#' specifying the statistics to use to diagnosis the goodness-of-fit of the
#' model.  They do not need to be in the model formula specified in
#' \code{formula}, and typically are not.  Currently supported terms are the
#' degree distribution (\dQuote{degree} for undirected graphs, or
#' \dQuote{idegree} and/or \dQuote{odegree} for directed graphs), geodesic
#' distances (\dQuote{distance}), shared partner distributions
#' (\dQuote{espartners} and \dQuote{dspartners}), the triad census
#' (\dQuote{triadcensus}), and the terms of the original model
#' (\dQuote{model}). The default formula for undirected networks is \code{~
#' degree + espartners + distance + model}, and the default formula for
#' directed networks is \code{~ idegree + odegree + espartners + distance +
#' model}.  By default a \dQuote{model} term is added to the formula.  It is a
#' very useful overall validity check and a reminder of the statistical
#' variation in the estimates of the mean value parameters.  To omit the
#' \dQuote{model} term, add \dQuote{- model} to the formula.
#' @param constraints A one-sided formula specifying one or more constraints on
#' the support of the distribution of the networks being modeled. See the help
#' for similarly-named argument in \code{\link{ergm}} for more information. For
#' \code{gof.formula}, defaults to unconstrained. For \code{gof.ergm}, defaults
#' to the constraints with which \code{object} was fitted.
#' @param control A list to control parameters, constructed using
#' \code{\link{control.gof.formula}} or \code{\link{control.gof.ergm}} (which
#' have different defaults).
#' @param unconditional logical; if \code{TRUE}, the simulation is
#' unconditional on the observed dyads.  if not \code{TRUE}, the simulation is
#' conditional on the observed dyads. This is primarily used internally when
#' the network has missing data and a conditional GoF is produced.
#' @param verbose Provide verbose information on the progress of the
#' simulation.
#' @return \code{\link{gof}}, \code{\link{gof.ergm}}, and
#' \code{\link{gof.formula}} return an object of class \code{gof.ergm}, which inherits from class `gof`.  This
#' is a list of the tables of statistics and \eqn{p}-values.  This is typically
#' plotted using \code{\link{plot.gof}}.
#' @seealso [ergm()], [network()], [simulate.ergm()], [summary.ergm()]
#' @keywords models
#' @examples
#' 
#' \donttest{
#' data(florentine)
#' gest <- ergm(flomarriage ~ edges + kstar(2))
#' gest
#' summary(gest)
#' 
#' # test the gof.ergm function
#' gofflo <- gof(gest)
#' gofflo
#' 
#' # Plot all three on the same page
#' # with nice margins
#' par(mfrow=c(1,3))
#' par(oma=c(0.5,2,1,0.5))
#' plot(gofflo)
#' 
#' # And now the log-odds
#' plot(gofflo, plotlogodds=TRUE)
#' 
#' # Use the formula version of gof
#' gofflo2 <-gof(flomarriage ~ edges + kstar(2), coef=c(-1.6339, 0.0049))
#' plot(gofflo2)
#' }
#' 
#' @export gof
gof <- function(object, ...){
      UseMethod("gof")
    }


#' @noRd
#' @importFrom utils methods
#' @export
gof.default <- function(object,...) {
  classes <- setdiff(gsub(pattern="^gof.",replacement="",as.vector(methods("gof"))), "default")
  stop("Goodness-of-Fit methods have been implemented only for class(es) ",
       paste.and(paste('"',classes,'"',sep="")), " in the packages loaded.")
}

#' @describeIn gof Perform simulation to evaluate goodness-of-fit for
#'   a specific [ergm()] fit.
#'
#' @note For \code{gof.ergm} and \code{gof.formula}, default behavior depends on the
#' directedness of the network involved; if undirected then degree, espartners,
#' and distance are used as default properties to examine.  If the network in
#' question is directed, \dQuote{degree} in the above is replaced by idegree
#' and odegree.
#'
#' @export
gof.ergm <- function (object, ..., 
                      coef=NULL,
                      GOF=NULL, 
                      constraints=NULL,
                      control=control.gof.ergm(),
                      verbose=FALSE) {
  check.control.class(c("gof.ergm","gof.formula"), "gof.ergm")
  control.toplevel(...)
  .gof.nw <- as.network(object$network)

  if(!is.null(object$response)) stop("GoF for valued ERGMs is not implemented at this time.")
  
  formula <- nonsimp_update.formula(object$formula, .gof.nw~., from.new=".gof.nw")
# paste("~",paste(unlist(dimnames(attr(terms(formula),"factors"))[-1]),collapse="+"),sep="")
  if(!is.network(.gof.nw)){
    stop("A network must be given as part of the network object.")
  }

  if(is.null(coef)) coef <- coef(object)

  control.transfer <- c("MCMC.burnin", "MCMC.prop.weights", "MCMC.prop.args", "MCMC.packagenames", "MCMC.init.maxedges","term.options")
  for(arg in control.transfer)
    if(is.null(control[[arg]]))
      control[arg] <- list(object$control[[arg]])

  # Rescale the interval by the ratio between the estimation sample size and the GOF sample size so that the total number of MCMC iterations would be about the same.
  NVL(control$MCMC.interval) <- max(ceiling(object$control$MCMC.interval*object$control$MCMC.samplesize/control$nsim),1)

  control <- set.control.class("control.gof.formula")
  
  if(is.null(constraints)) constraints <- object$constraints
  
  gof.formula(object=formula, coef=coef,
              GOF=GOF,
              constraints=constraints,
              control=control,
              verbose=verbose, ...)
}



#' @describeIn gof Perform simulation to evaluate goodness-of-fit for
#'   a model configuration specified by a [`formula`], coefficient,
#'   constraints, and other settings.
#' 
#' @export
gof.formula <- function(object, ..., 
                        coef=NULL,
                        GOF=NULL,
                        constraints=~.,
                        control=NULL,
			unconditional=TRUE,
                        verbose=FALSE) {
  if("response" %in% names(list(...))) stop("GoF for valued ERGMs is not implemented at this time.")

  if(!is.null(control$seed)) {set.seed(as.integer(control$seed))}
  if (verbose) 
    message("Starting GOF for the given ERGM formula.")
  # Unused code
  coefmissing <- NULL
  # get network
  lhs <- ERRVL(try(eval_lhs.formula(object)),
               stop("A network object on the RHS of the formula argument must be given"))
  if(is.ergm(lhs)){
    if(is.null(GOF)) GOF <- nonsimp_update.formula(object, ~.) # Remove LHS from formula.
    if(is.null(constraints)) constraints <- NULL
    if(is.null(control)) control <- control.gof.ergm()
    
    return(gof(lhs, GOF=GOF, coef=coef, control=control, unconditional=unconditional, verbose=verbose, ...)) # Kick it back to gof.ergm.
  }
  
  nw <- as.network(lhs)

  if(is.null(control)) control <- control.gof.formula()

  check.control.class(c("gof.formula","gof.ergm"), "ERGM gof.formula")
  control.toplevel(...)

  #Set up the defaults, if called with GOF==NULL
  if(is.null(GOF)){
    if(is.directed(nw))
      GOF<- ~idegree + odegree + espartners + distance + model
    else
      GOF<- ~degree + espartners + distance + model
  }
  # Add a model term, unless it is explicitly excluded
  GOFtrms <- list_rhs.formula(GOF)
  if(sum(attr(GOFtrms,"sign")[as.character(GOFtrms)=="model"])==0){ # either no "model"s or "-model"s don't outnumber "model"s
    #' @importFrom statnet.common nonsimp_update.formula
      GOF <- nonsimp_update.formula(GOF, ~ . + model)
  }
  
  all.gof.vars <- as.character(list_rhs.formula(GOF))

  # match variables

  for(i in seq(along=all.gof.vars)){
    all.gof.vars[i] <- match.arg(all.gof.vars[i],
                                 c('distance', 'espartners', 'dspartners', 'odegree', 'idegree', 
                                   'degree', 'triadcensus', 'model'
                                   )
                                 )
  }
  GOF <- as.formula(paste("~",paste(all.gof.vars,collapse="+")))
  
  m <- ergm_model(object, nw, term.options=control$term.options)

  proposal <- if(inherits(constraints, "ergm_proposal")) constraints
                else ergm_proposal(constraints,arguments=control$MCMC.prop.args,
                                   nw=nw, weights=control$MCMC.prop.weights, class="c"## ,reference=reference,response=response
                                   )

  if(is.null(coef)){
      coef <- numeric(nparam(m, canonical=FALSE))
      warning("No parameter values given, using 0.")
  }
# if(is.bipartite(nw)){
#     coef <- c(coef,-1)
# }

  # If missing simulate from the conditional model
  if(network.naedgecount(nw) & unconditional){
   if(verbose){message("Conditional simulations for missing fit")}
   if(is.null(coefmissing)){coefmissing <- coef}
   constraints.obs<-nonsimp_update.formula(constraints,~.+observed)
   SimCond <- gof(object=object, coef=coefmissing,
                  GOF=GOF, 
                  constraints=constraints.obs,
                  control=control,
                  unconditional=FALSE,
                  verbose=verbose)
  }

# test to see which of these is/are necessary
#  pval.model<-pval.triadcensus<-pval.dist<-pval.deg<-pval.espart<-pval.espart<-NULL
##
#  obs.model<-pobs.model<-sim.model<-psim.model<-pval.model<-bds.model<-NULL
#  obs.triadcensus<-pobs.triadcensus<-sim.triadcensus<-psim.triadcensus<-pval.triadcensus<-bds.triadcensus<-NULL
#  obs.dist<-pobs.dist<-sim.dist<-psim.dist<-pval.dist<-bds.dist<-NULL
#  obs.deg<-pobs.deg<-sim.deg<-psim.deg<-pval.deg<-bds.deg<-NULL
#  obs.espart<-pobs.espart<-sim.espart<-psim.espart<-pval.espart<-bds.espart<-NULL
#  obs.dspart<-pobs.dspart<-sim.dspart<-psim.dspart<-pval.dspart<-bds.dspart<-NULL
#
#  obs.ideg<-pobs.ideg<-sim.ideg<-psim.ideg<-pval.ideg<-bds.ideg<-pval.ideg<-NULL
#  obs.odeg<-pobs.odeg<-sim.odeg<-psim.odeg<-pval.odeg<-bds.odeg<-pval.odeg<-NULL

  n <- network.size(nw)

  # Calculate network statistics for the observed graph
  # Set up the output arrays of sim variables
  if(verbose)
    message("Calculating observed network statistics.")
  
  if ('model' %in% all.gof.vars) {
   if(!network.naedgecount(nw) | !unconditional){
    obs.model <- summary(object, term.options=control$term.options)
   }else{
    obs.model <- SimCond$obs.model
   }
   sim.model <- array(0,dim=c(control$nsim,length(obs.model)))
   dimnames(sim.model) <- list(paste(c(1:control$nsim)),names(obs.model))
  }

  if ('distance' %in% all.gof.vars) {
   if(!network.naedgecount(nw) | !unconditional){
    obs.dist <- ergm.geodistdist(nw)
    obs.dist[obs.dist==Inf] <- n
   }else{
    obs.dist <- SimCond$summary.dist[,"mean"]
   }
   sim.dist <-array(0,dim=c(control$nsim,n))
   dimnames(sim.dist)  <- list(paste(c(1:control$nsim)),paste(1:n))
  }

  if ('odegree' %in% all.gof.vars) {
   if(!network.naedgecount(nw) | !unconditional){
    mesp <- paste("c(",paste(0:(n-1),collapse=","),")",sep="")
    obs.odeg <- summary(as.formula(paste('nw ~ odegree(',mesp,')',sep="")))
   }else{
    obs.odeg <- SimCond$summary.odeg[,"mean"]
   }
   sim.odeg <- array(0,dim=c(control$nsim,n))
#  obs.odeg <- c(obs.odeg,rep(0,n-length(obs.odeg)))
   dimnames(sim.odeg)   <- list(paste(c(1:control$nsim)),paste(0:(n-1)))
   names(obs.odeg) <- dimnames(sim.odeg)[[2]]
  }

  if ('idegree' %in% all.gof.vars) {
   if(!network.naedgecount(nw) | !unconditional){
    mesp <- paste("c(",paste(0:(n-1),collapse=","),")",sep="")
    obs.ideg <- summary(as.formula(paste('nw ~ idegree(',mesp,')',sep="")))
   }else{
    obs.ideg <- SimCond$summary.ideg[,"mean"]
   }
   sim.ideg <- array(0,dim=c(control$nsim,n))
#  obs.ideg <- c(obs.ideg,rep(0,n-length(obs.ideg)))
   dimnames(sim.ideg)   <- list(paste(c(1:control$nsim)),paste(0:(n-1)))
   names(obs.ideg) <- dimnames(sim.ideg)[[2]]
  }

  if ('degree' %in% all.gof.vars) {
   if(!network.naedgecount(nw) | !unconditional){
    if(is.bipartite(nw)){
     obs.deg <- degreedist(nw, print=FALSE)$b2
     obs.deg <- c(obs.deg,rep(0,n-length(obs.deg)))
    }else{
     mesp <- paste("c(",paste(0:(n-1),collapse=","),")",sep="")
     obs.deg <- summary(as.formula(paste('nw ~ degree(',mesp,')',sep="")))
    }
   }else{
    obs.deg <- SimCond$summary.deg[,"mean"]
   }
   sim.deg <- array(0,dim=c(control$nsim,n))
   dimnames(sim.deg)   <- list(paste(c(1:control$nsim)),paste(0:(n-1)))
   names(obs.deg) <- dimnames(sim.deg)[[2]]
  }
 
  if ('espartners' %in% all.gof.vars) {
#  obs.espart <- espartnerdist(nw, print=verbose)
   if(!network.naedgecount(nw) | !unconditional){
    mesp <- paste("c(",paste(0:(network.size(nw)-2),collapse=","),")",sep="")
    obs.espart <- summary(as.formula(paste('nw ~ esp(',mesp,')',sep="")))
   }else{
    obs.espart <- SimCond$summary.espart[,"mean"]
   }
   sim.espart <- array(0,dim=c(control$nsim,n-1))
   dimnames(sim.espart) <- list(paste(c(1:control$nsim)),paste(0:(n-2)))
  }
 
  if ('dspartners' %in% all.gof.vars) {
   if(!network.naedgecount(nw) | !unconditional){
#   obs.dspart <- dspartnerdist(nw, print=verbose)
    mesp <- paste("c(",paste(0:(network.size(nw)-2),collapse=","),")",sep="")
    obs.dspart <- summary(as.formula(paste('nw ~ dsp(',mesp,')',sep="")))
   }else{
    obs.dspart <- SimCond$summary.dspart[,"mean"]
   }
   sim.dspart <- array(0,dim=c(control$nsim,n-1))
   dimnames(sim.dspart) <- list(paste(c(1:control$nsim)),paste(0:(n-2)))
  }

  if ('triadcensus' %in% all.gof.vars) {
   if(is.directed(nw)){
    triadcensus <- 0:15
    namestriadcensus <- c("003","012", "102", "021D", "021U", "021C",
      "111D", "111U", "030T",
      "030C", "201", "120D", "120U", "120C", "210", "300")
    triadcensus.formula <- "~ triadcensus(0:15)"
   }else{
    triadcensus <- 0:3
    namestriadcensus <- c("0","1","2", "3")
    triadcensus.formula <- "~ triadcensus(0:3)"
   }
   if(!network.naedgecount(nw) | !unconditional){
    obs.triadcensus <- summary(as.formula(paste('nw',triadcensus.formula,sep="")))
   }else{
    obs.triadcensus <- SimCond$summary.triadcensus[,"mean"]
   }
   sim.triadcensus <- array(0,dim=c(control$nsim,length(triadcensus)))
   dimnames(sim.triadcensus) <- list(paste(c(1:control$nsim)), namestriadcensus)
   names(obs.triadcensus) <- namestriadcensus
  }
 
  # Simulate an exponential family random graph model

#  SimNetworkSeriesObj <- simulate(object, control$nsim=control$nsim, seed=seed,
#                                  coef=coef,
#                                  burnin=burnin, interval=interval,
#                                  constraints=constraints,
#                                  control=control.simulate.formula(
#                                   prop.args=control$MCMC.prop.args,
#                                   prop.weights=control$MCMC.prop.weights,
#                                   summarizestats=control$summarizestats,
#                                   drop=control$drop),
#                                  verbose=verbose, basis=nw)
# New approach below avoids having to store gigantic unnecessary
# network.list object

  if(verbose)
    message("Starting simulations.")

  tempnet <- nw
  for (i in 1:control$nsim) {
    if(verbose){
      message("Sim ",i," of ",control$nsim,": ",appendLF=FALSE)
    }
    if(network.naedgecount(nw) & !unconditional){tempnet <- nw}
    tempnet <- simulate(m, nsim=1, coef=coef,
                        constraints=proposal,
                        control=set.control.class("control.simulate.formula",control),
                        basis=tempnet,
                        verbose=verbose)
    seed <- NULL # Don't re-seed after first iteration   
#    if(verbose){
#     cat(paste("...",i,sep=""))
#     if ((i %% 10 == 0) || (i==control$nsim)) cat("\n")
#    }
    if ('model' %in% all.gof.vars) {
     sim.model[i,] <- summary(nonsimp_update.formula(object,tempnet ~ ., from.new="tempnet"), term.options=control$term.options)
    }
    if ('distance' %in% all.gof.vars) {
     sim.dist[i,] <- ergm.geodistdist(tempnet)
    }
    if ('idegree' %in% all.gof.vars) {
     mesp <- paste("c(",paste(0:(n-1),collapse=","),")",sep="")
     gi <- tempnet
     sim.ideg[i,] <- summary(as.formula(paste('gi ~ idegree(',mesp,')',sep="")))
#    temp <- table(degreedist(tempnet, print=verbose)[1,])
#    sim.ideg[i,] <- c(temp, rep(0, n-length(temp)))
    }
    if ('odegree' %in% all.gof.vars) {
     mesp <- paste("c(",paste(0:(n-1),collapse=","),")",sep="")
     gi <- tempnet
     sim.odeg[i,] <- summary(as.formula(paste('gi ~ odegree(',mesp,')',sep="")))
#    temp <- table(degreedist(tempnet, print=verbose)[2,])
#    sim.odeg[i,] <- c(temp, rep(0, n-length(temp)))
    }
    if ('degree' %in% all.gof.vars) {
     gi <- tempnet
     if(is.bipartite(gi)){
      temp <- degreedist(gi, print=FALSE)$b2
      sim.deg[i,] <- c(temp,rep(0,n-length(temp)))
     }else{                                                
      mesp <- paste("c(",paste(0:(n-1),collapse=","),")",sep="")
      sim.deg[i,] <- summary(as.formula(paste('gi ~ degree(',mesp,')',sep="")))
     }
#    temp <- table(degreedist(tempnet, print=verbose))
#    sim.deg[i,] <- c(temp, rep(0, n-length(temp)))
    }
    if ('espartners' %in% all.gof.vars) {
#    sim.espart[i,] <- espartnerdist(tempnet,
#                                   print=verbose)
     gi <- tempnet
     mesp <- paste("c(",paste(0:(network.size(gi)-2),collapse=","),")",sep="")
     sim.espart[i,] <- summary(as.formula(paste('gi ~ esp(',mesp,')',sep="")))
    }
    if ('dspartners' %in% all.gof.vars) {
#    sim.espart[i,] <- dspartnerdist(tempnet,
#                                   print=verbose)
     gi <- tempnet
     mesp <- paste("c(",paste(0:(network.size(gi)-2),collapse=","),")",sep="")
     sim.dspart[i,] <- summary(as.formula(paste('gi ~ dsp(',mesp,')',sep="")))
    }
    if ('triadcensus' %in% all.gof.vars) {
     gi <- tempnet
     sim.triadcensus[i,] <- summary(as.formula(paste('gi',triadcensus.formula,sep="")))
    }
  }
  if(verbose){
    message("")
  }

  # calculate p-values
  
  returnlist <- list(network.size=n, GOF=GOF)

  calc_pvals <- function(gv, count=TRUE){
    sim <- get(paste("sim", gv, sep="."))
    obs <- get(paste("obs", gv, sep="."))

    pval <- apply(sim <= obs[col(sim)],2,mean)
    pval.top <- apply(sim >= obs[col(sim)],2,mean)
    pval <- cbind(obs,apply(sim, 2,min), apply(sim, 2,mean),
                  apply(sim, 2,max), pmin(1,2*pmin(pval,pval.top)))
    dimnames(pval)[[2]] <- c("obs","min","mean","max","MC p-value")
    pobs <- if(count) obs/sum(obs) else pval.top
    if(count){
      psim <- sweep(sim,1,apply(sim,1,sum),"/")
      psim[is.na(psim)] <- 1
    }else{
      psim <- apply(sim,2,rank)/nrow(sim)
      psim <- matrix(psim, ncol=ncol(sim)) # Guard against the case of sim having only one row.
    }
    bds <- apply(psim,2,quantile,probs=c(0.025,0.975))

    l <- list(pval = pval, summary = pval, pobs = pobs, psim = psim, bds = bds, obs = obs, sim = sim)
    setNames(l, paste(names(l), gv, sep="."))
  }

  GVMAP = list(model=c('model', FALSE),
               distance=c('dist', TRUE),
               idegree=c('ideg', TRUE),
               odegree=c('odeg', TRUE),
               degree=c('deg', TRUE),
               espartners=c('espart', TRUE),
               dspartners=c('dspart', TRUE),
               triadcensus=c('triadcensus', TRUE))
  for(gv in names(GVMAP))
    if (gv %in% all.gof.vars)
      returnlist <- modifyList(returnlist, calc_pvals(GVMAP[[gv]][1],GVMAP[[gv]][2]))

  class(returnlist) <- c("gof.ergm", "gof")
  returnlist
}



################################################################
# The <print.gof> function prints the summary matrices
# of each GOF term included in the build of the gof
#
# --PARAMETERS--
#   x  : a gof, as returned by one of the <gof.X> functions
#   ...: additional printing parameters; these are ignored
#
# --RETURNED--
#   NULL
#################################################################

#' @describeIn gof
#' 
#' \code{\link{print.gof}} summaries the diagnostics such as the
#' degree distribution, geodesic distances, shared partner
#' distributions, and reachability for the goodness-of-fit of
#' exponential family random graph models.  See \code{\link{ergm}} for
#' more information on these models. (`summary.gof` is a deprecated
#' alias that may be repurposed in the future.)
#' 
#' @param x an object of class \code{gof} for printing or plotting.
#' @export
print.gof <- function(x, ...){
  all.gof.vars <- as.character(list_rhs.formula(x$GOF))
  # match variables
  goftypes <- matrix( c(
      "model", "model statistics", "summary.model",
      "distance", "minimum geodesic distance", "summary.dist",
      "idegree", "in-degree", "summary.ideg",
      "odegree", "out-degree", "summary.odeg",
      "degree", "degree", "summary.deg",
      "espartners", "edgewise shared partner", "summary.espart",
      "dspartners", "dyadwise shared partner", "summary.dspart",
      "triadcensus", "triad census", "summary.triadcensus"), 
                      byrow=TRUE, ncol=3)
  for(i in seq(along=all.gof.vars)){
    all.gof.vars[i] <- match.arg(all.gof.vars[i], goftypes[,1])
  }
  for(statname in all.gof.vars){
    r <- match(statname, goftypes[,1])  # find row in goftypes matrix
    cat("\nGoodness-of-fit for", goftypes[r, 2],"\n\n")
    m <- x[[goftypes[r, 3] ]] # get summary statistics
    zerorows <- m[,"obs"]==0 & m[,"min"]==0 & m[,"max"]==0
    print(m[!zerorows,])
  }
  invisible()
}

###################################################################
# The <plot.gof> function plots the GOF diagnostics for each
# term included in the build of the gof
#
# --PARAMETERS--
#   x          : a gof, as returned by one of the <gof.X>
#                functions
#   ...        : additional par arguments to send to the native R
#                plotting functions
#   cex.axis   : the magnification of the text used in axis notation;
#                default=0.7
#   plotlogodds: whether the summary results should be presented
#                as their logodds; default=FALSE
#   main       : the main title; default="Goodness-of-fit diagnostics"
#   verbose    : this parameter is ignored; default=FALSE
#   normalize.reachibility: whether to normalize the distances in
#                the 'distance' GOF summary; default=FALSE
#
# --RETURNED--
#   NULL
#
###################################################################



#' @describeIn gof
#' 
#' \code{\link{plot.gof}} plots diagnostics such as the degree
#' distribution, geodesic distances, shared partner distributions, and
#' reachability for the goodness-of-fit of exponential family random graph
#' models.  See \code{\link{ergm}} for more information on these models.
#' 
#' @param cex.axis Character expansion of the axis labels relative to that for
#' the plot.
#' @param plotlogodds Plot the odds of a dyad having given characteristics
#' (e.g., reachability, minimum geodesic distance, shared partners). This is an
#' alternative to the probability of a dyad having the same property.
#' @param main Title for the goodness-of-fit plots.
#' @param normalize.reachability Should the reachability proportion be
#' normalized to make it more comparable with the other geodesic distance
#' proportions.
#' @keywords graphs
#' 
#' @importFrom graphics boxplot points lines mtext plot
#' @export
plot.gof <- function(x, ..., 
         cex.axis=0.7, plotlogodds=FALSE,
         main="Goodness-of-fit diagnostics", 
         normalize.reachability=FALSE,
         verbose=FALSE) {

 color <- "gray75"
#par(oma=c(0.5,2,1,0.5))

#statsno <- (sum(stats=='deg')>0) + (sum(stats=='espart')>0) + (sum(stats=='d
 all.gof.vars <- as.character(list_rhs.formula(x$GOF))
 statsno <- length(all.gof.vars)

# match variables

 for(i in seq(along=all.gof.vars)){
   all.gof.vars[i] <- match.arg(all.gof.vars[i],
    c('distance', 'triadcensus', 'espartners', 'dspartners', 'odegree', 'idegree', 
      'degree', 'model'
     )
                               )
 }
 GOF <- as.formula(paste("~",paste(all.gof.vars,collapse="+")))

 if(statsno==0){
  stop("The gof object does not contain any statistics!\n")
 }
 n <- x$network.size

#attach(x) 
  
  gofcomp <- function(tag, unit, idx=c("finite","infinite","nominal")){
    idx <- match.arg(idx)

    pval <- x[[paste0("pval.",tag)]]
    sim <- x[[paste0("sim.",tag)]]
    obs <- x[[paste0("obs.",tag)]]
    psim <- x[[paste0("psim.",tag)]]
    pobs <- x[[paste0("pobs.",tag)]]
    bds <- x[[paste0("bds.",tag)]]
    
    if(idx == "nominal"){
      ngs <- nrow(pval)
      i <- seq_len(ngs)
    }else{
      ngs <- min(n-1, nrow(pval))
      
      if( min(pval[,"MC p-value"]) <1) {
        pval.max <- max((1:ngs)[pval[1:ngs, "MC p-value"] < 1]) + 3
      }
      else {
        pval.max <- max((1:ngs)[obs[1:ngs] > 0]) + 3
      }
      
      i <- seq_len(if(is.finite(pval.max) & pval.max < n) pval.max
                   else ngs+1)
    }
    if(idx == "infinite"){
      i <- c(i, NA, n)
    }

    if (plotlogodds) {
      odds <- psim
      odds[!is.na(odds)] <- log(odds[!is.na(odds)]/(1 - odds[!is.na(odds)]))
      odds.obs <- pobs
      odds.obs[!is.na(odds.obs)] <- log(odds.obs[!is.na(odds.obs)]/(1 - odds.obs[!is.na(odds.obs)]))
      odds.bds <- bds
      odds.bds[!is.na(odds.bds)] <- log(odds.bds[!is.na(odds.bds)]/(1 - odds.bds[!is.na(odds.bds)]))
      mininf <- min(min(odds[is.finite(odds)]),min(odds.obs[is.finite(odds.obs)]),min(odds.bds[is.finite(odds.bds)]))
      maxinf <- max(max(odds[is.finite(odds)]),max(odds.obs[is.finite(odds.obs)]),max(odds.bds[is.finite(odds.bds)]))

      odds[!is.finite(odds)&odds>0] <- maxinf
      odds[!is.finite(odds)&odds<0] <- mininf
      odds.obs[!is.finite(odds.obs)&odds.obs>0] <- maxinf
      odds.obs[!is.finite(odds.obs)&odds.obs<0] <- mininf
      odds.bds[!is.finite(odds.bds)&odds.bds>0] <- maxinf
      odds.bds[!is.finite(odds.bds)&odds.bds<0] <- mininf

      out <- odds
      out.obs <- odds.obs
      out.bds <- odds.bds

      ylab <- paste0("log-odds for a ", unit)
    }
    else {
      out <- psim
      out.obs <- pobs
      out.bds <- bds
      ylab <- paste0("proportion of ", unit, "s")
    }
    pnames <- NVL(colnames(sim), i-1)[i]

    list(out=out, pnames=pnames, out.obs=out.obs, out.bds=out.bds, i=i, ylab=ylab)
  }

  gofplot <- function(gofstats, xlab){
    out <- gofstats$out
    out.obs <- gofstats$out.obs
    out.bds <- gofstats$out.bds
    i <- gofstats$i
    ylab <- gofstats$ylab
    pnames <- gofstats$pnames

    ymin <- min(min(out,na.rm=TRUE),min(out.obs,na.rm=TRUE))
    ymax <- max(max(out,na.rm=TRUE),max(out.obs,na.rm=TRUE))

    boxplot(data.frame(out[, na.omit(i), drop=FALSE]), xlab = xlab,
            ylab = ylab, names = na.omit(pnames), cex.axis = cex.axis, outline=FALSE,
            ylim=c(ymin,ymax), ...)

    points(cumsum(!is.na(i)), out.bds[1,i], pch = 1,cex=0.75)
    points(cumsum(!is.na(i)), out.bds[2,i], pch = 1,cex=0.75)
    lines(cumsum(!is.na(i)), out.bds[1,i], pch = 18,lty=1,lwd=1,col=color)
    lines(cumsum(!is.na(i)), out.bds[2,i], pch = 18,lty=1,lwd=1,col=color)
    points(cumsum(!is.na(i)), out.obs[i], pch = 16,cex=0.75)
    lines(cumsum(!is.na(i)), out.obs[i], lty = 1,lwd=3)
    points(cumsum(!is.na(i)), colMeans(out[, i, drop=FALSE]), pch=18, cex=2, col="blue")
  }

 ###model####

 for(statname in all.gof.vars){

  if ('model' == statname) {
    gofcomp("model", "statistic", "nominal") %>% gofplot("model statistics")
  }

  if ('degree' == statname) {
    gofcomp("deg", "node") %>% gofplot("degree")
  }

  if ('odegree' == statname) {
    gofcomp("odeg", "node") %>% gofplot("out degree")
  }

  if ('idegree' == statname) {
    gofcomp("ideg", "node") %>% gofplot("in degree")
  }

  if ('espartners' == statname) {
    gofcomp("espart", "edge") %>% gofplot("edge-wise shared partners")
  }

  if ('dspartners' == statname) {
    gofcomp("dspart", "dyad") %>% gofplot("dyad-wise shared partners")
  }

  if ('triadcensus' == statname) {
    gofcomp("triadcensus", "triad", "nominal") %>% gofplot("triad census")
  }

  if ('distance' == statname) {
    gc <- gofcomp("dist", "dyad", "infinite")
    ult(gc$pnames) <- "NR"
    if(normalize.reachability){
      gc <- within(gc,
      {
        mi <- max(i,na.rm=TRUE)
        totrange <- range(out.bds[1,][out.bds[1,] > out.bds[1,mi]],
                          out.bds[2,][out.bds[2,] < out.bds[2,mi]])
        out[,mi] <- (out[,mi]-out.bds[1,mi]) * 
          diff(totrange) / diff(out.bds[,mi]) + totrange[1]
        out.obs[mi] <- (out.obs[mi]- out.bds[1,mi]) *
          diff(totrange) / diff(out.bds[,mi]) + totrange[1]
        out.bds[,mi] <- totrange
      }
      )
    }
    gofplot(gc, "minimum geodesic distance")
  }
 }
   mtext(main,side=3,outer=TRUE,cex=1.5,padj=2)
   invisible()
}





