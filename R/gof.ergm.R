#  File R/gof.ergm.R in package ergm, part of the Statnet suite of packages for
#  network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################

#' Conduct Goodness-of-Fit Diagnostics on a Exponential Family Random Graph
#' Model
#' 
#' [gof()] calculates \eqn{p}-values for geodesic distance, degree,
#' and reachability summaries to diagnose the goodness-of-fit of exponential
#' family random graph models.  See [ergm()] for more information on
#' these models.
#' 
#' A sample of graphs is randomly drawn from the specified model.  The first
#' argument is typically the output of a call to [ergm()] and the
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
#' @param object Either a formula or an [`ergm`] object.
#' See documentation for [ergm()].
#' @param \dots Additional arguments, to be passed to lower-level functions.
#' @param coef When given either a formula or an object of class ergm,
#' \code{coef} are the parameters from which the sample is drawn. By default
#' set to a vector of 0. \matchnames{coefficient}
#' @param GOF formula; an formula object, of the form \code{~ <model terms>}
#' specifying the statistics to use to diagnosis the goodness-of-fit of the
#' model.  They do not need to be in the model formula specified in
#' \code{formula}, and typically are not.  Currently supported terms are the
#' degree distribution (\dQuote{degree} for undirected graphs,
#' \dQuote{idegree} and/or \dQuote{odegree} for directed graphs, and \dQuote{b1degree} and \dQuote{b2degree} for bipartite undirected graphs), geodesic
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
#' for similarly-named argument in [ergm()] for more information. For
#' \code{gof.formula}, defaults to unconstrained. For \code{gof.ergm}, defaults
#' to the constraints with which \code{object} was fitted.
#'
#' @templateVar mycontrols [control.gof.formula()] or [control.gof.ergm()]
#' @template control2
#' @template verbose
#'
#' @param unconditional logical; if \code{TRUE}, the simulation is
#' unconditional on the observed dyads.  if not \code{TRUE}, the simulation is
#' conditional on the observed dyads. This is primarily used internally when
#' the network has missing data and a conditional GoF is produced.
#' @template basis
#' @return [gof()], [gof.ergm()], and
#' [gof.formula()] return an object of class \code{gof.ergm}, which inherits from class `gof`.  This
#' is a list of the tables of statistics and \eqn{p}-values.  This is typically
#' plotted using [plot.gof()].
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
                      coef = coefficients(object),
                      GOF = NULL,
                      constraints = object$constraints,
                      control = control.gof.ergm(),
                      verbose = FALSE) {
  check_dots_used(error = unused_dots_warning)
  check.control.class(c("gof.ergm","gof.formula"), "gof.ergm")
  handle.control.toplevel("gof.ergm", ...)

  if(is.valued(object)) stop("GoF for valued ERGMs is not implemented at this time.")
  
  # If both the passed control and the object's control are NULL (such as if MPLE was estimated), overwrite with simulate.formula()'s defaults.
  formula.control <- control.simulate.formula()
  for(arg in STATIC_MCMC_CONTROLS)
    if(is.null(control[[arg]]))
      control[arg] <- list(NVL(object$control[[arg]], formula.control[[arg]]))

  MCMC.interval.set <- !is.null(control$MCMC.interval)
  for(arg in SCALABLE_MCMC_CONTROLS)
    if(is.null(control[[arg]]))
      control[arg] <- list(EVL(object$control[[arg]]*control$MCMC.scale, formula.control[[arg]]))

  # Rescale the interval by the ratio between the estimation sample size and the GOF sample size so that the total number of MCMC iterations would be about the same.
  if(!MCMC.interval.set) control$MCMC.interval <- max(ceiling(control$MCMC.interval*EVL(object$control$MCMC.samplesize/control$nsim,1)),1)

  control <- set.control.class("control.gof.formula")

  gof.formula(object=object$formula, coef=coef,
              GOF=GOF,
              constraints=constraints,
              control=control,
              basis=object$network,
              verbose=verbose, ...)
}

GOF_VALID_VARS <- c('distance', 'espartners', 'dspartners', 'odegree', 'idegree',
                    'degree', 'triadcensus', 'model', 'b1degree', 'b2degree')

#' @describeIn gof Perform simulation to evaluate goodness-of-fit for
#'   a model configuration specified by a [`formula`], coefficient,
#'   constraints, and other settings.
#' 
#' @export
gof.formula <- function(object, ..., 
                        coef=NULL,
                        GOF=NULL,
                        constraints=~.,
                        basis=eval_lhs.formula(object),
                        control=NULL,
                        unconditional=TRUE,
                        verbose=FALSE) {
  if("response" %in% ...names()) stop("GoF for valued ERGMs is not implemented at this time.")

  if(!is.null(control$seed)){
    set.seed(as.integer(control$seed))
  }
  if (verbose) message("Starting GOF for the given ERGM formula.")

  if(is.ergm(basis)){ # Kick it back to gof.ergm().
    NVL(GOF) <- nonsimp_update.formula(object, ~.) # Remove LHS from formula.
    NVL(control) <- control.gof.ergm()
    
    return(
      gof(basis, GOF = GOF, coef = coef, control = control, unconditional = unconditional, verbose = verbose, ...)
    )
  }

  # Otherwise, LHS/basis must be a network.
  nw <- ensure_network(basis)
  NVL(control) <- control.gof.formula()

  check_dots_used(error = unused_dots_warning)
  check.control.class(c("gof.formula","gof.ergm"), "ERGM gof.formula")
  handle.control.toplevel("gof.formula", ...)

  #Set up the defaults, if called with GOF==NULL
  if(is.null(GOF)){
    GOF <-
      if(is.directed(nw)) ~idegree + odegree + espartners + distance + model
      else if(is.bipartite(nw)) ~b1degree + b2degree + espartners + distance + model
      else ~degree + espartners + distance + model
  }

  # Add a model term, unless it is explicitly excluded
  GOFtrms <- list_rhs.formula(GOF)
  if(sum(attr(GOFtrms,"sign")[as.character(GOFtrms)=="model"])==0){ # either no "model"s or "-model"s don't outnumber "model"s
    GOFtrms <- c(GOFtrms[as.character(GOFtrms)!="model"], list(as.name("model")))
  }

  # match variables
  all.gof.vars <- as.character(GOFtrms[attr(GOFtrms,"sign")>0]) %>%
    sapply(match.arg, GOF_VALID_VARS)

  GOF <- as.formula(paste("~",paste(all.gof.vars,collapse="+")), baseenv())

  m <- ergm_model(object, nw, term.options=control$term.options)

  proposal <- if(inherits(constraints, "ergm_proposal")) constraints
                else ergm_proposal(constraints,arguments=control$MCMC.prop.args,
                                   nw=nw, weights=control$MCMC.prop.weights, class="c", term.options=control$term.options## ,reference=reference
                                   )

  if(is.null(coef)){
      coef <- numeric(nparam(m, canonical=FALSE))
      warning("No parameter values given, using 0.")
  }

  # If missing simulate from the conditional model
  if(network.naedgecount(nw) & unconditional){
   if(verbose){message("Conditional simulations for missing fit")}
   constraints.obs<-nonsimp_update.formula(constraints,~.+observed)
   SimCond <- gof(object=object, coef=coef,
                  GOF=GOF, 
                  constraints=constraints.obs,
                  control=control,
                  basis=basis,
                  unconditional=FALSE,
                  verbose=verbose, ...)
  }

  n <- network.size(nw)

  # Calculate network statistics for the observed graph
  # Set up the output arrays of sim variables
  if(verbose)
    message("Calculating observed network statistics.")

  summ_form <- function(term, range){
    as.formula(call('~',call(term, range)))
  }

  if(is.bipartite(nw)){
    nb1 <- nw %n% "bipartite"
    nb2 <- n-nb1
  }else{
    nb1 <- nb2 <- n
  }

  if(is.directed(nw)){
    triadcensus <- 0:15
    namestriadcensus <- c("003","012", "102", "021D", "021U", "021C",
                          "111D", "111U", "030T",
                          "030C", "201", "120D", "120U", "120C", "210", "300")
  }else{
    triadcensus <- 0:3
    namestriadcensus <- c("0","1","2", "3")
  }

  GVMAP <- list(model=list('model', NULL, object),
                distance=list('dist', 1:n, function(x){o <- ergm.geodistdist(x); o[o==Inf]<-n; o}),
                odegree=list('odeg', 0:(n-1), summ_form('odegree', 0:(n-1))),
                idegree=list('ideg', 0:(n-1), summ_form('idegree', 0:(n-1))),
                degree=list('deg', 0:(n-1), summ_form('degree', 0:(n-1))),
                b1degree=list('b1deg', 0:nb2, summ_form('b1degree', 0:nb2)),
                b2degree=list('b2deg', 0:nb1, summ_form('b2degree', 0:nb1)),
                espartners=list('espart', 0:(n-2), summ_form('esp', 0:(n-2))),
                dspartners=list('dspart', 0:(n-2), summ_form('dsp', 0:(n-2))),
                triadcensus=list('triadcensus', namestriadcensus, summ_form('triadcensus', triadcensus)))

  GVMAP <- GVMAP[names(GVMAP)%in%all.gof.vars]

  # If gv[[3]] is a formula, preinitialize the model and have the
  # function execute it on the network.
  for(i in seq_along(GVMAP))
    if(is(GVMAP[[i]][[3]], "formula"))
      GVMAP[[i]][[3]] <- local({
        m <- ergm_model(GVMAP[[i]][[3]], nw, term.options = control$term.options)
        function(x) summary(m, nw = x)
      })

  calc_obs_stat <- function(gv, names, calc){
    simname <- paste("sim", gv, sep=".")
    obsname <- paste("obs", gv, sep=".")

    obs <- if(!network.naedgecount(nw) | !unconditional) calc(nw)
           else colMeans(SimCond[[simname]])
    assign(obsname, obs, parent.frame())

    sim <- array(0,
                 dim = c(control$nsim,length(obs)),
                 dimnames = list(paste(c(1:control$nsim)), NVL3(names, paste(.), names(obs))))
    assign(simname, sim, parent.frame())
  }

  for(gv in GVMAP)
    calc_obs_stat(gv[[1]], gv[[2]], gv[[3]])

  if(verbose)
    message("Starting simulations.")

  myenv <- environment()

  calc_sim_stat <- function(nw, gv, calc, i){
    simname <- paste("sim", gv, sep=".")
    sim <- get(simname)
    sim[i,] <- calc(nw)
    assign(simname, sim, myenv)
  }

  summfun <- function(state, iter, ...)
    for(gv in GVMAP) calc_sim_stat(as.network(state), gv[[1]], gv[[3]], iter)

  simulate(m, nsim=control$nsim, coef=coef,
           constraints=proposal,
           control=set.control.class("control.simulate.formula",control),
           output=summfun,
           basis=nw,
           verbose=verbose, ...)

  # calculate p-values
  
  returnlist <- list(network.size=n, GOF=GOF)

  calc_pvals <- function(gv, names){
    sim <- get(paste("sim", gv, sep="."))
    obs <- get(paste("obs", gv, sep="."))

    pval <- apply(sim <= obs[col(sim)],2,mean)
    pval.top <- apply(sim >= obs[col(sim)],2,mean)
    pval <- cbind(obs,apply(sim, 2,min), apply(sim, 2,mean),
                  apply(sim, 2,max), pmin(1,2*pmin(pval,pval.top)))
    dimnames(pval)[[2]] <- c("obs","min","mean","max","MC p-value")
    pobs <- if(!is.null(names)) obs/sum(obs) else pval.top
    if(!is.null(names)){
      psim <- sweep(sim,1,apply(sim,1,sum),"/")
      psim %[f]% is.na <- 1
    }else{
      psim <- apply(sim,2,rank)/nrow(sim)
      psim <- matrix(psim, ncol=ncol(sim)) # Guard against the case of sim having only one row.
    }
    bds <- apply(psim,2,quantile,probs=c(0.025,0.975))

    l <- list(pval = pval, summary = pval, pobs = pobs, psim = psim, bds = bds, obs = obs, sim = sim)
    setNames(l, paste(names(l), gv, sep="."))
  }

  for(gv in GVMAP)
    returnlist <- modifyList(returnlist, calc_pvals(gv[[1]],gv[[2]]))

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
#' [print.gof()] summaries the diagnostics such as the
#' degree distribution, geodesic distances, shared partner
#' distributions, and reachability for the goodness-of-fit of
#' exponential family random graph models. (`summary.gof` is a deprecated
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
      "b1degree", "bipartition 1 degree", "summary.b1deg",
      "b2degree", "bipartition 2 degree", "summary.b2deg",
      "espartners", "edgewise shared partner", "summary.espart",
      "dspartners", "dyadwise shared partner", "summary.dspart",
      "triadcensus", "triad census", "summary.triadcensus"), 
                      byrow=TRUE, ncol=3)
  all.gof.vars <- sapply(all.gof.vars,
                         match.arg, goftypes[,1])
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
#' [plot.gof()] plots diagnostics such as the degree
#' distribution, geodesic distances, shared partner distributions, and
#' reachability for the goodness-of-fit of exponential family random graph
#' models.
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

  # match variables
  all.gof.vars <- as.character(list_rhs.formula(x$GOF)) %>%
    sapply(match.arg, GOF_VALID_VARS)

  if(length(all.gof.vars) == 0) stop("The gof object does not contain any statistics!")

  GOF <- as.formula(paste("~",paste(all.gof.vars,collapse="+")))

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
      ngs <- nrow(pval) - (idx == "infinite")
      
      pval.max <- 3 + # padding
        if(min(pval[,"MC p-value"]) < 1) max(which(pval[1:ngs, "MC p-value"] < 1))
        else max(which(obs[1:ngs] > 0))

      i <- seq_len(if(is.finite(pval.max) && pval.max <= ngs) pval.max
                   else ngs)
    }
    if(idx == "infinite"){
      i <- c(i, NA, ngs+1)
    }

    if (plotlogodds) {
      odds <- logit(psim)
      odds.obs <- logit(pobs)
      odds.bds <- logit(bds)
      mininf <- min(c(odds, odds.obs, odds.bds) %[f]% is.finite)
      maxinf <- max(c(odds, odds.obs, odds.bds) %[f]% is.finite)

      odds %[.]% (!is.finite(.) & . > 0) <- maxinf
      odds %[.]% (!is.finite(.) & . < 0) <- mininf
      odds.obs %[.]% (!is.finite(.) & . > 0) <- maxinf
      odds.obs %[.]% (!is.finite(.) & . < 0) <- mininf
      odds.bds %[.]% (!is.finite(.) & . > 0) <- maxinf
      odds.bds %[.]% (!is.finite(.) & . < 0) <- mininf

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
    pnames <- colnames(sim)[i]

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

  GVMAP <- list(model = list('model', 'statistic', 'n', 'model statistics', identity),
                degree = list('deg', 'node', 'f', 'degree', identity),
                b1degree = list('b1deg', 'node', 'f', 'b1degree', identity),
                b2degree = list('b2deg', 'node', 'f', 'b2degree', identity),
                odegree = list('odeg', 'node', 'f', 'odegree', identity),
                idegree = list('ideg', 'node', 'f', 'idegree', identity),
                espartners = list('espart', 'edge', 'f', 'edge-wise shared partners', identity),
                dspartners = list('dspart', 'dyad', 'f', 'dyad-wise shared partners', identity),
                triadcensus = list('triadcensus', 'triad', 'n', 'triad census', identity),
                distance = list('dist', 'dyad', 'i', 'minimum geodesic distance', function(gc){
                  ult(gc$pnames) <- "NR"
                  if(normalize.reachability){
                    gc <- within(gc,
                    {
                      mi <- max(gc$i,na.rm=TRUE)
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
                  gc
                }))

  GVMAP <- GVMAP[names(GVMAP)%in%all.gof.vars]
  for(gv in GVMAP)
    gofcomp(gv[[1]], gv[[2]], gv[[3]]) %>% (gv[[5]]) %>% gofplot(gv[[4]])

  mtext(main,side=3,outer=TRUE,cex=1.5,padj=2)
  invisible()
}
