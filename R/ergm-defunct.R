#  File R/ergm-defunct.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################
# The last home for functions to removed from ergm.

#' @name ergm-defunct
#' @usage sociality(object, ...)
#' @title Functions that have been removed from this package
#' @description Functions that have been removed after a period of deprecation.
#' @param x,minsize,center,cov,inverted,...,object,statistics,formula,init,nsim,burnin,interval,constraints,prop.weights,prop.args,seed,drop,ninflast,V,tol,g,print,nw,radius,probs,n,cols,control,statistic,H Arguments to defunct functions.
#' @keywords internal

NULL

# The following were defunct-ed on 2018-04-07.

#' @rdname ergm-defunct
robust.inverse <- function (H, tol = sqrt(.Machine$double.eps)) .Defunct("MASS::ginv")

# The following were defunct-ed on 2019-03-07.
#' @rdname ergm-defunct
plot.network.ergm <- function(x,
    attrname=NULL,
    label=network.vertex.names(x),
    coord=NULL,
    jitter=TRUE,
    thresh=0,
    usearrows=TRUE,
    mode="fruchtermanreingold",
    displayisolates=TRUE,
    interactive=FALSE,
    xlab=NULL,
    ylab=NULL,
    xlim=NULL,
    ylim=NULL,
    pad=0.2,
    label.pad=0.5,
    displaylabels=FALSE,
    boxed.labels=TRUE,
    label.pos=0,
    label.bg="white",
    vertex.sides=8,
    vertex.rot=0,
    arrowhead.cex=1,
    label.cex=1,
    loop.cex=1,
    vertex.cex=1,
    edge.col=1,
    label.col=1,
    vertex.col=2,
    label.border=1,
    vertex.border=1,
    edge.lty=1,
    label.lty=NULL,
    vertex.lty=1,
    edge.lwd=0,
    label.lwd=par("lwd"),
    edge.len=0.5,
    edge.curve=0.1,
    edge.steps=50,
    loop.steps=20,
    object.scale=0.01,
    uselen=FALSE,
    usecurve=FALSE,
    suppress.axes=TRUE,
    vertices.last=TRUE,
    new=TRUE,
    layout.par=NULL,
    cex.main=par("cex.main"),
    cex.sub=par("cex.sub"),
    seed=NULL,
    latent.control=list(maxit=500,trace=0,dyadsample=10000,
               penalty.sigma=c(5,0.5), nsubsample=200),
    colornames="rainbow",
    verbose=FALSE, latent=FALSE, ...)
  .Defunct("latentnet::plot.ergmm()")

#' @rdname ergm-defunct
ergm.getterms<-function(formula) .Defunct("statnet.common::list_rhs.formula() and statnet.common::eval_lhs.formula()")

#' @rdname ergm-defunct
plot.mcmc.list.ergm <- function(...) .Defunct("ergm_plot.mcmc.list()")

#' @rdname ergm-defunct
plot.ergm <- function (x, ...) .Defunct("mcmc.diagnostics(x,...)")

#' @rdname ergm-defunct
summary.statistics <- function(...) .Defunct("summary_formula()")

#' @rdname ergm-defunct
ergm.checkargs <- function(fname, arglist, varnames=NULL, vartypes=NULL,
                           defaultvalues=list(), required=NULL) .Defunct("check.ErgmTerm")

#' @rdname ergm-defunct
ergm.checkbipartite <- function(fname, nw.bipartiteflag, requirement,
                               extramessage="") .Defunct("check.ErgmTerm")

#' @rdname ergm-defunct
ergm.checkdirected <- function(fname, nw.directedflag, requirement,
                               extramessage="") .Defunct("check.ErgmTerm")

#' @rdname ergm-defunct
summary.gof <- function(object, ...) .Defunct("print.gof")

#' @rdname ergm-defunct
ergm.getMCMCsample <- function(nw, model, MHproposal, eta0, control, 
                               verbose=FALSE, response=NULL, update.nws = TRUE,...) .Defunct("ergm_MCMC_sample")

#' @rdname ergm-defunct
ergm.MHP.table <- function(...) .Defunct("ergm_proposal_table()")

#' @rdname ergm-defunct
MHproposal <- function(...) .Defunct("ergm_proposal()")

#' @rdname ergm-defunct
MHproposal.character <- function(...) .Defunct("ergm_proposal()")

#' @rdname ergm-defunct
MHproposal.ergm <- function(...) .Defunct("ergm_proposal()")

#' @rdname ergm-defunct
MHproposal.formula <- function(...) .Defunct("ergm_proposal()")

#' @rdname ergm-defunct
ergm.init.methods <- function(...) .Defunct(msg="Function ergm.init.methods() has been deprecated in favor of specifying init_methods in InitErgmReference.*() functions, and has no effect.")

#' @rdname ergm-defunct
ergm.ConstraintImplications <- function(...) .Defunct(msg="Function ergm.ConstraintImplications() has been deprecated in favor of specifying the implications in the InitErgmConstraint.*() functions, and has no effect.")

#' @rdname ergm-defunct
ergm.mcmcslave <- function(Clist,MHproposal,eta0,control,verbose,...,prev.run=NULL, burnin=NULL, samplesize=NULL, interval=NULL, maxedges=NULL) .Defunct("ergm_MCMC_slave")
