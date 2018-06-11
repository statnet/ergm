#  File R/ergm-defunct.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2018 Statnet Commons
#######################################################################
# The last home for functions to removed from ergm.

#' @name ergm-defunct
#' @usage sociality(object, ...)
#' @title Functions that have been removed from this package
#' @description Functions that have been removed after a period of deprecation.
#' @param x,minsize,center,cov,inverted,...,object,statistics,formula,init,nsim,burnin,interval,constraints,prop.weights,prop.args,seed,drop,ninflast,V,tol,g,print,nw,radius,probs,n,cols,control,statistic,H Arguments to defunct functions.
#' @keywords internal

NULL

# The following were defunct-ed on 2017-05-16.

#' @rdname ergm-defunct
delete.isolates<-function(x) .Defunct()

#' @rdname ergm-defunct
largest.components<-function(x, minsize=4) .Defunct()

#' @rdname ergm-defunct
central.network<-function(x) .Defunct()

#' @rdname ergm-defunct
ergm.mahalanobis <- function(x, center, cov, inverted=FALSE, ...) .Defunct('stats::mahalanobis')


# sociality() is a special case to avoid an alias conflict with ergm term "sociality".
sociality <- function(object, ...)
  UseMethod("sociality")

#' @rdname ergm-defunct
sociality.default <- function(object,...) .Defunct()

#' @rdname ergm-defunct
sociality.network <- function(object, ..., statistics=NULL) .Defunct("summary.formula and the ergm term 'sociality'")

#' @rdname ergm-defunct
sociality.formula <- function (formula, ..., init, nsim=100,
                               burnin=100, interval=100,
                               constraints=~.,
                               prop.weights="default",
                               prop.args=list(),
                               seed=NULL,  drop=FALSE,
                               statistics=NULL
                               ) .Defunct("summary.formula and the ergm term 'sociality'")

#' @rdname ergm-defunct

sociality.ergm <- function (object, ..., nsim=100,
                            burnin=100, interval=100,
                            constraints=NULL, prop.weights="default", prop.args =list(),
                            seed=NULL, drop=FALSE,
                            statistics=NULL) .Defunct("summary.formula and the ergm term 'sociality'")

#' @rdname ergm-defunct
ostar2deg <- function(object, ninflast=TRUE) .Defunct()

#' @rdname ergm-defunct
is.invertible <- function(V, tol=1e-12) .Defunct("rcond")

#' @rdname ergm-defunct
espartnerdist <- function(g, print=TRUE) .Defunct("summary.formula with 'esp' term")

#' @rdname ergm-defunct
dspartnerdist <- function(g, print=TRUE) .Defunct("summary.formula with 'dsp' term")

#' @rdname ergm-defunct
twopathdist <- function(g, print=TRUE) .Defunct()

#' @rdname ergm-defunct
rspartnerdist <- function (g, print = TRUE) .Defunct("summary.formula with 'esp' and 'dsp' terms")

#' @rdname ergm-defunct
invert.network <- function(nw) .Defunct(msg="!.network")

#' @rdname ergm-defunct
drawpie <- function(center,radius,probs,n=50,cols=1:length(probs),...) .Defunct("latentnet::ergmm.drawpie")

#' @rdname ergm-defunct
mvmodel <- function(object, ...) UseMethod("mvmodel")
mvmodel.default <- function(object,...) stop("Either a network, an ergm object or a formula argument must be given")
mvmodel.formula <- function (formula, ..., init, nsim=100,
                             burnin=10000, interval=1000,
                             constraints=NULL,
                             control=control.simulate.formula(),
                             seed=NULL, 
                             statistic=NULL
		      ) .Defunct(new = 'simulate.formula')

#' @rdname ergm-defunct
mvmodel.ergm <- function (object, ..., nsim=100,
                          burnin=10000, interval=1000,
                          constraints=NULL,
                          seed=NULL,
                          control=control.simulate.ergm(),
                          statistic=NULL) .Defunct(new = 'simulate.ergm')

#' @rdname ergm-defunct
degreedistfactor <- function(g,x) .Defunct("summary.formula with 'degree' terms")

# The following were defunct-ed on 2018-04-07.

#' @rdname ergm-defunct
robust.inverse <- function (H, tol = sqrt(.Machine$double.eps)) .Defunct("MASS::ginv")
