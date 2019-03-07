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
