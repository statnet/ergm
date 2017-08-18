# The last home for functions to removed from ergm.

# The following were defunct-ed on 2017-05-16.

delete.isolates<-function(x) .Defunct()

largest.components<-function(x, minsize=4) .Defunct()

central.network<-function(x) .Defunct()

ergm.mahalanobis <- function(x, center, cov, inverted=FALSE, ...) .Defunct('stats::mahalanobis')

sociality <- function(object, ...)
UseMethod("sociality")

sociality.default <- function(object,...) .Defunct()

sociality.network <- function(object, ..., statistics=NULL) .Defunct("summary.formula and the ergm term 'sociality'")

sociality.formula <- function (formula, ..., init, nsim=100,
                               burnin=100, interval=100,
                               constraints=~.,
                               prop.weights="default",
                               prop.args=list(),
                               seed=NULL,  drop=FALSE,
                               statistics=NULL
                               ) .Defunct("summary.formula and the ergm term 'sociality'")


sociality.ergm <- function (object, ..., nsim=100,
                            burnin=100, interval=100,
                            constraints=NULL, prop.weights="default", prop.args =list(),
                            seed=NULL, drop=FALSE,
                            statistics=NULL) .Defunct("summary.formula and the ergm term 'sociality'")

ostar2deg <- function(object, ninflast=TRUE) .Defunct()

is.invertible <- function(V, tol=1e-12) .Defunct("rcond")

espartnerdist <- function(g, print=TRUE) .Defunct("summary.formula with 'esp' term")

dspartnerdist <- function(g, print=TRUE) .Defunct("summary.formula with 'dsp' term")

twopathdist <- function(g, print=TRUE) .Defunct()

rspartnerdist <- function (g, print = TRUE) .Defunct("summary.formula with 'esp' and 'dsp' terms")

invert.network <- function(nw) .Defunct(msg="!.network")

drawpie <- function(center,radius,probs,n=50,cols=1:length(probs),...) .Defunct("latentnet::ergmm.drawpie")

mvmodel <- function(object, ...) UseMethod("mvmodel")
mvmodel.default <- function(object,...) stop("Either a network, an ergm object or a formula argument must be given")
mvmodel.formula <- function (formula, ..., init, nsim=100,
                             burnin=10000, interval=1000,
                             constraints=NULL,
                             control=control.simulate.formula(),
                             seed=NULL, 
                             statistic=NULL
		      ) .Defunct(new = 'simulate.formula')

mvmodel.ergm <- function (object, ..., nsim=100,
                          burnin=10000, interval=1000,
                          constraints=NULL,
                          seed=NULL,
                          control=control.simulate.ergm(),
                          statistic=NULL) .Defunct(new = 'simulate.ergm')

degreedistfactor <- function(g,x) .Defunct("summary.formula with 'degree' terms")

