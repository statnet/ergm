library(ergm)
set.seed(0)
data(sampson)

times <- 2
truth <- summary(samplike~edges)*times

stopifnot(summary(samplike~.edges_times, term.options=list(times=2))==truth)

truth <- coef(ergm(samplike~edges))/times
e1 <- ergm(samplike~.edges_times, control=control.ergm(term.options=list(times=2)))
stopifnot(isTRUE(all.equal(coef(e1),truth,check.attributes=FALSE)))


e2 <- ergm(samplike~.edges_times, control=control.ergm(force.main=TRUE,term.options=list(times=2)))
stopifnot(isTRUE(all.equal(coef(e2),truth,check.attributes=FALSE,tolerance=.005)))
stopifnot(isTRUE(all.equal(logLik(e2),logLik(e1),check.attributes=FALSE,tolerance=.005)))

gof(e2)

options(ergm.eval.loglik=FALSE)
e3 <- ergm(samplike~.edges_times, target.stats=as.vector(summary(samplike~edges)/2), control=control.ergm(term.options=list(times=2)))

## stopifnot(isTRUE(all.equal(coef(e2),truth,check.attributes=FALSE,tolerance=.005)))
## stopifnot(isTRUE(all.equal(logLik(e2),logLik(e1),check.attributes=FALSE,tolerance=.005)))

## Now, set globally.
options(ergm.term=list(times=2))
truth <- summary(samplike~edges)*times
stopifnot(summary(samplike~.edges_times)==truth)
## Argument overrides global setting.
options(ergm.term=list(times=1))
stopifnot(summary(samplike~.edges_times, term.options=list(times=2))==truth)
