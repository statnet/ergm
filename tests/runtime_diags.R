library(ergm)
data(florentine)

gest <- ergm(flomarriage ~ kstar(1:2) + absdiff("wealth") + triangle,
             control=control.ergm(MCMC.burnin=200, MCMC.burnin.retries=10, MCMC.runtime.traceplot=TRUE))

gest <- ergm(flomarriage ~ kstar(1:2) + absdiff("wealth") + triangle,
             control=control.ergm(MCMC.burnin=100, MCMC.burnin.retries=1, MCMC.runtime.traceplot=TRUE))

gest <- ergm(flomarriage ~ kstar(1:2) + absdiff("wealth") + triangle,
             control=control.ergm(MCMC.runtime.traceplot=TRUE))
