library(ergm)
data(florentine)

gest <- ergm(flomarriage ~ kstar(1:2) + absdiff("wealth") + triangle,
             control=control.ergm(burnin.retries=10,runtime.traceplot=TRUE),burnin=200)

gest <- ergm(flomarriage ~ kstar(1:2) + absdiff("wealth") + triangle,
             control=control.ergm(burnin.retries=1,runtime.traceplot=TRUE),burnin=100)

gest <- ergm(flomarriage ~ kstar(1:2) + absdiff("wealth") + triangle,
             control=control.ergm(runtime.traceplot=TRUE))
