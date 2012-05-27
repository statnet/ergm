#  File ergm/tests/runtime_diags.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
library(ergm)
data(florentine)

gest <- ergm(flomarriage ~ kstar(1:2) + absdiff("wealth") + triangle,
             control=control.ergm(MCMC.burnin=200, MCMC.burnin.retries=10, MCMC.runtime.traceplot=TRUE))

gest <- ergm(flomarriage ~ kstar(1:2) + absdiff("wealth") + triangle,
             control=control.ergm(MCMC.burnin=100, MCMC.burnin.retries=1, MCMC.runtime.traceplot=TRUE))

gest <- ergm(flomarriage ~ kstar(1:2) + absdiff("wealth") + triangle,
             control=control.ergm(MCMC.runtime.traceplot=TRUE))
