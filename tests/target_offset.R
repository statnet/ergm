#  File tests/target_offset.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
library(ergm)

data(florentine)

out <- ergm(flomarriage~offset(edges)+edges+degree(1)+offset(degree(0)),target.stats=summary(flomarriage~edges+degree(1)),
            offset.coef=c(0,-0.25),control=control.ergm(init=c(0,-1.47,0.462,-0.25)))
warnings() -> w

if(length(w)!=1 || names(w)!="Using target.stats for a model with offset terms may produce an inaccurate estimate of the log-likelihood and derived quantities (deviance, AIC, BIC, etc.), because some of the target stats must be imputed.") stop('unexpected warning in target + offset test', w)

summary(out)
mcmc.diagnostics(out)

set.seed(11)

out <- ergm(flomarriage~offset(edges)+edges+gwdegree(fix=FALSE)+degree(0)+offset(degree(1)),target.stats=summary(flomarriage~edges+gwdegree(fix=FALSE)+degree(0)),
            offset.coef=c(0,-0.25))
warnings() -> w

if(length(w)!=1 || names(w)!="Using target.stats for a model with offset terms may produce an inaccurate estimate of the log-likelihood and derived quantities (deviance, AIC, BIC, etc.), because some of the target stats must be imputed.") stop('unexpected warning in target + offset test', w)

summary(out)
mcmc.diagnostics(out)
