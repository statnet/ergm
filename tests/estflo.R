#  File tests/estflo.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
#######################################################################
library(statnet.common)
opttest({
library(ergm)
data(florentine)
# a markov graph fit to the Florentine data
gest <- ergm(flomarriage ~ edges + kstar(2), control=control.ergm(seed=16124))
gest
summary(gest)
#anova(gest)

#Newton-Raphson iterations:  4
#MCMC sample of size 1000 based on:
#   edges     star2
#-1.66463   0.01181
#
#Monte Carlo MLE Coefficients:
#    edges      star2
#-1.622292   0.006467

# While we are at it, test the constrainted version.
# (The edges term will be ignored because the constraint makes it irrelevant.)
# XXX uncomment:
#gest <- ergm(flomarriage ~ edges + kstar(2), constraints=~edges, control=control.ergm(seed=16124))
gest <- ergm(flomarriage ~ kstar(2), constraints=~edges, control=control.ergm(seed=16124))
gest
summary(gest)
}, "Florentine")
