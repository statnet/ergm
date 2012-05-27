#  File ergm/tests/ergmMPLEtest.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
library(ergm)
data(faux.mesa.high)
formula <- faux.mesa.high ~ nodematch("Sex")

mplesetup <- ergmMPLE(formula)
ord <- order(mplesetup$weights)
m <- cbind(mplesetup$weights, mplesetup$response, mplesetup$predictor)[ord,]
if (!all(m - matrix(c(71, 132, 10284, 10423, 1, 1, 0, 0, 0, 1, 1, 0), 4,3) == 0)) {
  stop("Failed first test of ergmMPLE")
}

modelfit <- ergmMPLE(formula, fitmodel=TRUE)
alt <- glm(mplesetup$response ~ mplesetup$predictor - 1, 
           weights = mplesetup$weights, family="binomial")
if(!all(abs(modelfit$coef - alt$coefficients)<1e-4)) {
  stop("Failed second test of ergmMPLE")
}




