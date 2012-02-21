library(ergm)
data(sampson)

total.theta <- coef(ergm(samplike~edges))
offset.theta <- pi

stopifnot(all.equal(total.theta, coef(ergm(samplike~edges+offset(edges), offset.coef=c(pi)))[1],
          total.theta-offset.theta, tolerance=0.00001))

stopifnot(all.equal(total.theta, coef(ergm(samplike~offset(edges)+edges, offset.coef=c(pi)))[1], 
          total.theta-offset.theta, tolerance=0.00001))
