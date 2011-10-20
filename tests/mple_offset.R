library(ergm)
data(sampson)

total.theta <- coef(ergm(samplike~edges))
offset.theta <- pi

stopifnot(all.equal(total.theta, coef(ergm(samplike~edges+offset(edges), theta0=c(NA,pi)))[1],total.theta-offset.theta, tolerance=0.00001))

stopifnot(all.equal(total.theta, coef(ergm(samplike~offset(edges)+edges, theta0=c(pi,NA)))[1],total.theta-offset.theta, tolerance=0.00001))
