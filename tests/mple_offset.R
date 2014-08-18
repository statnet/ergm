library(statnet.common)
opttest({
library(ergm)
data(sampson)

total.theta <- coef(ergm(samplike~edges))
offset.theta <- pi

stopifnot(all.equal(total.theta, coef(ergm(samplike~edges+offset(edges), offset.coef=c(pi)))[1],
          total.theta-offset.theta, tolerance=0.00001))

stopifnot(all.equal(total.theta, coef(ergm(samplike~offset(edges)+edges, offset.coef=c(pi)))[1], 
          total.theta-offset.theta, tolerance=0.00001))

data(florentine)
boo<-flomarriage
boo[1:3,]<-0
foo <- suppressWarnings(
  ergm(flomarriage~edges+offset(edgecov(boo))+gwesp(0.25,fixed=T),offset.coef=20))
if (max(abs(foo$coef), na.rm=T) > 20) stop("MPLE + offset error")

}, "MPLE + offset")
