library(ergm)

data(florentine)

out <- tryCatch(ergm(flomarriage~offset(edges)+edges+degree(1)+offset(degree(0)),target.stats=summary(flomarriage~edges+degree(1)),
              offset.coef=c(0,-0.25),control=control.ergm(init=c(0,-1.47,0.462,-0.25))) ,
         error=function(e) stop('error in target + offset test', e), 
         warning=function(w) stop('unexpected warning in target + offset test', w))

summary(out)
mcmc.diagnostics(out)

set.seed(11)

out <- tryCatch(ergm(flomarriage~offset(edges)+edges+gwdegree(fix=FALSE)+degree(0)+offset(degree(1)),target.stats=summary(flomarriage~edges+gwdegree(fix=FALSE)+degree(0)),
              offset.coef=c(0,-0.25)),
         error=function(e) stop('error in target + offset test', e), 
         warning=function(w) stop('unexpected warning in target + offset test', w))

summary(out)
mcmc.diagnostics(out)
