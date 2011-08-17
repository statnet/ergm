library(ergm)
library(Rglpk)

####Load the data (provided in the package):
data(ecoli)
form <- ecoli2 ~ edges + degree(2:5) + gwdegree(0.25, fixed = T)

m2<-ergm(formula=form, seed=12345, 
        control=control.ergm(style="Stepping", stepMCMCsize=100, gridsize=10000,
        metric="lognormal"), 
        maxit=1, MCMCsamplesize=1000, burnin=1e+4, interval=1000, verb=FALSE)
if (m2$iterations <5 || m2$iterations > 25) stop("Something fishy in stepping test: Iterations = ", m2$iterations)

