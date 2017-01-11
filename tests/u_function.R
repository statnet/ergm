library(ergm)
nw <- network.initialize(4, dir=TRUE)
nw[1,2] <- 1
nw[1,3] <- 1
nw[3,2] <- 1
s <- summary(nw~edges+test.abs.edges.minus.5+test.abs.edges.minus.5(FALSE))
stopifnot(all(abs(s[1]-5)==s[2]) && all(abs(s[1]-5)==s[3]))
s <- simulate(nw~edges,monitor=~test.abs.edges.minus.5+test.abs.edges.minus.5(FALSE), statsonly=TRUE,coef=0, nsim=1000, control=control.simulate.formula(MCMC.burnin=0,MCMC.interval=1))
stopifnot(all(abs(s[,1]-5)==s[,2]) && all(abs(s[,1]-5)==s[,3]))
