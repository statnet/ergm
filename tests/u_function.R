#  File tests/u_function.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################
library(ergm)
n <- 4
nw <- network.initialize(n, dir=TRUE)
nw[1,2] <- 1
nw[1,3] <- 1
nw[3,2] <- 1

# Private storage, auxiliaries, and auxiliaries of auxiliaries

s <- summary(nw~edges+test.abs.edges.minus.5+test.abs.edges.minus.5(FALSE)+sociomatrix+discord.sociomatrix+discord.inter.union.net(nw, implementation="Network"))
stopifnot(all(abs(s[1]-5)==s[2]) && all(abs(s[1]-5)==s[3]))
stopifnot(all(s[3+seq_len(n^2)]==c(as.matrix(nw))))
stopifnot(all(s[3+n^2+seq_len(n^2)]==0))
sim <- simulate(nw~edges,monitor=~test.abs.edges.minus.5+test.abs.edges.minus.5(FALSE)+sociomatrix+discord.sociomatrix+discord.inter.union.net(nw, implementation="Network"), coef=0, nsim=20, control=control.simulate.formula(MCMC.burnin=0,MCMC.interval=1))
s <- attr(sim, "stats")
stopifnot(all(abs(s[,1]-5)==s[,2]) && all(abs(s[,1]-5)==s[,3]))
sim.dyads <- t(sapply(lapply(sim, as.matrix), c))
stopifnot(all(sim.dyads==s[,3+seq_len(n^2)]))
stopifnot(all(t(t(sim.dyads)!=c(as.matrix(nw)))==s[,3+n^2+seq_len(n^2)]))

# Multiple auxiliaries in one term

data(florentine)
floempty <- flomarriage
floempty[,] <- 0

sim <- simulate(flomarriage~edges,monitor=~discord.inter.union.net(floempty, implementation="Network"), coef=0, nsim=100, control=control.simulate.formula(MCMC.burnin=0,MCMC.interval=1), output="stats")
stopifnot(all(sim[,2:4]^2-sim[,5:7]==0))

sim <- simulate(flomarriage~edges,monitor=~discord.inter.union.net(flomarriage, implementation="Network"), coef=0, nsim=100, control=control.simulate.formula(MCMC.burnin=0,MCMC.interval=1), output="stats")
stopifnot(all(sim[,2:4]^2-sim[,5:7]==0))

sim <- simulate(flomarriage~edges,monitor=~discord.inter.union.net(flobusiness, implementation="Network"), coef=0, nsim=100, control=control.simulate.formula(MCMC.burnin=0,MCMC.interval=1), output="stats")
stopifnot(all(sim[,2:4]^2-sim[,5:7]==0))

sim <- simulate(flomarriage~edges,monitor=~discord.inter.union.net(floempty, implementation="DyadSet"), coef=0, nsim=100, control=control.simulate.formula(MCMC.burnin=0,MCMC.interval=1), output="stats")
stopifnot(all(sim[,2:4]^2-sim[,5:7]==0))

sim <- simulate(flomarriage~edges,monitor=~discord.inter.union.net(flomarriage, implementation="DyadSet"), coef=0, nsim=100, control=control.simulate.formula(MCMC.burnin=0,MCMC.interval=1), output="stats")
stopifnot(all(sim[,2:4]^2-sim[,5:7]==0))

sim <- simulate(flomarriage~edges,monitor=~discord.inter.union.net(flobusiness, implementation="DyadSet"), coef=0, nsim=100, control=control.simulate.formula(MCMC.burnin=0,MCMC.interval=1), output="stats")
stopifnot(all(sim[,2:4]^2-sim[,5:7]==0))

# Multiple auxiliaries in one term: multiplicitous proposal

data(florentine)
floempty <- flomarriage
floempty[,] <- 0

sim <- simulate(flomarriage~edges,monitor=~discord.inter.union.net(floempty), coef=0, nsim=100, control=control.simulate.formula(MCMC.burnin=0,MCMC.interval=1), output="stats", constraints=~degrees)
stopifnot(all(sim[,2:4]^2-sim[,5:7]==0))

sim <- simulate(flomarriage~edges,monitor=~discord.inter.union.net(flomarriage), coef=0, nsim=100, control=control.simulate.formula(MCMC.burnin=0,MCMC.interval=1), output="stats", constraints=~degrees)
stopifnot(all(sim[,2:4]^2-sim[,5:7]==0))

sim <- simulate(flomarriage~edges,monitor=~discord.inter.union.net(flobusiness), coef=0, nsim=100, control=control.simulate.formula(MCMC.burnin=0,MCMC.interval=1), output="stats", constraints=~degrees)
stopifnot(all(sim[,2:4]^2-sim[,5:7]==0))

# Private storage: valued

if(require(ergm.count) && packageVersion("ergm.count") >= "4.0"){
nw <- network.initialize(4, dir=TRUE)
nw[1,2,names.eval="v",add.edges=TRUE] <- 1
nw[1,3,names.eval="v",add.edges=TRUE] <- 1
nw[3,2,names.eval="v",add.edges=TRUE] <- 1
nw %ergmlhs% "response" <- "v"
s <- summary(nw~sum+test.abs.sum.minus.5+test.abs.sum.minus.5(FALSE)+test.abs.sum.minus.5(FALSE,TRUE)+sociomatrix)
stopifnot(all(abs(s[1]-5)==s[2]) && all(abs(s[1]-5)==s[3]) && all(abs(s[1]-5)==s[4]))
stopifnot(all(s[-(1:4)]==c(as.matrix(nw, attrname="v"))))
sim <- simulate(nw~sum,monitor=~test.abs.sum.minus.5+test.abs.sum.minus.5(FALSE)+test.abs.sum.minus.5(FALSE,TRUE)+sociomatrix, reference=~Poisson, coef=0, nsim=100, control=control.simulate.formula(MCMC.burnin=0,MCMC.interval=1))
s <- attr(sim, "stats")
stopifnot(all(abs(s[,1]-5)==s[,2]) && all(abs(s[,1]-5)==s[,3]) && all(abs(s[,1]-5)==s[,4]))
sim.dyads <- t(sapply(lapply(sim, as.matrix, attrname="v"), c))
stopifnot(all(sim.dyads==s[,-(1:4)]))

# Multiplicitous proposal
sim <- simulate(nw~sum,monitor=~test.abs.sum.minus.5+test.abs.sum.minus.5(FALSE)+test.abs.sum.minus.5(FALSE,TRUE)+sociomatrix, reference=~DiscUnif(-1,2), coef=0, nsim=100, control=control.simulate.formula(MCMC.burnin=0,MCMC.interval=1,MCMC.prop.weights="random2"))
s <- attr(sim, "stats")
stopifnot(all(abs(s[,1]-5)==s[,2]) && all(abs(s[,1]-5)==s[,3]) && all(abs(s[,1]-5)==s[,4]))
sim.dyads <- t(sapply(lapply(sim, as.matrix, attrname="v"), c))
stopifnot(all(sim.dyads==s[,-(1:4)]))
}
