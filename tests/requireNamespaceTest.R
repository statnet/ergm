library(network)
data(flo)
requireNamespace('ergm')   #load the namespace, but don't attach the package

# run a summary
# gives error Error: The term density does not exist for this type of ERGM. Are you sure you have the right name?
ergm::summary.statistics.formula(as.network(flo)~density)

# try a user-defined ergm term
InitErgmTerm.myedges<-ergm:::InitErgmTerm.edges
ergm::summary.statistics.formula(as.network(flo)~myedges)

# try loading a term from another package namespace
requireNamespace('ergm.userterms')  # NOTE: this now loads ergm normally because it id a Depends of eut
ergm::summary.statistics.formula(network.initialize(2,directed=FALSE)~mindegree(2))
ergm::summary.statistics.formula(network.initialize(2,directed=FALSE)~edges)

requireNamespace('tergm')
par(ask=FALSE)
n<-30
g0<-network.initialize(n,dir=FALSE)

#                     edges, degree(1), mean.age
target.stats<-c(      n*1/2,    n*0.6,        20)

dynfit<-tergm::stergm(g0,formation = ~edges+degree(1), dissolution = ~edges,
               targets = ~edges+degree(1)+mean.age,
               target.stats=target.stats, estimate="EGMME",
               control=tergm::control.stergm(SA.plot.progress=TRUE, SA.restart.on.err=FALSE))

par(ask=TRUE)
ergm::mcmc.diagnostics(dynfit)
summary(dynfit)

# use a summary formula to display number of isolates and edges
                                        # at discrete time points
my.nD <- simulate(dynfit, nsim=1, time.slices=100, output="networkDynamic")
summary(my.nD~isolates+edges+mean.age, at=1:10)
ergm::summary.statistics.formula(my.nD~mean.age)

# try a term that needs to mysteriously access an environment variable
data(sampson, package="ergm")
ergm::summary.statistics.formula(samplike~hammingmix("group"))


# check that we get an appropriate error if no term exists
# should generate an 'unable to locate termed named ... ' error
# ergm::summary.statistics.formula(as.network(flo)~foobar)
