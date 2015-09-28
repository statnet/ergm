library(network)
data(flo)
requireNamespace('ergm')   #load the namespace, but don't attach the package

# run a summary
# gives error Error: The term density does not exist for this type of ERGM. Are you sure you have the right name?
ergm::summary.statistics.formula(as.network(flo)~density)



# try loading a term from another package namespace
 
library(ergm.userterms)  # NOTE: this now loads ergm normally because it id a Depends of eut

summary.statistics.formula(network.initialize(2,directed=FALSE)~mindegree(2))
summary.statistics.formula(network.initialize(2,directed=FALSE)~edges)

# try a term that needs to mysteriously access an environment variable
library(ergm)
data(sampson)
s.a <- summary(samplike~hammingmix("group"))
