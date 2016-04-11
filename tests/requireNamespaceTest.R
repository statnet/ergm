library(network)
data(flo)
requireNamespace('ergm')   #load the namespace, but don't attach the package

# run a summary
ergm::summary.statistics.formula(as.network(flo)~density)

# try a user-defined ergm term
InitErgmTerm.myedges<-ergm:::InitErgmTerm.edges
ergm::summary.statistics.formula(as.network(flo)~myedges)

# try a term that needs to mysteriously access an environment variable
data(sampson, package="ergm")
ergm::summary.statistics.formula(samplike~hammingmix("group"))


# check that we get an appropriate error if no term exists
# should generate an 'unable to locate termed named ... ' error
# ergm::summary.statistics.formula(as.network(flo)~foobar)
