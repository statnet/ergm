library(ergm)

while(exists("test.network"))
  rm(test.network)

lhs.subst.san <- function(rhs,meanstats) {
        test.network <- network.initialize(n=10,directed=F)
        form <- ergm.update.formula(rhs,test.network~.)
        san(form,meanstats=meanstats)
        }



stopifnot(summary(lhs.subst.san(~edges,20)~edges)==20)
