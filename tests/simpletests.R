# Simulate a network with a high number of nodes with outdegree=3 and a low number with indegree=3:
library(ergm)
data(sampson)
m <- simulate.formula(samplike~odegree(3)+idegree(3), theta0=c(100,-100))
s <- summary(m~odegree(3)+idegree(3))
if (diff(s) >=0) stop("failed odegree and idegree simulation test.")
 

