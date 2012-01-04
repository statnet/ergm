library(ergm)
data(florentine)

efit <- ergm(flomarriage~edges, constraints=~edges)
