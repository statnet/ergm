library(ergm)
opttest({
temp.func <- function(nw,deg) {
  fit <- ergm(nw~edges+degree(deg))
  return(fit)
}
mynet <- network.initialize(100,directed=F)
mynet <- simulate(mynet~edges,coef=-2)
aaa <- temp.func(mynet,3)
},"scoping issues (Ticket #193)")
