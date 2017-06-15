library(ergm)

logit <- function(p) log(p/(1-p))

nw0 <- network.initialize(5, dir=FALSE)
nw1 <- simulate(nw0~edges, coef=-.5)
nw2 <- simulate(nw0~edges, coef=+.5)

stopifnot(isTRUE(all.equal(summary(Layer(nw1, nw2) ~ OnLayer(~edges, ~1) + OnLayer(~edges, ~2) + OnLayer(~edges, ~1+2)), c(summary(nw1~edges),summary(nw2~edges),summary(nw1~edges)+summary(nw2~edges)),check.attributes=FALSE)))

stopifnot(isTRUE(all.equal(coef(ergm(Layer(nw1, nw2) ~ OnLayer(~edges, ~1) + OnLayer(~edges, ~2))), c(coef(ergm(nw1~edges)), coef(ergm(nw2~edges))),check.attributes=FALSE)))

stopifnot(isTRUE(all.equal(coef(ergm(Layer(nw1, nw2) ~ OnLayer(~edges, ~1+2))), logit((network.edgecount(nw1)+network.edgecount(nw2))/(network.dyadcount(nw1)+network.dyadcount(nw2))), check.attributes=FALSE)))
