library(ergm)

logit <- function(p) log(p/(1-p))

nw0 <- network.initialize(5, dir=FALSE)
(nw1 <- simulate(nw0~edges, coef=-.5))
(nw2 <- simulate(nw0~edges, coef=+.5))

(layer <- summary(Layer(nw1, nw2) ~
                    OnLayer(~edges, ~1) +
                    OnLayer(~edges, ~2) +
                    OnLayer(~edges, ~1+2) +
                    OnLayer(~density) +
                    OnLayer(~meandeg) +
                    OnLayer(~edges, ~1&2) +
                    OnLayer(~edges, ~1||2) +
                    OnLayer(~edges, ~(!1)&2) +
                    OnLayer(~edges, ~(1!=2)) +
                    OnLayer(~edges, ~xor(1,2))
                  ))
(logic <- c(summary(nw1~edges),
            summary(nw2~edges),
            summary(nw1~edges)+summary(nw2~edges),
            summary(nw1~density)+summary(nw2~density),
            summary(nw1~meandeg)+summary(nw2~meandeg),
            summary((nw1&nw2)~edges),
            summary((nw1|nw2)~edges),
            summary(((!nw1)&nw2)~edges),
            summary(((nw1&!nw2)|(!nw1&nw2))~edges),
            summary(((nw1&!nw2)|(!nw1&nw2))~edges)
            ))

stopifnot(isTRUE(all.equal(layer, logic, check.attributes=FALSE)))

(layer <- coef(ergm(Layer(nw1, nw2) ~ OnLayer(~edges, ~1) + OnLayer(~edges, ~2))))
(logic <- c(coef(ergm(nw1~edges)), coef(ergm(nw2~edges))))
stopifnot(isTRUE(all.equal(layer, logic, check.attributes=FALSE)))

(layer <- coef(ergm(Layer(nw1, nw2) ~ OnLayer(~edges, ~1+2))))
(logic <- logit((network.edgecount(nw1)+network.edgecount(nw2))/(network.dyadcount(nw1)+network.dyadcount(nw2))))
stopifnot(isTRUE(all.equal(layer, logic, check.attributes=FALSE)))


data(samplk)

(layer <- summary(Layer(samplk1, samplk2)~mutual(layer1=~1,layer2=~2)))
m1 <- as.matrix(samplk1)
m2 <- as.matrix(samplk2)
(logic <- (sum(m1*t(m2)+m2*t(m1))/2))
stopifnot(isTRUE(all.equal(layer, logic, check.attributes=FALSE)))

(layer <- summary(Layer(samplk1, samplk2)~mutual(layer1=~1,layer2=~2&1)))
m1 <- as.matrix(samplk1)
m2 <- as.matrix(samplk2) * as.matrix(samplk1)
(logic <- (sum(m1*t(m2)+m2*t(m1))/2))
stopifnot(isTRUE(all.equal(layer, logic, check.attributes=FALSE)))
