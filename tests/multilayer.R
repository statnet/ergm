library(ergm)

logit <- function(p) log(p/(1-p))

nw0 <- network.initialize(10, dir=FALSE)
(nw1 <- simulate(nw0~edges, coef=-.5))
(nw2 <- simulate(nw0~edges, coef=+.5))

(layer <- summary(Layer(A=nw1, B=nw2) ~
                    L(~edges, -1.5~.) +
                    L(~edges, ~A) +
                    L(~edges, ~`2`) +
                    L(~edges, 2~A) +
                    L(~edges, c(~`1`,~B)) +
                    L(~edges, c(1.5~`1`,-.5~B)) +
                    L(~edges, c(~`1`,-.5~B)) +
                    L(~density) +
                    L(~meandeg) +
                    L(~edges, ~A&B) +
                    L(~edges, ~`1`||`2`) +
                    L(~edges, ~(!A)&`2`) +
                    L(~edges, ~(`1`!=`2`)) +
                    L(~edges, ~xor(`1`,B))
                  ))
(logic <- c((summary(nw1~edges)+summary(nw2~edges))*-1.5,
            summary(nw1~edges),
            summary(nw2~edges),
            summary(nw1~edges)*2,
            summary(nw1~edges)+summary(nw2~edges),
            summary(nw1~edges)*1.5+summary(nw2~edges)*-.5,
            summary(nw1~edges)+summary(nw2~edges)*-.5,
            summary(nw1~density)+summary(nw2~density),
            summary(nw1~meandeg)+summary(nw2~meandeg),
            summary((nw1&nw2)~edges),
            summary((nw1|nw2)~edges),
            summary(((!nw1)&nw2)~edges),
            summary(((nw1&!nw2)|(!nw1&nw2))~edges),
            summary(((nw1&!nw2)|(!nw1&nw2))~edges)
            ))

stopifnot(isTRUE(all.equal(layer, logic, check.attributes=FALSE)))

(layer <- coef(ergm(Layer(nw1, nw2) ~ L(~edges, ~`1`) + L(~edges, ~`2`))))
(logic <- c(coef(ergm(nw1~edges)), coef(ergm(nw2~edges))))
stopifnot(isTRUE(all.equal(layer, logic, check.attributes=FALSE)))

(layer <- coef(ergm(Layer(nw1, nw2) ~ L(~edges, c(~`1`,~`2`)))))
(logic <- logit((network.edgecount(nw1)+network.edgecount(nw2))/(network.dyadcount(nw1)+network.dyadcount(nw2))))
stopifnot(isTRUE(all.equal(layer, logic, check.attributes=FALSE)))


data(samplk)

(layer <- summary(Layer(samplk1, samplk2)~mutual(Ls=c(~`1`,~`2`))))
m1 <- as.matrix(samplk1)
m2 <- as.matrix(samplk2)
(logic <- (sum(m1*t(m2)+m2*t(m1))/2))
stopifnot(isTRUE(all.equal(layer, logic, check.attributes=FALSE)))

(layer <- summary(Layer(samplk1, samplk2)~mutual(Ls=c(~`1`,~`2`&`1`))))
m1 <- as.matrix(samplk1)
m2 <- as.matrix(samplk2) * as.matrix(samplk1)
(logic <- (sum(m1*t(m2)+m2*t(m1))/2))
stopifnot(isTRUE(all.equal(layer, logic, check.attributes=FALSE)))

(layer <- summary(Layer(samplk1, samplk2, samplk3)~lCMB))
m1 <- as.matrix(samplk1)
m2 <- as.matrix(samplk2)
m3 <- as.matrix(samplk3)
msum <- m1 + m2 + m3
diag(msum) <- NA
(logic <- sum(lfactorial(msum) + lfactorial(3-msum) - lfactorial(3), na.rm=TRUE))
stopifnot(isTRUE(all.equal(layer, logic, check.attributes=FALSE)))

data(florentine)
(layer <- summary(Layer(m=flomarriage, b=flobusiness)~ddsp(0:10,Ls.path=c(~b,~b))))
(logic <- summary(flobusiness~dsp(0:10)))
stopifnot(isTRUE(all.equal(layer, logic, check.attributes=FALSE)))
(layer <- summary(Layer(m=flomarriage, b=flobusiness)~desp(0:10,Ls.path=c(~b,~b),L.base=~b)))
(logic <- summary(flobusiness~esp(0:10)))
stopifnot(isTRUE(all.equal(layer, logic, check.attributes=FALSE)))
(layer <- summary(Layer(m=flomarriage, b=flobusiness)~dnsp(0:10,Ls.path=c(~b,~b),L.base=~b)))
(logic <- summary(flobusiness~nsp(0:10)))
stopifnot(isTRUE(all.equal(layer, logic, check.attributes=FALSE)))


nw0 <- network.initialize(3, dir=FALSE)
nw1 <- nw0
nw1[1,2] <- 1
nw1[2,3] <- 1

nw2 <- nw0
nw2[1,2] <- 1
nw2[1,3] <- 1

layer_and_MLE <- function(nw1, nw2){
  nd <- network.dyadcount(nw1)
  ne <- summary(nw1&nw2~edges)
  log(3*ne/nd)-log(1-ne/nd)
}

stopifnot(isTRUE(all.equal(layer_and_MLE(nw1,nw2),coef(ergm(Layer(nw1,nw2)~L(~edges, ~`1`&`2`))),check.attributes=FALSE,tolerance=.1)))

# Heterogeneous directedness

nwu <- network.initialize(2, dir=FALSE)
nwu[1,2] <- 1
nwd <- network.initialize(2, dir=TRUE)
nwd[2,1] <- 1

lnw <- Layer(nwu, nwd)
stopifnot(summary(lnw~L(~edges,~`1`)+L(~edges,~`2`)+lCMB)==c(2,1,-sum(lchoose(2,as.matrix(nwu)+as.matrix(nwd)))))
