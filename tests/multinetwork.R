library(ergm)
library(statnet.common)

data(samplk)

logit <- function(p) log(p/(1-p))
samplk1%n%"t" <- 1
samplk2%n%"t" <- 2
samplk3%n%"t" <- 3
samplkl <- list(samplk1, samplk2, samplk3)
samplks <- Networks(samplk1, samplk2, samplk3)

summ.N <- summary(samplks~N(~edges+nodematch("cloisterville"), ~1+t))
summ.l <- unlist(lapply(samplkl, function(nw) summary(nonsimp_update.formula(~edges+nodematch("cloisterville"), nw~., from.new="nw"))))
stopifnot(isTRUE(all.equal(summ.l, summ.N, check.attributes=FALSE)))


ergmfit <- ergm(samplks~N(~edges+nodematch("cloisterville"), ~1+t))

pl <- lapply(samplkl, function(nw) ergmMPLE(nonsimp_update.formula(~edges+nodematch("cloisterville"), nw~., from.new="nw")))
nr <- sapply(lapply(pl, `[[`, "response"),length)

y <- unlist(lapply(pl, `[[`, "response"))
x <- do.call(rbind,lapply(pl, `[[`, "predictor"))
x <- as.data.frame(cbind(x,t=rep(1:3, nr)))
w <- unlist(lapply(pl, `[[`, "weights"))
glmfit <- glm(y~t*nodematch.cloisterville,data=x,weights=w,family="binomial")
stopifnot(all.equal(coef(glmfit),coef(ergmfit),check.attributes=FALSE))
