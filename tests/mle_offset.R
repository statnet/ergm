library(statnet.common)
opttest({
library(ergm)
data(florentine)
 
fit1 <- ergm(flomarriage ~ offset(edges) + kstar(2:3), offset.coef=-1, control=control.ergm(MCMLE.maxit=3))

print(summary(fit1))

fit2 <- ergm(flomarriage ~ edges + offset(kstar(2:3)), offset.coef=c(1,-1), control=control.ergm(MCMLE.maxit=3))

print(summary(fit2))
}, "MLE + offset")
