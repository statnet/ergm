library(ergm)
data(sampson)


degreedist(samplike)

outdegrees <- apply(as.matrix(samplike,m="a"),1,sum)
table(outdegrees)

efit <- ergm(samplike ~ edges + triangle, MPLEonly=T)
summary(efit)

#
# This fit holds the out degrees fixed
#
efit <- ergm(samplike ~ edges + triangle,
 maxit=3,
 constraints=~outdegrees,
 MCMCsamplesize=10000)
summary(efit)

#
# This fit estimates them post estimation of the triangles
#
efit <- ergm(samplike ~ edges + triangle,
 maxit=3, mixed=TRUE, verbose=TRUE,
 constraints=~outdegrees,
 MCMCsamplesize=10000)
summary(efit)
