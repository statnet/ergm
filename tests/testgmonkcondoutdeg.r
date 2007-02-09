library(ergm)#, lib.loc="/net/home/handcock/Projects/library")
data(gmonk)

summary(gmonk)

degreedist(gmonk)

outdegrees <- apply(sociomatrix(gmonk),1,sum)
table(outdegrees)

boundDeg <- list(maxout=outdegrees, minout=outdegrees)

efit <- ergm(gmonk ~ edges + triangle, MPLEonly=T)
summary(efit)

#
# This fit holds the out degrees fixed
#
efit <- ergm(gmonk ~ edges + triangle,
 maxit=3,
 proposaltype="condoutdegswap",
 MCMCsamplesize=10000)
summary(efit)

#
# This fit estimates them post estimation of the triangles
#
efit <- ergm(gmonk ~ edges + triangle,
 maxit=3, mixed=TRUE, verbose=TRUE,
 proposaltype="condoutdegswap",
 MCMCsamplesize=10000)
summary(efit)
