library(ergm)
data(gflo)
# a markov graph fit to the Florentine data
gest <- ergm(gflo ~ edges + kstar(2), 
	randseed=16124, startatMPLE=TRUE)
gest
summary(gest)
#anova(gest)

#Newton-Raphson iterations:  4
#MCMC sample of size 1000 based on:
#   edges     star2
#-1.66463   0.01181
#
#Monte Carlo MLE Coefficients:
#    edges      star2
#-1.622292   0.006467
