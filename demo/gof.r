par(ask=TRUE)
#
data(florentine)
#
# test the gof.ergm function
#
gest <- ergm(flomarriage ~ edges + kstar(2))
gest
summary(gest)

#
# Plot the probabilities first
#
gofflo <- gof(gest)
gofflo
#
# Place all three on the same page
# with nice margins
#
par(mfrow=c(1,3))
par(oma=c(0.5,2,1,0.5))
#
plot(gofflo)
#
# And now the odds 
#
plot(gofflo, plotlogodds=TRUE)
#
# Use the formula version
#
plot(gof(flomarriage ~ edges + kstar(2), theta0=c(-1.6339, 0.0049)))
