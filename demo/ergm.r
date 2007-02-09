pause <- function(){readline(prompt="Pause. Press <Enter> to continue...");invisible()}
#
# Load a network object of the Florentine data
#
data(florentine)
#
# Fit a model where the propensity to form ties between
# families depends on the absolute difference in wealth
#
gest <- ergm(flomarriage ~ edges + absdiff("wealth"))
summary(gest)
#
pause()
#
# add terms for the propensity to form 2-stars and triangles
# of families 
#
gest <- ergm(flomarriage ~ kstar(1:2) + absdiff("wealth") + triangle)
summary(gest)
#
pause()

# import synthetic network that looks like a molecule
#
data(molecule)
# Add a attribute to it to mimic the atomic type
molecule %v% "atomic type" <- c(1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3,3)
#
# measure tendency to match within each atomic type
gest <- ergm(molecule ~ edges + kstar(2) + triangle + nodematch("atomic type"),
  MCMCsamplesize=10000)
summary(gest)
pause()

# compare it to differential homophily by atomic type
gest <- ergm(molecule ~ edges + kstar(2) + triangle + nodematch("atomic type",diff=TRUE),
  MCMCsamplesize=10000)
summary(gest)
