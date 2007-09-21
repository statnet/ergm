#
pause <- function(){readline(prompt="Pause. Press <Enter> to continue...");invisible()}
#
# Use 'data(package = "ergm")' to list the data sets in it
#
data(package="ergm")
#
# load the Florentine marriage network
#
data(florentine)
pause()
#
# create a plot of the social network
#
plot(flomarriage)
pause()
#
# now make the vertex size proportional to their wealth
#
plot(flomarriage, vertex.cex="wealth", main="Marriage Ties")
pause()
# import synthetic network that looks like a molecule
data(molecule)
# Add a attribute to it to mimic the atomic type
molecule %v% "atomic type" <- c(1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3,3)
pause()
#
# create a plot of the social network
# colored by atomic type
#
plot(molecule, vertex.col="atomic type",vertex.cex=3)
