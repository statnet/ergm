#
pause <- function(){readline(prompt="Pause. Press <Enter> to continue...");invisible()}
#
# load the Florentine marriage data matrix
#
data(flo)
#
# attach the sociomatrix for the Florentine marriage data
# This is not yet a network object.
#
flo
pause()
#
# Create a network object out of the adjacency matrix
#
flomarriage <- network(flo,directed=FALSE)
flomarriage
pause()
#
# print out the sociomatrix for the Florentine marriage data
#
flomarriage[,]
pause()
#
# create a vector indicating the wealth of each family (in thousands of lira) 
# and add it as a covariate to the network object
#
flomarriage %v% "wealth" <- c(10,36,27,146,55,44,20,8,42,103,48,49,10,48,32,3)
flomarriage
