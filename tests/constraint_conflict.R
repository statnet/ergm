library(ergm)
data(florentine)

if(!inherits(try(efit <- ergm(flomarriage~edges, constraints=~edges)), 
   "try-error"))
  stop("Should have had an error here.")
