library(ergm)
# Test undirected network
data(faux.mesa.high)
m <- matrix(c(75, 0, 0, 1, 1, 1, 0, 33, 2, 4, 2, 1,
 0, 2, 23, 7, 6, 4, 1, 4, 7, 9, 1, 5, 1, 2, 6, 1, 17, 5,
 1, 1, 4, 5, 5, 6), 6, 6)  # Correct answer!
if (max(abs(m - mixingmatrix(faux.mesa.high, "Grade")$mat)) > 0)
  stop ("mixingmatrix failed test on undirected network faux.mesa.high")
if (max(abs(m[upper.tri(m, diag=T)] - summary(faux.mesa.high ~ nodemix("Grade")))) > 0)
  stop ("nodemix failed test on undirected network faux.mesa.high")


# directed network
data(sampson)
m2 <- matrix(c(30, 9, 7, 5, 23, 1, 1, 2, 10), 3, 3) # Correct answer!
if (max(abs(m2 - mixingmatrix(samplike, "group")$mat)) > 0)
  stop ("mixingmatrix failed test on directed network samplike")
if (max(abs(as.vector(m2) - summary(samplike ~ nodemix("group")))) > 0)
  stop ("nodemix failed test on directed network samplike")

# bipartite network
el <- cbind( c(17, 20, 22, 26, 19, 24, 16, 22, 18, 23, 28, 20,
               22, 23, 17, 21, 25, 21, 27, 16, 19, 18, 23),
           c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 9, 10, 
             10, 11, 11))
mynw <- network(el, bipartite=15, directed=FALSE) 
mynw %v% "names" <- rep(letters[1:3], c(10,10,10))
m3 <- matrix(c(0, 0, 0, 9, 1, 0, 12, 1, 0), 3, 3) # Correct answer!
if (max(abs(m3 - mixingmatrix(mynw, "names")$mat)) > 0)
  stop ("mixingmatrix failed test on bipartite network")
if (max(abs(as.vector(m3) - summary(mynw ~ nodemix("names")))) > 0)
  stop ("nodemix failed test on bipartite network")



