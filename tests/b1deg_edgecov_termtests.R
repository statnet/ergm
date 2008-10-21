library(ergm)
data(sampson)

# edgecov
madeup.cov <- matrix(0,18,18)
madeup.cov[1:9,10:18] <- 1
madeup.net <- network(madeup.cov)
s1 <- summary(samplike~edgecov(madeup.cov))
e1 <- ergm(samplike ~ edgecov(madeup.cov))
s2 <- summary(samplike~edgecov(madeup.net))
e2 <- ergm(samplike ~ edgecov(madeup.net))
if (s1!=13 || s2!=13 || round(e1$coef,3)!=-1.655 || round(e2$coef,3)!=-1.655) {
 print(list(s1=s1,e1=e1,s2=s2,e2=e2))
 stop("Failed edgecov test")
}else{
 print("Passed edgecov test")
}

#b1degree
el <- cbind( c(17, 20, 22, 26, 19, 24, 16, 22, 18, 23, 28, 20,
               22, 23, 17, 21, 25, 21, 27, 16, 19, 18, 23),
           c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 9, 10, 
             10, 11, 11))
mynw <- network(el, bipartite=15, directed=FALSE) 
s <- summary(mynw~b1degree(0:5))
#e <- ergm(mynw~edges+b1degree(1:3),MPLEonly=T)
e <- ergm(mynw~b1degree(1:3), MPLEonly=T)
if (!all(s==c(4,2,6,3,0,0)) || 
# !all(round(e$coef,3)==c(-6.270,3.012,8.589,12.867))) {
 !all(round(e$coef,3)==c(-0.693, 0.560, -0.134))) {
 print(list(s=s,e=e))
 stop("Failed b1degree test")
} else {
 print("Passed b1degree test")
}
# NB:  The commented-out lines in the b1degree example above 
# appear to give a horrible MPLE example in which
# the MPLE does not actually exist but there is no indication of this
# fact given by the algorithm!  This is not a bug, except to the extent
# that the glm function in R does not check for this type of behavior.
# I haven't tried to understand this fully, but decreasing the size of
# of the epsilon in glm gives very different results (!)

#b2degree
s <- summary(mynw~b2degree(0:5))
e <- ergm(mynw~edges+b2degree(1:2),MPLEonly=T)
if (!all(s==c(0,5,6,2,0,0)) || 
 !all(round(e$coef,3)==c(-1.367,2.099,1.534))) {
 print(list(s=s,e=e))
 stop("Failed b2degree test")
} else {
 print("Passed b2degree test")
}

