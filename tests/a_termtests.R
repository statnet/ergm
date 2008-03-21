library(ergm)
data(faux.mesa.high)
data(sampson)

# absdiff
s <- summary(faux.mesa.high ~ absdiff("Grade"))
e <- ergm(faux.mesa.high ~ absdiff("Grade"))
if (s-79 != 0 || round(e$coef + 4.354,3) != 0) {
 stop("Failed absdiff test")
}else{
 print("Passed absdiff test")
}

# absdiffcat
s <- summary(faux.mesa.high ~ absdiffcat("Grade", base=4:5))
e <- ergm(faux.mesa.high ~ absdiffcat("Grade"))
if (!all(s==c(15,15,7)) || 
 !all(round(e$coef+c(6.005,5.788,6.063,6.891,6.611),3)==0)) {
 stop("Failed absdiffcat test")
} else {
 print("Passed absdiffcat test")
}

#altkstar
s <- summary(faux.mesa.high~altkstar(1,fixed=TRUE))
e <- ergm(faux.mesa.high~altkstar(.5, fixed=TRUE), MPLEonly=TRUE)
if (s != 258 || round(e$coef+4.166,3) != 0) {
 stop("Failed altkstar test")
} else {
 print("Passed altkstar test")
}

#asymmetric
s <- summary(samplike~asymmetric)
e <- ergm(samplike~asymmetric, MPLEonly=TRUE)
if (s != 32 || round(e$coef+1.33,3)!=0) {
 stop("Failed asymmetric test")
} else {
 print("Passed asymmetric test")
}

#nodefactor
s <- summary(faux.mesa.high~nodefactor("Grade", base=0))
fmh <- faux.mesa.high
fmh[146,192]<-0 # eliminate only "Other"-race tie to test drop=TRUE
e<-ergm(fmh ~ nodefactor("Race",base=0)+nodefactor("Sex")+nodefactor("Grade",base=3:4))
if (!all(s==c(153, 75, 65, 36, 49, 28)) ||
    !all(round(e$coef[!is.infinite(e$coef)]-c(-1.500, -2.613, -2.203,  -2.163, 
         -0.386, 0.488, 0.181, 0.244, 0.381),3)==0) ||
    e$coef[4]!=-Inf) {
  stop("Failed nodefactor test")
} else {
  print("Passed nodefactor test")
}

#nodeifactor and nodeofactor
m <- samplike
m[,] <- 0
m[1:3,12:16] <- 1 # making a new directed network
s <- summary(m ~ nodeifactor("group", base=0) + nodeofactor("group", base=-(1:3)))
e <- ergm(m ~ nodeifactor("group") + nodeofactor("group"))
if (!all(s==c(9, 0, 6, 0, 15, 0)) ||
    !all(round(e$coef[!is.infinite(e$coef)]-c(-1.719, -0.992),3)==0) ||
    e$coef[1]!=-Inf || e$coef[4]!=-Inf) {
  stop("Failed nodeifactor and nodeofactor test")
} else {
  print("Passed nodeifactor and nodeofactor test")
}


