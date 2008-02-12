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


