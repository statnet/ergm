library(ergm)
data(sampson)
set.seed(42)
set.edge.attribute(samplike, "YearsTrusted", rbinom(88, 4, .5))
set.seed(296)
set.vertex.attribute(samplike, "YearsServed", rbinom(18, 10, .5))
samplike %v% "Trinity" <- c("F", "S", "H")

num.passed.tests=0
num.tests=0


# transitiveties, directed
 num.tests=num.tests + 1
 s.0 <- summary(samplike~transitiveties)
e.0 <- ergm(samplike~transitiveties, estimate="MPLE")
s.a <- summary(samplike~transitiveties("group"))
e.a <- ergm(samplike~transitiveties("group"), estimate="MPLE")                
s.ad <- summary(samplike~transitiveties("Trinity", diff=TRUE))
e.ad <- ergm(samplike~transitiveties("Trinity", diff=TRUE), estimate="MPLE")   
if (s.0 != 62 || round(e.0$coef + .06997, 3) != 0 ||
    s.a != 18 || round(e.a$coef - .06354, 3) != 0 ||
    !all(s.ad==c(2,0,0)) ||
    !all(round(e.ad$coef + c(.70278, .44099), 3) == 0)) { 
 print(list(s.0=s.0, e.0=e.0, s.a=s.a, e.a=e.a, s.ad=s.ad, e.ad=e.ad))
 stop("Failed transitiveties term test")
} else {
 num.passed.tests=num.passed.tests+1
 print("Passed transitiveties term test")
}

