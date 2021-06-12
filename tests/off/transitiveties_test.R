#  File tests/off/transitiveties_test.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################
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
if (s.0 != 62 || round(coef(e.0) + .06997, 3) != 0 ||
    s.a != 18 || round(coef(e.a) - .06354, 3) != 0 ||
    !all(s.ad==c(2,0,0)) ||
    !all(round(coef(e.ad) + c(.70278, .44099), 3) == 0)) { 
 print(list(s.0=s.0, e.0=e.0, s.a=s.a, e.a=e.a, s.ad=s.ad, e.ad=e.ad))
 stop("Failed transitiveties term test")
} else {
 num.passed.tests=num.passed.tests+1
 print("Passed transitiveties term test")
}

