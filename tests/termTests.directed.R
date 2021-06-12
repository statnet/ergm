#  File tests/termTests.directed.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################

library(ergm)

# a directed nw
load("sampson.wrong.RData") # Old (wrong) version of sampson's monks
set.seed(42)
set.edge.attribute(samplike, "YearsTrusted", rbinom(88, 4, .5))
set.seed(296)
set.vertex.attribute(samplike, "YearsServed", rbinom(18, 10, .5))
samplike %v% "Trinity" <- c("F", "S", "H")

num.passed.tests=0
num.tests=0


#asymmetric, directed
num.tests=num.tests + 1
s.0 <- summary(samplike~asymmetric)
e.0 <- ergm(samplike~asymmetric, estimate="MPLE")
s.a <- summary(samplike~asymmetric("group"))
e.a <- ergm(samplike~asymmetric("group"), estimate="MPLE")
s.ad <- summary(samplike~asymmetric(function(x) x %v% "group", diff=TRUE))
e.ad <- ergm(samplike~asymmetric(~group, diff=TRUE), estimate="MPLE")
s.ak <- summary(samplike~asymmetric("group", levels=3))
e.ak <- ergm(samplike~asymmetric(function(x) x %v% "group", levels=3), estimate="MPLE")
s.adk <- summary(samplike~asymmetric(~group, diff=TRUE, keep=1:2))
e.adk <- ergm(samplike~asymmetric("group", diff=TRUE, keep=c(1,3)), estimate="MPLE")
if (s.0 != 32 || round(coef(e.0)+1.33,3)!=0 ||
    s.a != 17 || round(coef(e.a)+.6008,3)!=0  ||
    !all(s.ad==c(7,2,8)) ||
    !all(round(coef(e.ad)+c(.6931, .6931, .4855),3)==0) ||
    s.ak!=8 || round(coef(e.ak)+.4855,3)!=0 ||
    !all(s.adk==c(7,2)) ||
    !all(round(coef(e.adk)+c(.6931, .4855),3)==0)) {
 print(list(s.0=s.0, e.0=e.0, s.a=s.a, e.a=e.a, s.ad=s.ad, e.ad=e.ad, s.ak=s.ak,
            e.ak=e.ak, s.adk=s.adk, e.adk=e.adk))
 stop("Failed asymmetric term test")
} else {
 print("Passed asymmetric term test")
 num.passed.tests=num.passed.tests+1
}




# ctripe=ctriad, directed
num.tests=num.tests + 1
s.0 <- summary(samplike~ctriple)
e.0 <- ergm(samplike~ctriad, estimate="MPLE")
s.a <- summary(samplike~ctriple(function(x) x %v% "group"))
e.a <- ergm(samplike~ctriple("group"), estimate="MPLE")
s.ad <- summary(samplike~ctriad(~group, diff=TRUE))
e.ad <- ergm(samplike~ctriple(function(x) x %v% "group", diff=TRUE), estimate="MPLE")
if (s.0 != 39 || round(coef(e.0) + .3522, 3) != 0 ||
    s.a != 34 || round(coef(e.a) - .1217, 3) != 0 ||
    !all(s.ad==c(8,4,22)) ||
    !all(round(coef(e.ad)+c(.1949, -.6931, -.2023),3)==0)) {
 print(list(s.0=s.0, e.0=e.0, s.a=s.a, e.a=e.a, s.ad=s.ad, e.ad=e.ad))
 stop("Failed ctriple term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed ctriple term test")
}



# cyclicalties, directed
num.tests=num.tests+1
s.0 <- summary(samplike~cyclicalties)
e.0 <- ergm(samplike~cyclicalties, estimate="MPLE")
s.b <- summary(samplike~cyclicalties(attr=function(x) x %v% "group"))
e.b <- ergm(samplike~cyclicalties(attr="group"), estimate="MPLE")
if (s.0 != 62 || round(coef(e.0) + 0.4154, 3) != 0 ||
		!all(s.b==55) ||
		!all(round(coef(e.b) - 0.2289 , 3) == 0)) {
	print(list(s.0=s.0, e.0=e.0, s.b=s.b, e.b=e.b))
	stop("Failed cyclicalties term test")
} else {
	num.passed.tests = num.passed.tests+1
	print("Passed cyclicalties term test")
}



# idegrange, directed
num.tests=num.tests + 1
s.0 <- summary(samplike~idegrange(5:8))
e.0 <- ergm(samplike~idegrange(5:8), estimate="MPLE")
s.h <- summary(samplike~idegrange(5:8, by="group", homophily=TRUE))
e.h <- ergm(samplike~idegrange(5:8, by=~group, homophily=TRUE), estimate="MPLE")
if (!all(s.0==c(9, 6, 4, 3)) || round(coef(e.0) + c(-0.1431, 1.0986, 1.1451, 0.2231 ))!= 0 ||!all(s.h==c(5, 3, 0, 0)) || round(coef(e.h) + c(-0.5108, -2.1972 , Inf, Inf ))!= 0) {
	print(list(s.0=s.0, e.0=e.0, s.h=s.h, e.h=e.h))
	stop("Failed idegrange term test")
} else {
	print("Passed idegrange term test")
	num.passed.tests = num.passed.tests+1
}





# gwidegree, directed
num.tests=num.tests + 1
s.d <- summary(samplike~gwidegree())
e.d <- ergm(samplike~gwidegree(.4, fixed=TRUE), estimate="MPLE")
s.df <- summary(samplike~gwidegree(.3, fixed=TRUE))
e.df <- ergm(samplike~gwidegree(.2, fixed=TRUE), estimate="MPLE")
s.dfa <- summary(samplike~gwidegree(.1, TRUE, "group"))
e.dfa <- ergm(samplike~gwidegree(.5, TRUE, function(x) x %v% "group"), estimate="MPLE")
if (!all(head(s.d)==c(0,3,5,1,3,2)) ||
    round(coef(e.d) + 5.783202, 3) != 0 ||
    round(s.df - 23.89614, 3) != 0 ||
    round(coef(e.df) + 2.247936, 3) != 0 ||
    !all(round(s.dfa-c(7.715119, 4.408762, 7.734290),3)==0) ||
    !all(round(coef(e.dfa)+c(5.460448, 5.754111, 6.144961),3)==0)) {
 print(list(s.d=head(s.d), e.d=e.d, s.df=s.df, e.df=e.df, e.da=e.da, s.dfa=s.dfa, e.dfa=e.dfa))
 stop("Failed gwidegree term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed gwidegree term test")
}


# gwodegree, directed
num.tests=num.tests + 1
s.d <- summary(samplike~gwodegree())
e.d <- ergm(samplike~gwodegree(.4, fixed=TRUE), estimate="MPLE")
s.df <- summary(samplike~gwodegree(.3, fixed=TRUE))
e.df <- ergm(samplike~gwodegree(.2, fixed=TRUE), estimate="MPLE")
s.dfa <- summary(samplike~gwodegree(.1, TRUE, ~group))
e.dfa <- ergm(samplike~gwodegree(.5, TRUE, "group"), estimate="MPLE")
if (!all(head(s.d)==c(0,0,1,5,7,5)) ||
    round(coef(e.d) + 1.990492, 3) != 0 ||
    round(s.df - 24.23040, 3) != 0 ||
    round(coef(e.df) - 43.61801, 3) != 0 ||
    !all(round(s.dfa-c(7.735906, 4.419631, 7.736070),3)==0) ||
    !all(round(coef(e.dfa)+c(4.1860720, 5.9706455, 0.4921623),3)==0)) {
 print(list(s.d=s.d, e.d=e.d, s.df=s.df, e.df=e.df, e.da=e.da, s.dfa=s.dfa, e.dfa=e.dfa))
 stop("Failed gwodegree term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed gwodegree term test")
}


# idegree, directed
num.tests=num.tests + 1
s.d <- summary(samplike~idegree(2:3))
e.d <- ergm(samplike~idegree(2), estimate="MPLE")
s.db <- summary(samplike~idegree(1:3, "group"))
e.db <- ergm(samplike~idegree(3, function(x) x %v% "group"), estimate="MPLE")
s.dbh <- summary(samplike~idegree(4:5, "group", TRUE))
e.dbh <- ergm(samplike~idegree(2, ~group, TRUE), estimate="MPLE")
if (!all(s.d==c(3,5)) || round(coef(e.d) - 1.223775, 3) != 0 ||
    !all(s.db==c(0,2,1,0,1,2,0,0,2)) ||
    !all(round(coef(e.db)-c(-0.6931472, 0.8183103, 17.4324836),3)==0) ||
    !all(s.dbh==c(3,2)) || round(coef(e.dbh) - .3448,3) !=0) {
 print(list(s.d=s.d, e.d=e.d, s.db=s.db, e.db=e.db, s.dbh=s.dbh, e.dbh=e.dbh))
 stop("Failed idegree term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed idegree term test")
}


# intransitive, directed
num.tests=num.tests + 1
s.0 <- summary(samplike~intransitive)
e.0 <- ergm(samplike~intransitive, estimate="MPLE")
if (s.0 != 224 || round(coef(e.0) + .21, 3) != 0) {
 print(list(s.0=s.0, e.0=e.0))
 stop("Failed intransitive term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed intransitive term test")
}


                     
# idegree1.5, directed
num.tests=num.tests + 1
s.0 <- summary(samplike~idegree1.5)
e.0 <- ergm(samplike~idegree1.5, estimate="MPLE")
if (round(s.0-214.6543,3) != 0 || round(coef(e.0) + .2387, 3) != 0) {
 print(list(s.0=s.0, e.0=e.0))
 stop("Failed idegree1.5 term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed idegree1.5 term test")
}


# istar, directed
num.tests=num.tests + 1
s.k <- summary(samplike~istar(1:3))
e.k <- ergm(samplike~istar(c(2,4)), estimate="MPLE")
s.ka <- summary(samplike~istar(2, function(x) x %v% "group"))
e.ka <- ergm(samplike~istar(2, "group"), estimate="MPLE")
if (!all(s.k == c(88,233,455)) ||
    round(coef(e.k) - c(-.28615, .02477), 3) != 0 ||
    s.ka != 100 || round(coef(e.ka) - .2401, 3) != 0) {
 print(list(s.k=s.k, e.k=e.k, s.ka=s.ka, e.ka=e.ka))
 stop("Failed istar term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed istar term test")
}


# m2star, directed
num.tests=num.tests + 1
s.0 <- summary(samplike~m2star)
e.0 <- ergm(samplike~m2star, estimate="MPLE")
if (s.0 != 378 || round(coef(e.0) + .1028, 3) != 0) {
 print(list(s.0=s.0, e.0=e.0))
 stop("Failed m2star term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed m2star term test")
}


                
# mutual, directed
num.tests=num.tests + 1
s.0 <- summary(samplike~mutual)
e.0 <- ergm(samplike~mutual, estimate="MPLE")
s.s <- summary(samplike~mutual(same=function(x) x %v% "group"))
e.s <- ergm(samplike~mutual(same="group"), estimate="MPLE")
s.b <- summary(samplike~mutual(by=~Trinity))
e.b <- ergm(samplike~mutual(by="Trinity"), estimate="MPLE")
s.sd <- summary(samplike~mutual(same="group", diff=TRUE))
e.sd <- ergm(samplike~mutual(same=function(x) x %v% "group", diff=TRUE), estimate="MPLE")
s.sk <- summary(samplike~mutual(same="group", levels=2))
e.sk <- ergm(samplike~mutual(same=~group, levels=1), estimate="MPLE")
s.bk <- summary(samplike~mutual(by="Trinity", keep=2))
e.bk <- ergm(samplike~mutual(by="Trinity", keep=2:3), estimate="MPLE")
if (s.0 != 28 || round(coef(e.0) - .5596, 3) != 0 ||
    s.s != 23 || round(coef(e.s) - .9954, 3) != 0 ||
    !all(s.b==c(17,18,21)) ||
    !all(round(coef(e.b)-c(.0157, .0667, .8077),3)==0) ||
    !all(s.sd==c(8,4,11)) ||
    !all(round(coef(e.sd)-c(.8267, 1.3863, 1.0116),3)==0) ||
    s.sk != 4 || round(coef(e.sk) - .8266, 3) != 0 ||
    s.bk != 18 || !all(round(coef(e.bk) - c(.0714, .8108)) == 0)){
 print(list(s.0=s.0, e.0=e.0,s.s=s.s, e.s=e.s, s.b=s.b, e.b=e.b,
            s.sd=s.sd, e.sd=e.sd, s.sk=s.sk, e.sk=e.sk,
            s.bk=s.bk, e.bk=e.bk))
 stop("Failed mutual term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed mutual term test")
}


# nearsimmelian, directed
num.tests=num.tests + 1
s.0 <- summary(samplike~nearsimmelian)
e.0 <- ergm(samplike~nearsimmelian, estimate="MPLE")
if (s.0 != 18 || round(coef(e.0) + .4366483, 3) != 0) {
 print(list(s.0=s.0, e.0=e.0))
 stop("Failed nearsimmelian term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed nearsimmelian term test")
}


#nodeicov, directed
num.tests=num.tests + 1
s.a <- summary(samplike~nodeicov("YearsServed"))
e.a <- ergm(samplike~nodeicov(function(x) x %v% "YearsServed"), estimate="MPLE")
s.at <- summary(samplike~nodeicov(~YearsServed^2))
e.at <- ergm(samplike~nodeicov(~(.%v%"YearsServed")^2), estimate="MPLE")
s.att <- summary(samplike~nodeicov(~poly(YearsServed,2,raw=TRUE)))
if (s.a != 439 || round(coef(e.a) + .1739, 3) != 0 ||
    s.at != 2345 || round(coef(e.at) + .02805, 3) != 0 ||
    any(s.att != c(439,2345))) {
 print(list(s.a=s.a, e.a=e.a, s.at=s.at, e.at=e.at, s.att=s.att))
 stop("Failed nodeicov term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed nodeicov term test")
}

                
#nodeifactor, directed
num.tests=num.tests + 1
s.a <- summary(samplike~nodeifactor("group"))
e.a <- ergm(samplike~nodeifactor(~group), estimate="MPLE")
s.ab <- summary(samplike~nodeifactor(function(x) x %v% "Trinity", levels=TRUE))
e.ab <- ergm(samplike~nodeifactor("Trinity", base=(2:3)), estimate="MPLE")
if (!all(s.a==c(13,46)) ||
    !all(round(coef(e.a)+c(1.4424, .4618),3)==0) ||
    !all(s.ab==c(28, 29, 31)) ||
    round(coef(e.ab)+ .9719, 3) != 0) {
  print(list(s.a=s.a,e.a=e.a, s.ab=s.ab, e.ab=e.ab))
  stop("Failed nodeifactor term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed nodeifactor term test")
}



#nodeocov, directed
num.tests=num.tests + 1
s.a <- summary(samplike~nodeocov("YearsServed"))
e.a <- ergm(samplike~nodeocov(function(x) x %v% "YearsServed"), estimate="MPLE")
s.at <- summary(samplike~nodeocov(~YearsServed^2))
e.at <- ergm(samplike~nodeocov(~(.%v%"YearsServed")^2), estimate="MPLE")
s.att <- summary(samplike~nodeocov(~poly(YearsServed,2,raw=TRUE)))
if (s.a != 467 || round(coef(e.a) + .1581, 3) != 0 ||
    s.at != 2691 || round(coef(e.at) + .02243, 3) != 0 ||
    any(s.att != c(467,2691))) {
 print(list(s.a=s.a, e.a=e.a, s.at=s.at, e.at=e.at, s.att=s.att))
 stop("Failed nodeocov term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed nodeocov term test")
}

                
#nodeofactor, directed
num.tests=num.tests + 1
s.a <- summary(samplike~nodeofactor("group"))
e.a <- ergm(samplike~nodeofactor(~group), estimate="MPLE")
s.ab <- summary(samplike~nodeofactor("Trinity", levels=TRUE))
e.ab <- ergm(samplike~nodeofactor(function(x) x %v% "Trinity", base=(2:3)), estimate="MPLE")
if (!all(s.a==c(18,36)) ||
    !all(round(coef(e.a)+c(1.0217, .8353),3)==0) ||
    !all(s.ab==c(31,30,27)) ||
    round(coef(e.ab)+ .8287, 3) != 0) {
  print(list(s.a=s.a,e.a=e.a, s.ab=s.ab, e.ab=e.ab))
  stop("Failed nodeofactor term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed nodeofactor term test")
}


# odegree, directed
num.tests=num.tests + 1
s.d <- summary(samplike~odegree(2:3))
e.d <- ergm(samplike~odegree(3), estimate="MPLE")
s.db <- summary(samplike~odegree(1:3, function(x) x %v% "group"))
e.db <- ergm(samplike~odegree(4, "group"), estimate="MPLE")
s.dbh <- summary(samplike~odegree(4:5, ~group, TRUE))
e.dbh <- ergm(samplike~odegree(2, "group", TRUE), estimate="MPLE")
if (!all(s.d==c(0,1)) || round(coef(e.d) + .1625189, 3) != 0 ||
    !all(s.db==c(0,0,0,0,0,1,0,0,0)) ||
    !all(round(coef(e.db)+c(-1.6292, 0.1112, 0.1625),3)==0) ||
    !all(s.dbh==c(6,3)) || round(coef(e.dbh) + 1.344,3) !=0) {
 print(list(s.d=s.d, e.d=e.d, s.db=s.db, e.db=e.db,
            s.dbh=s.dbh, e.dbh=e.dbh))
 stop("Failed odegree term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed odegree term test")
}



# ostar, directed
num.tests=num.tests + 1
s.k <- summary(samplike~ostar(1:3))
e.k <- ergm(samplike~ostar(c(2,4)), estimate="MPLE")
s.ka <- summary(samplike~ostar(2, "group"))
e.ka <- ergm(samplike~ostar(2, function(x) x %v% "group"), estimate="MPLE")
if (!all(s.k == c(88,178, 191)) ||
    round(coef(e.k) - c(.1224, -.1986), 3) != 0 ||
    s.ka != 88 || round(coef(e.ka) - .1466, 3) != 0) {
 print(list(s.k=s.k, e.k=e.k, s.ka=s.ka, e.ka=e.ka))
 stop("Failed ostar term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed ostar term test")
}



                     
# odegree1.5, directed
num.tests=num.tests + 1
s.0 <- summary(samplike~odegree1.5)
e.0 <- ergm(samplike~odegree1.5, estimate="MPLE")
if (round(s.0-196.9432,3) != 0 || round(coef(e.0) + 0.2909, 3) != 0) {
 print(list(s.0=s.0, e.0=e.0))
 stop("Failed odegree1.5 term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed odegree1.5 term test")
}


# odegrange, directed
num.tests=num.tests + 1
s.0 <- summary(samplike~odegrange(5:8))
e.0 <- ergm(samplike~odegrange(5:8), estimate="MPLE")
s.h <- summary(samplike~odegrange(5:8, by=function(x) x %v% "group", homophily=TRUE))
e.h <- ergm(samplike~odegrange(5:8, by="group", homophily=TRUE), estimate="MPLE")
if (!all(s.0==c(12, 5, 0, 0)) || round(coef(e.0) + c(0.619, -1.030, Inf, Inf))!= 0 ||!all(s.h==c(3, 0, 0, 0)) || round(coef(e.h) + c(-0.2231, Inf , Inf, Inf ))!= 0) {
	print(list(s.0=s.0, e.0=e.0, s.h=s.h, e.h=e.h))
	stop("Failed odegrange term test")
} else {
	print("Passed odegrange term test")
	num.passed.tests = num.passed.tests+1
}


                
# receiver, directed
num.tests=num.tests + 1
s.0 <- summary(samplike~receiver)
e.0 <- ergm(samplike~receiver, estimate="MPLE")
s.b <- summary(samplike~receiver(nodes=-(2:16)))
e.b <- ergm(samplike~receiver(nodes=-(3:18)), estimate="MPLE")
if (!all(s.0==c(8, 4, 2, 5, 3, 5, 7, 11, 10, 6, 3, 6, 3, 5, 3, 2, 3)) ||
    !all(round(coef(e.0)-c(-0.1178,-1.1787,-2.0149,-0.8755,-1.5404,-0.8755,
                          -0.3567, 0.6061, 0.3567,-0.6061,-1.5404,-0.6061,
                          -1.5404,-0.8755,-1.5404,-2.0149,-1.5404),3) == 0) ||
    !all(s.b==c(2,2,3)) ||
    !all(round(coef(e.b)- c(-2.0149, -0.1178),3)==0)) {
 print(list(s.0=s.0, e.0=e.0, s.b=s.b, e.b=e.b))
 stop("Failed receiver term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed receiver term test")
}



# sender, directed
num.tests=num.tests + 1
s.0 <- summary(samplike~sender)
e.0 <- ergm(samplike~sender, estimate="MPLE")
s.b <- summary(samplike~sender(nodes=-(2:16)))
e.b <- ergm(samplike~sender(nodes=-(3:18)), estimate="MPLE")
if (!all(s.0==c(5, 4, 4, 4, 5, 6, 4, 6, 5, 5, 6, 5, 5, 3, 5, 4, 6)) ||
    !all(round(coef(e.0)+c(0.8755,1.1787,1.1787,1.1787,0.8755,0.6061,1.1787,
                          0.6061,0.8755,0.8755,0.6061,0.8755,0.8755,1.5404,
                          0.8755,1.1787,0.6061), 3) == 0) ||
    !all(s.b==c(6,4,6)) ||
    !all(round(coef(e.b)+c(.6061, .8755),3)==0)) {
 print(list(s.0=s.0, e.0=e.0, s.b=s.b, e.b=e.b))
 stop("Failed sender term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed sender term test")
}



                
# simmelian, directed
num.tests=num.tests + 1
s.0 <- summary(samplike~simmelian)
e.0 <- ergm(samplike~simmelian, estimate="MPLE")
if (s.0 != 8 || round(coef(e.0) - .6069, 3) != 0) {
  print(list(s.0=s.0, e.0=e.0))
  stop("Failed simmelian term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed simmelian term test")
}


# simmelianties, directed
num.tests=num.tests + 1
s.0 <- summary(samplike~simmelianties)
e.0 <- ergm(samplike~simmelianties, estimate="MPLE")
if (s.0 != 32 || round(coef(e.0) - .1984, 3) != 0) {
  print(list(s.0=s.0, e.0=e.0))
  stop("Failed simmelianties term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed simmelianties term test")
}



# transitive, directed
num.tests=num.tests + 1
s.0 <- summary(samplike~transitive)
e.0 <- ergm(samplike~transitive, estimate="MPLE")
if (s.0 != 154 || round(coef(e.0) + .07745, 3) != 0) {
 print(list(s.0=s.0, e.0=e.0))
 stop("Failed transitive term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed transitive term test")
}



# transitiveties, directed
num.tests=num.tests + 1
s.0 <- summary(samplike~transitiveties)
e.0 <- ergm(samplike~transitiveties, estimate="MPLE")
if (s.0 != 69 || round(coef(e.0) + 0.4116, 3) != 0) {
	print(list(s.0=s.0, e.0=e.0))
	stop("Failed transitiveties term test")
} else {
	num.passed.tests=num.passed.tests+1
	print("Passed transitiveties term test")
}



# ttriple=ttriad, directed
num.tests=num.tests + 1
s.0 <- summary(samplike~ttriple)
e.0 <- ergm(samplike~ttriad, estimate="MPLE")
s.a <- summary(samplike~ttriple("group"))
e.a <- ergm(samplike~ttriple(~group), estimate="MPLE")
s.ad <- summary(samplike~ttriad("group", diff=TRUE))
e.ad <- ergm(samplike~ttriple(function(x) x %v% "group", diff=TRUE), estimate="MPLE")
if (s.0 != 154 || round(coef(e.0) + .07745, 3) != 0 ||
    s.a != 121 || round(coef(e.a) - .09518, 3) != 0 ||
    !all(s.ad==c(26,14,81)) ||
    !all(round(coef(e.ad)-c(-.05078, .38935, .13469),3)==0)) {
 print(list(s.0=s.0, e.0=e.0, s.a=s.a, e.a=e.a, s.ad=s.ad, e.ad=e.ad))
 stop("Failed ttriple term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed ttriple term test")
}



if (num.passed.tests==num.tests)
  print("Passed all directed term tests")
