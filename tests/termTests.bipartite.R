#  File tests/termTests.bipartite.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################

library(ergm)

# a bipartite nw
set.seed(143)
b1 <- floor(runif(60, 1,100))
b2 <- floor(runif(60, 101, 130))
exbip.el <- cbind(b1,b2)
bipnw <- as.network(exbip.el, matrix.type="edgelist", bipartite=100, directed=FALSE)
bipnw %v% "Letter" <- letters[1:3]
bipnw %v% "Cost" <- c(3,2,1)
                          

# another bipartite nw with more ties and 2 attributes
set.seed(258)
b1 <- floor(runif(150, 1,200))
b2 <- floor(runif(150, 201, 400))
exbip.el <- cbind(b1,b2)
bipnw2 <- as.network(exbip.el, matrix.type="edgelist", bipartite=200, directed=FALSE)
bipnw2 %v% "Letter" <- letters[1:2]
color <- rbinom(400, 1, .4)
color[color ==1] <- "Purple"
color[color ==0] <- "Gold"
bipnw2 %v% "Color" <- color


num.passed.tests=0
num.tests=0


#b1concurrent, bipartite, undirected
num.tests=num.tests+1
s.0 <- summary(bipnw~b1concurrent)
e.0 <- ergm(bipnw~b1concurrent, estimate="MPLE")
s.b <- summary(bipnw~b1concurrent("Letter"))
e.b <- ergm(bipnw~b1concurrent(function(x) x %v% "Letter"), estimate="MPLE")
if (s.0 != 12 || round(coef(e.0) + 3.961, 3) != 0 ||
    !all(s.b==c(4,5,3)) ||
    !all(round(coef(e.b)+c(4.143, 3.62, 4.105),3)==0)) {
 print(list(s.0=s.0, e.0=e.0, s.b=s.b, e.b=e.b))
 stop("Failed b1concurrent term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed b1concurrent term test")
}


#b1cov, bipartite, undirected
num.tests=num.tests+1
s.a <- summary(bipnw~b1cov("Cost"))
e.a <- ergm(bipnw~b1cov(~Cost), estimate="MPLE")
s.at <- summary(bipnw~b1cov(~poly(Cost,2,raw=TRUE)))
if (!all(s.a==c(121)) ||
    !all(round(coef(e.a)+c(2.212),3)==0) ||
    !all(s.at==c(121,283))) {
 print(list(s.a=s.a, e.a=e.a, s.at=s.at))
 stop("Failed b1cov term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed b1cov term test")
}


#b1degrange, bipartite, undirected
num.tests=num.tests+1
s.d <- summary(bipnw~b1degrange(from=c(1,2),to=c(Inf,Inf)))
e.d <- ergm(bipnw~b1degrange(from=c(1,2),to=c(Inf,Inf)), estimate="MPLE")
s.anh <- summary(bipnw~b1degrange(from=c(1,2),to=c(Inf,Inf),by="Letter",homophily=FALSE))
e.anh <- ergm(bipnw~b1degrange(from=c(1,2),to=c(Inf,Inf),by=function(x) x %v% "Letter",homophily=FALSE), estimate="MPLE")
s.dh <- summary(bipnw~b1degrange(from=c(1,2),to=c(Inf,Inf),by="Letter",homophily=TRUE))
e.dh <- ergm(bipnw~b1degrange(from=c(1,2),to=c(Inf,Inf),by=function(x) x %v% "Letter",homophily=TRUE), estimate="MPLE")
if (!all(s.d==c(42,12)) ||
		!all(round(coef(e.d)+c(4.027, 3.961),3)==0) ||
        !all(s.anh==c(13,4,13,5,16,3)) ||
        !all(round(coef(e.anh)+c(4.215, 4.143, 4.284, 3.620, 3.636, 4.105),3)==0) ||
		!all(s.dh==c(19,3)) ||
		!all(round(coef(e.dh)+c(3.891, 3.143 ),3)==0)) {
	print(list(s.d=s.d, e.d=e.d, s.dh=s.dh, e.db=e.db))
	stop("Failed b1degree term test")
} else {
	num.passed.tests=num.passed.tests+1
	print("Passed b1degree term test")
}

#b1degree, bipartite, undirected
num.tests=num.tests+1
s.d <- summary(bipnw~b1degree(1:3))
e.d <- ergm(bipnw~b1degree(1:3), estimate="MPLE")
s.db <- summary(bipnw~b1degree(2:4, by="Letter"))
e.db <- ergm(bipnw~b1degree(2, by=~Letter), estimate="MPLE")

if (!all(s.d==c(30,8,2)) ||
    !all(round(coef(e.d)+c(2.991, 5.442, 6.484),3)==0) ||
    !all(s.db==c(2,1,1,3,1,1,3,0,0)) ||
    !all(round(coef(e.db)+c(1.481, .959, 1.431),3)==0)) {
 print(list(s.d=s.d, e.d=e.d, s.db=s.db, e.db=e.db))
 stop("Failed b1degree term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed b1degree term test")
}

#b1dsp, bipartite
num.tests=num.tests+1
s.d0 <- summary(bipnw~b1dsp(0))
e.d0 <- ergm(bipnw~b1dsp(0), estimate="MPLE")
s.d1 <- summary(bipnw~b1dsp(1:3))
e.d1 <- ergm(bipnw~b1dsp(1:3), estimate="MPLE")

if (s.d0 != 4900 ||
    !all(s.d1 == c(49,1,0)) ||
    round(coef(e.d0) - 2.11343,3) !=0 ||    
    !all(round(coef(e.d1)[1:2] + c(2.096629, 3.0399271),3)==0) ||
    !is.infinite(coef(e.d1)[3])) {
 print(list(s.d0=s.d0, e.d0=e.d0, s.d1=s.d1, e.d1=e.d1))
 stop("Failed b1dsp term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed b1dsp term test")
}





#b1factor, bipartite, undirected
num.tests=num.tests+1
s.a <- summary(bipnw~b1factor("Letter"))
e.a <- ergm(bipnw~b1factor(function(x) x %v% "Letter"), estimate="MPLE")
s.ab <- summary(bipnw~b1factor(~Letter, levels=-3))
e.ab <- ergm(bipnw~b1factor("Letter", base=2), estimate="MPLE")
if (!all(s.a==c(21,19)) ||
    !all(round(coef(e.a)+c(3.797, 3.899),3)==0) ||
    !all(s.ab==c(20,21)) ||
    !all(round(coef(e.ab)+c(3.877, 3.899),3)==0)) {
 print(list(s.a=s.a, e.a=e.a, s.ab=s.ab, e.ab=e.ab))
 stop("Failed b1factor term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed b1factor term test")
}




#b1mindegree, bipartite, undirected
num.tests=num.tests+1
s.d <- summary(bipnw~b1mindegree(1:3))
e.d <- ergm(bipnw~b1mindegree(1:3), estimate="MPLE")
if (!all(s.d==c(42,12,4)) ||
		!all(round(coef(e.d)+c(4.027, 3.961, 3.584),3)==0)) {
	print(list(s.d=s.d, e.d=e.d))
	stop("Failed b1mindegree term test")
} else {
	num.passed.tests=num.passed.tests+1
	print("Passed b1mindegree term test")
}

#b1sociality, bipartite, undirected
num.tests=num.tests+1
s.d <- summary(bipnw~b1sociality(nodes=94:96))
e.d <- ergm(bipnw~b1sociality(nodes=94:96), estimate="MPLE")
if (!all(s.d == c(4,4,2)) ||
    !all(round(coef(e.d) + c(1.833, 1.833, 2.603),3) == 0)) {
 print(list(s.d=s.d, e.d=e.d))
 stop("Failed b1sociality term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed b1sociality term test")
}


#b1star, bipartite, undirected
num.tests=num.tests+1
s.k <- summary(bipnw~b1star(1:2))
e.k <- ergm(bipnw~b1star(1:2), estimate="MPLE")
s.ka <- summary(bipnw~b1star(2:3, function(x) x %v% "Letter"))
e.ka <- ergm(bipnw~b1star(2:2, ~Letter), estimate="MPLE")
if (!all(s.k==c(60,26)) ||
    !all(round(coef(e.k)+c(4.0823, -.3179),3)==0) ||
    !all(s.ka==c(3,0)) ||
    round(coef(e.ka)+ 3.157, 3) != 0) {
 print(list(s.k=s.k, e.k=e.k, s.ka=s.ka, e.ka=e.ka))
 stop("Failed b1star term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed b1star term test")
}




#b1starmix, bipartite, undirected
num.tests=num.tests+1
s.ka <- summary(bipnw2~b1starmix(2, "Letter"))
e.ka <- ergm(bipnw2~b1starmix(2, "Letter"), estimate="MPLE")
s.kab <- summary(bipnw2~b1starmix(1, "Letter", base=2))
e.kab <- ergm(bipnw2~b1starmix(1, "Letter", base=2:3), estimate="MPLE")
s.kad <- summary(bipnw2~b1starmix(1, "Letter", diff=FALSE))
e.kad <- ergm(bipnw2~b1starmix(1, "Letter", diff=FALSE), estimate="MPLE")
s.kabd <- summary(bipnw2~b1starmix(1, "Letter", base=2, diff=FALSE))
e.kabd <- ergm(bipnw2~b1starmix(1, "Letter", base=2:3, diff=FALSE), estimate="MPLE")
if (!all(s.ka==c(9,4,7,4)) ||
    !all(round(coef(e.ka)+c(4.870, 5.802, 5.267, 5.962),3)==0) ||
    !all(s.kab==c(36, 39, 40)) ||
    !all(round(coef(e.kab)+c(5.613, 5.497),3)==0) ||
    !all(s.kad==c(75, 75)) ||
    !all(round(coef(e.kad)+c(5.567, 5.567),3)==0) ||
    s.kabd != 75 || round(coef(e.kabd)+ 5.567,3)!=0)  {
 print(list(s.ka=s.ka, e.ka=e.ka, s.kab=s.kab, e.kab=e.kab,
            s.kad=s.kad, e.kad=e.kad, s.kabd=s.kabd, e.kabd=e.kabd))
 stop("Failed b1starmix term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed b1starmix term test")
}



#b1twostar, bipartite, undirected
num.tests=num.tests+1
s.a <- summary(bipnw2~b1twostar("Letter"))
e.a <- ergm(bipnw2~b1twostar(function(x) x %v% "Letter"), estimate="MPLE")
s.aa <- summary(bipnw2~b1twostar("Letter", "Color"))
e.aa <- ergm(bipnw2~b1twostar(~Letter, "Color"), estimate="MPLE")
s.ab <- summary(bipnw2~b1twostar(function(x) x %v% "Letter", levels2=-(2:4)))
e.ab <- ergm(bipnw2~b1twostar("Letter", levels2=-c(1,3,5)), estimate="MPLE")
s.aab <- summary(bipnw2~b1twostar(~Letter, "Color", base=(2:4)))
e.aab <- ergm(bipnw2~b1twostar("Letter", "Color", base=c(1,3,5)), estimate="MPLE")
if (!all(s.a==c(9,4,15,17,7,4)) ||
    !all(round(coef(e.a)+c(4.523, 5.22, 4.773, 4.593, 4.881, 5.446),3)==0) ||
    !all(s.aa==c(13,2,13,15,5,8)) ||
    !all(round(coef(e.aa)+c(4.548, 6.281, 4.882, 4.702, 4.758, 4.364),3)==0) ||
    !all(s.ab==c(9,7,4)) ||
    !all(round(coef(e.ab)+c(5.22, 4.593, 5.446),3)==0) ||
    !all(s.aab==c(13,5,8)) ||
    !all(round(coef(e.aab)+c(6.281, 4.702, 4.364),3)==0)) {
 print(list(s.a=s.a, e.a=e.a, s.aa=s.aa, e.aa=e.aa, s.ab=s.ab, e.ab=e.ab,
            s.aab=s.aab, e.aab=e.aab))
 stop("Failed b1twostar term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed b1twostar term test")
}


  
#b2concurrent, bipartite, undirected
num.tests=num.tests+1
s.0 <- summary(bipnw~b2concurrent)
e.0 <- ergm(bipnw~b2concurrent, estimate="MPLE")
s.b <- summary(bipnw~b2concurrent(function(x) x %v% "Letter"))
e.b <- ergm(bipnw~b2concurrent("Letter"), estimate="MPLE")
if (s.0 != 20 || round(coef(e.0) + 3.497, 3) != 0 ||
    !all(s.b==c(8,6,6)) ||
    !all(round(coef(e.b)+c(2.803, 4.190, 2.803),3)==0)) {
 print(list(s.0=s.0, e.0=e.0, s.b=s.b, e.b=e.b))
 stop("Failed b2concurrent term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed b2concurrent term test")
}


#b2cov, bipartite, undirected
num.tests=num.tests+1
s.a <- summary(bipnw~b2cov("Cost"))
e.a <- ergm(bipnw~b2cov(~Cost), estimate="MPLE")
s.at <- summary(bipnw~b2cov(~poly(Cost,2,raw=TRUE)))
if (!all(s.a==c(129)) ||
    !all(round(coef(e.a)+c(2.191),3)==0) ||
    !all(s.at==c(129,317))) {
 print(list(s.a=s.a, e.a=e.a, s.at=s.at))
 stop("Failed b2cov term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed b2cov term test")
}


#b2degrange, bipartite, undirected
num.tests=num.tests+1
s.d <- summary(bipnw~b2degrange(from=c(1,2),to=c(Inf,Inf)))
e.d <- ergm(bipnw~b2degrange(from=c(1,2),to=c(Inf,Inf)), estimate="MPLE")
s.anh <- summary(bipnw~b2degrange(from=c(1,2),to=c(Inf,Inf),by=function(x) x %v% "Letter",homophily=FALSE))
e.anh <- ergm(bipnw~b2degrange(from=c(1,2),to=c(Inf,Inf),by=~Letter,homophily=FALSE), estimate="MPLE")
s.dh <- summary(bipnw~b2degrange(from=c(1,2),to=c(Inf,Inf),by=function(x) x %v% "Letter",homophily=TRUE))
e.dh <- ergm(bipnw~b2degrange(from=c(1,2),to=c(Inf,Inf),by=~Letter,homophily=TRUE), estimate="MPLE")
if (!all(s.d==c(26,20)) ||
		!all(round(coef(e.d)+c(3.912,3.497),3)==0) ||
        !all(s.anh==c(9,8,10,6,7,6)) ||
        !all(round(coef(e.anh)+c(-13.566, 2.803, -14.365, 4.190, 5.704, 2.803),3)==0) ||
		!all(s.dh==c(19,3)) ||
		!all(round(coef(e.dh)+c(3.03, 4.46 ),3)==0)) {
	print(list(s.d=s.d, e.d=e.d, s.dh=s.dh, e.db=e.db))
	stop("Failed b2degree term test")
} else {
	num.passed.tests=num.passed.tests+1
	print("Passed b2degree term test")
}

#b2degree, bipartite, undirected
num.tests=num.tests+1
s.d <- summary(bipnw~b2degree(1:3))
e.d <- ergm(bipnw~b2degree(1:3), estimate="MPLE")
s.db <- summary(bipnw~b2degree(2:4, by="Letter"))
e.db <- ergm(bipnw~b2degree(2, by=~Letter), estimate="MPLE")
if (!all(s.d==c(6,9,8)) ||
    !all(round(coef(e.d)-c(1.7203, 1.4941, .6768),3)==0) ||
    !all(s.db==c(3,2,3,3,3,0,3,3,0)) ||
    !all(round(coef(e.db)-c(1.0498, -.3001, 1.0217),3)==0)) {
 print(list(s.d=s.d, e.d=e.d, s.db=s.db, e.db=e.db))
 stop("Failed b2degree term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed b2degree term test")
}


#b2dsp, bipartite
num.tests=num.tests+1
s.d0 <- summary(bipnw~b2dsp(0))
e.d0 <- ergm(bipnw~b2dsp(0), estimate="MPLE")
s.d1 <- summary(bipnw~b2dsp(1:3))
e.d1 <- ergm(bipnw~b2dsp(1:3), estimate="MPLE")

if (s.d0 != 381 ||
    !all(s.d1 == c(24,1,0)) ||
    round(coef(e.d0) - 2.829767 ,3) !=0 ||    
    !all(round(coef(e.d1)[1:2] + c(2.804156, 5.140782),3)==0) ||
    !is.infinite(coef(e.d1)[3])) {
 print(list(s.d0=s.d0, e.d0=e.d0, s.d1=s.d1, e.d1=e.d1))
 stop("Failed b2dsp term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed b2dsp term test")
}


#b2mindegree, bipartite, undirected
num.tests=num.tests+1
s.d <- summary(bipnw~b2mindegree(1:3))
e.d <- ergm(bipnw~b2mindegree(1:3), estimate="MPLE")
if (!all(s.d==c(26,20,11)) ||
		!all(round(coef(e.d)+c(3.912,3.497,3.604  ),3)==0)) {
	print(list(s.d=s.d, e.d=e.d))
	stop("Failed b2mindegree term test")
} else {
	num.passed.tests=num.passed.tests+1
	print("Passed b2mindegree term test")
}


#b2factor, bipartite, undirected
num.tests=num.tests+1
s.a <- summary(bipnw~b2factor("Letter"))
e.a <- ergm(bipnw~b2factor(function(x) x %v% "Letter"), estimate="MPLE")
s.ab <- summary(bipnw~b2factor(~Letter, levels=-3))
e.ab <- ergm(bipnw~b2factor("Letter", base=2), estimate="MPLE")
if (!all(s.a==c(19,16)) ||
    !all(round(coef(e.a)+c(3.944, 4.119),3)==0) ||
    !all(s.ab==c(25,19)) ||
    !all(round(coef(e.ab)+c(3.555, 4.119),3)==0)) {
 print(list(s.a=s.a, e.a=e.a, s.ab=s.ab, e.ab=e.ab))
 stop("Failed b2factor term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed b2factor term test")
}


#b2sociality, bipartite, undirected
num.tests=num.tests+1
s.d <- summary(bipnw~b2sociality(nodes=2:6))
e.d <- ergm(bipnw~b2sociality(nodes=2:6), estimate="MPLE")
if (!all(s.d == c(2, 3, 2, 3, 1)) ||
    !all(round(coef(e.d) + c(3.892, 3.476, 3.892, 3.476, 4.595),3) == 0)) {
 print(list(s.d=s.d, e.d=e.d))
 stop("Failed b2sociality term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed b2sociality term test")
}

#b2star, bipartite, undirected
num.tests=num.tests+1
s.k <- summary(bipnw~b2star(1:2))
e.k <- ergm(bipnw~b2star(1:2), estimate="MPLE")
s.ka <- summary(bipnw~b2star(2:3, function(x) x %v% "Letter"))
e.ka <- ergm(bipnw~b2star(2:2, ~Letter), estimate="MPLE")
if (!all(s.k==c(60,51)) ||
    !all(round(coef(e.k)+c(3.3457, .2724),3)==0) ||
    !all(s.ka==c(3,0)) ||
    round(coef(e.ka)+ 4.464, 3) != 0) {
 print(list(s.k=s.k, e.k=e.k, s.ka=s.ka, e.ka=e.ka))
 stop("Failed b2star term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed b2star term test")
}




#b2starmix, bipartite, undirected
num.tests=num.tests+1
s.ka <- summary(bipnw2~b2starmix(2, "Letter"))
e.ka <- ergm(bipnw2~b2starmix(1, "Letter"), estimate="MPLE")
s.kab <- summary(bipnw2~b2starmix(1, "Letter", base=2))
e.kab <- ergm(bipnw2~b2starmix(1, "Letter", base=2:3), estimate="MPLE")
s.kad <- summary(bipnw2~b2starmix(1, "Letter", diff=FALSE))
e.kad <- ergm(bipnw2~b2starmix(1, "Letter", diff=FALSE), estimate="MPLE")
s.kabd <- summary(bipnw2~b2starmix(1, "Letter", base=2, diff=FALSE))
e.kabd <- ergm(bipnw2~b2starmix(1, "Letter", base=2:3, diff=FALSE), estimate="MPLE")
if (!all(s.ka==c(6, 8, 3, 6)) ||
    !all(round(coef(e.ka)+c(5.613, 5.641, 5.523, 5.497),3)==0) ||
    !all(s.kab==c(36, 39, 40)) ||
    !all(round(coef(e.kab)+c(5.613, 5.497),3)==0) ||
    !all(s.kad==c(71,79)) ||
    !all(round(coef(e.kad)+c(5.627, 5.510),3)==0) ||
    s.kabd != 71 || round(coef(e.kabd)+ 5.627,3)!=0)  {
 print(list(s.ka=s.ka, e.ka=e.ka, s.kab=s.kab, e.kab=e.kab,
            s.kad=s.kad, e.kad=e.kad, s.kabd=s.kabd, e.kabd=e.kabd))
 stop("Failed b2starmix term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed b2starmix term test")
}



#b2twostar, bipartite, undirected
num.tests=num.tests+1
s.a <- summary(bipnw2~b2twostar("Letter"))
e.a <- ergm(bipnw2~b2twostar(~Letter), estimate="MPLE")
s.aa <- summary(bipnw2~b2twostar(function(x) x %v% "Letter", "Color"))
e.aa <- ergm(bipnw2~b2twostar("Letter", "Color"), estimate="MPLE")
s.ab <- summary(bipnw2~b2twostar(~Letter, levels2=-(2:4)))
e.ab <- ergm(bipnw2~b2twostar(function(x) x %v% "Letter", base=c(1,3,5)), estimate="MPLE")
s.aab <- summary(bipnw2~b2twostar(~Letter, "Color", base=(2:4)))
e.aab <- ergm(bipnw2~b2twostar("Letter", "Color", levels2=-c(1,3,5)), estimate="MPLE")
if (!all(s.a==c(6,3,16,16,8,6)) ||
    !all(round(coef(e.a)+c(5, 5.754, 4.603, 4.780, 4.462, 5.055),3)==0) ||
    !all(s.aa==c(8, 1, 17, 15, 4, 10)) ||
    !all(round(coef(e.aa)+c(4.823, 6.739, 4.690, 4.702, 5.351, 4.37),3)==0) ||
    !all(s.ab==c(6,8,6)) ||
    !all(round(coef(e.ab)+c(5.754, 4.78, 5.055),3)==0) ||
    !all(s.aab==c(8,4,10)) ||
    !all(round(coef(e.aab)+c(6.739, 4.702, 4.370),3)==0)) {
 print(list(s.a=s.a, e.a=e.a, s.aa=s.aa, e.aa=e.aa, s.ab=s.ab, e.ab=e.ab,
            s.aab=s.aab, e.aab=e.aab))
 stop("Failed b2twostar term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed b2twostar term test")
}





#coincidence, bipartite, undirected, no test on fitting due to the number of coef is large
num.tests=num.tests+1
s.c <- table(summary(bipnw~coincidence,active=0))
#e.c <- table(round(ergm(bipnw~coincidence, estimate="MPLE")$coef,0))
if (!all(s.c==c(381,24,1))) {
	print(list(s.c=s.c, e.c=e.c))
	stop("Failed coincidence term test")
} else {
	num.passed.tests=num.passed.tests+1
	print("Passed coincidence term test")
}

# gwb1degree, bipartite
num.tests=num.tests+1
s.d <- summary(bipnw~gwb1degree())
e.d <- ergm(bipnw~gwb1degree(.4, fixed=TRUE), estimate="MPLE")
s.df <- summary(bipnw~gwb1degree(.3, fixed=TRUE))
e.df <- ergm(bipnw~gwb1degree(.2, fixed=TRUE), estimate="MPLE")
s.da <- summary(bipnw~gwb1degree(.1, fixed=TRUE, attr="Letter"))
e.da <- ergm(bipnw~gwb1degree(.1, attr=function(x) x %v% "Letter", fixed=TRUE), estimate="MPLE")
s.dfa <- summary(bipnw~gwb1degree(.1, TRUE, ~Letter))
e.dfa <- ergm(bipnw~gwb1degree(.1, TRUE, "Letter"), estimate="MPLE")
if (round(coef(e.d) + 6.979, 3) != 0 ||
    round(s.df -45.4137, 3) != 0 ||
    round(coef(e.df) + 8.057, 3) != 0 ||
    !all(round(coef(e.da)+c(6.729, 6.762, 6.418),3)==0) ||
    !all(round(s.dfa-c(13.39962, 13.49479, 16.28549),3)==0) ||
    !all(round(coef(e.dfa)+c(6.729, 6.762, 6.418),3)==0)) {
 print(list(e.d=e.d, s.df=s.df, e.df=e.df, e.da=e.da, s.dfa=s.dfa, e.dfa=e.dfa))
 stop("Failed gwb1degree term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed gwb1degree term test")
}

# gwb1dsp, bipartite
num.tests=num.tests+1
s.d0 <- summary(bipnw~gwb1dsp)
s.d1 <- summary(bipnw~gwb1dsp(.3))
s.d2 <- summary(bipnw~gwb1dsp(.3, TRUE))
e.d2 <- ergm(bipnw~gwb1dsp(.3, TRUE), estimate="MPLE")
s.d3 <- summary(bipnw~gwb1dsp(.3, TRUE, 1))
e.d3 <- ergm(bipnw~gwb1dsp(.3, TRUE, 1), estimate="MPLE")

if (!all(s.d0 == c(49,1,rep(0,27))) ||
    !all(s.d1 == c(49,1,rep(0,27))) ||
    round(s.d2 - 50.25918, 3) != 0 || 
    round(s.d3 - 49, 3) != 0 || 
    round(coef(e.d2) + 2.105815, 3) != 0 ||
    round(coef(e.d3) + 2.072566, 3) != 0) {
 print(list(s.d0=s.d0, s.d1=s.d1, s.d2=s.d2, e.d2=e.d2, s.d3=s.d3, e.d3=e.d3))
 stop("Failed gwb1dsp term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed gwb1dsp term test")
}


# gwb2degree, bipartite
num.tests=num.tests+1
s.d <- summary(bipnw~gwb2degree())
e.d <- ergm(bipnw~gwb2degree(.4, fixed=TRUE), estimate="MPLE")
s.df <- summary(bipnw~gwb2degree(.3, fixed=TRUE))
e.df <- ergm(bipnw~gwb2degree(.2, fixed=TRUE), estimate="MPLE")
s.da <- summary(bipnw~gwb2degree(.1, fixed=TRUE, attr="Letter"))
e.da <- ergm(bipnw~gwb2degree(.1, attr=function(x) x %v% "Letter", fixed=TRUE), estimate="MPLE")
s.dfa <- summary(bipnw~gwb2degree(.1, TRUE, ~Letter))
e.dfa <- ergm(bipnw~gwb2degree(.1, TRUE, "Letter"), estimate="MPLE")
if (!all(summary(bipnw~gwb2degree())==summary(bipnw~b2degree(1:30))) ||
    round(coef(e.d) + 25.99385, 3) != 0 ||
    round(s.df -31.97479, 3) != 0 ||
    round(coef(e.df) + 32.78813, 3) != 0 ||
    !all(round(coef(e.da)+c(33.82191, 24.76756, 34.28992),3)==0) ||
    !all(round(s.dfa-c(9.809166, 10.598143, 7.598143),3)==0) ||
    !all(round(coef(e.dfa)+c(33.82191, 24.76756, 34.28992),3)==0)) {
 print(list(e.d=e.d, s.df=s.df, e.df=e.df, e.da=e.da, s.dfa=s.dfa, e.dfa=e.dfa))
 stop("Failed gwb2degree term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed gwb2degree term test")
}

# gwb2dsp, bipartite
num.tests=num.tests+1
s.d0 <- summary(bipnw~gwb2dsp)
s.d1 <- summary(bipnw~gwb2dsp(.3))
s.d2 <- summary(bipnw~gwb2dsp(.3, TRUE))
e.d2 <- ergm(bipnw~gwb2dsp(.3, TRUE), estimate="MPLE")
s.d3 <- summary(bipnw~gwb2dsp(.3, TRUE, 1))
e.d3 <- ergm(bipnw~gwb2dsp(.3, TRUE, 1), estimate="MPLE")

if (!all(s.d0 == c(24,1,rep(0,28))) ||
    !all(s.d1 == c(24,1,rep(0,28))) ||
    round(s.d2 - 25.25918, 3) != 0 || 
    round(s.d3 - 24, 3) != 0 || 
    round(coef(e.d2) + 2.875923, 3) != 0 ||
    round(coef(e.d3) + 2.220758, 3) != 0) {
 print(list(s.d0=s.d0, s.d1=s.d1, s.d2=s.d2, e.d2=e.d2, s.d3=s.d3, e.d3=e.d3))
 stop("Failed gwb2dsp term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed gwb2dsp term test")
}


if(num.passed.tests==num.tests)
  print("Passed all bipartite term tests")
