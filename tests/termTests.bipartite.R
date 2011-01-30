
library(ergm)
#library(network)  This is not necessary; ergm already depends on network

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
bipnw2 <- as.network(exbip.el, matrix.type="edgelist", bipartite=100, directed=FALSE)
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
e.0 <- ergm(bipnw~b1concurrent, MPLEonly=TRUE)
s.b <- summary(bipnw~b1concurrent("Letter"))
e.b <- ergm(bipnw~b1concurrent("Letter"), MPLEonly=TRUE)
if (s.0 != 12 || round(e.0$coef + 3.961, 3) != 0 ||
    !all(s.b==c(4,5,3)) ||
    !all(round(e.b$coef+c(4.143, 3.62, 4.105),3)==0)) {
 print(list(s.0=s.0, e.0=e.0, s.b=s.b, e.b=e.b))
 stop("Failed b1concurrent term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed b1concurrent term test")
}


#b1degree, bipartite, undirected
num.tests=num.tests+1
s.d <- summary(bipnw~b1degree(1:3))
e.d <- ergm(bipnw~b1degree(1:3), MPLEonly=TRUE)
s.db <- summary(bipnw~b1degree(2:4, by="Letter"))
e.db <- ergm(bipnw~b1degree(2, by="Letter"), MPLEonly=TRUE)
if (!all(s.d==c(30,8,2)) ||
    !all(round(e.d$coef+c(2.991, 5.442, 6.484),3)==0) ||
    !all(s.db==c(2,1,1,3,1,1,3,0,0)) ||
    !all(round(e.db$coef+c(1.481, .959, 1.431),3)==0)) {
 print(list(s.d=s.d, e.d=e.d, s.db=s.db, e.db=e.db))
 stop("Failed b1degree term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed b1degree term test")
}




#b1factor, bipartite, undirected
num.tests=num.tests+1
s.a <- summary(bipnw~b1factor("Letter"))
e.a <- ergm(bipnw~b1factor("Letter"), MPLEonly=TRUE)
s.ab <- summary(bipnw~b1factor("Letter", base=3))
e.ab <- ergm(bipnw~b1factor("Letter", base=2), MPLEonly=TRUE)
if (!all(s.a==c(21,19)) ||
    !all(round(e.a$coef+c(3.797, 3.899),3)==0) ||
    !all(s.ab==c(20,21)) ||
    !all(round(e.ab$coef+c(3.877, 3.899),3)==0)) {
 print(list(s.a=s.a, e.a=e.a, s.ab=s.ab, e.ab=e.ab))
 stop("Failed b1factor term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed b1factor term test")
}



#b1star, bipartite, undirected
num.tests=num.tests+1
s.k <- summary(bipnw~b1star(1:2))
e.k <- ergm(bipnw~b1star(1:2), MPLEonly=TRUE)
s.ka <- summary(bipnw~b1star(2:3, "Letter"))
e.ka <- ergm(bipnw~b1star(2:2, "Letter"), MPLEonly=TRUE)
if (!all(s.k==c(60,26)) ||
    !all(round(e.k$coef+c(4.0823, -.3179),3)==0) ||
    !all(s.ka==c(3,0)) ||
    round(e.ka$coef+ 3.157, 3) != 0) {
 print(list(s.k=s.k, e.k=e.k, s.ka=s.ka, e.ka=e.ka))
 stop("Failed b1star term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed b1star term test")
}




#b1starmix, bipartite, undirected
num.tests=num.tests+1
s.ka <- summary(bipnw2~b1starmix(2, "Letter"))
e.ka <- ergm(bipnw2~b1starmix(2, "Letter"), MPLEonly=TRUE)
s.kab <- summary(bipnw2~b1starmix(1, "Letter", base=2))
e.kab <- ergm(bipnw2~b1starmix(1, "Letter", base=2:3), MPLEonly=TRUE)
s.kad <- summary(bipnw2~b1starmix(1, "Letter", diff=FALSE))
e.kad <- ergm(bipnw2~b1starmix(1, "Letter", diff=FALSE), MPLEonly=TRUE)
s.kabd <- summary(bipnw2~b1starmix(1, "Letter", base=2, diff=FALSE))
e.kabd <- ergm(bipnw2~b1starmix(1, "Letter", base=2:3, diff=FALSE), MPLEonly=TRUE)
if (!all(s.ka==c(2,2,5,2)) ||
    !all(round(e.ka$coef+c(5.915, 6.445, 5.084, 6.244),3)==0) ||
    !all(s.kab==c(14,18,18)) ||
    !all(round(e.kab$coef+c(6.275, 6.016),3)==0) ||
    !all(s.kad==c(32,39)) ||
    !all(round(e.kad$coef+c(6.138, 5.940),3)==0) ||
    s.kabd != 32 || round(e.kabd$coef+ 6.138,3)!=0)  {
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
e.a <- ergm(bipnw2~b1twostar("Letter"), MPLEonly=TRUE)
s.aa <- summary(bipnw2~b1twostar("Letter", "Color"))
e.aa <- ergm(bipnw2~b1twostar("Letter", "Color"), MPLEonly=TRUE)
s.ab <- summary(bipnw2~b1twostar("Letter", base=2:4))
e.ab <- ergm(bipnw2~b1twostar("Letter", base=c(1,3,5)), MPLEonly=TRUE)
s.aab <- summary(bipnw2~b1twostar("Letter", "Color", base=2:4))
e.aab <- ergm(bipnw2~b1twostar("Letter", "Color", base=c(1,3,5)), MPLEonly=TRUE)
if (!all(s.a==c(9,4,15,17,7,4)) ||
    !all(round(e.a$coef+c(5.688, 6.012, 5.223, 5.463, 4.635, 5.688),3)==0) ||
    !all(s.aa==c(12,3,18,10,8,5)) ||
    !all(round(e.aa$coef+c(5.174, 6.789, 5.155, 5.299, 4.996, 5.587),3)==0) ||
    !all(s.ab==c(9,7,4)) ||
    !all(round(e.ab$coef+c(6.012, 5.463, 5.688),3)==0) ||
    !all(s.aab==c(12,8,5)) ||
    !all(round(e.aab$coef+c(6.789, 5.299, 5.587),3)==0)) {
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
e.0 <- ergm(bipnw~b2concurrent, MPLEonly=TRUE)
s.b <- summary(bipnw~b2concurrent("Letter"))
e.b <- ergm(bipnw~b2concurrent("Letter"), MPLEonly=TRUE)
if (s.0 != 20 || round(e.0$coef + 3.497, 3) != 0 ||
    !all(s.b==c(8,6,6)) ||
    !all(round(e.b$coef+c(2.803, 4.190, 2.803),3)==0)) {
 print(list(s.0=s.0, e.0=e.0, s.b=s.b, e.b=e.b))
 stop("Failed b2concurrent term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed b2concurrent term test")
}


#b2degree, bipartite, undirected
num.tests=num.tests+1
s.d <- summary(bipnw~b2degree(1:3))
e.d <- ergm(bipnw~b2degree(1:3), MPLEonly=TRUE)
s.db <- summary(bipnw~b2degree(2:4, by="Letter"))
e.db <- ergm(bipnw~b2degree(2, by="Letter"), MPLEonly=TRUE)
if (!all(s.d==c(6,9,8)) ||
    !all(round(e.d$coef-c(1.7203, 1.4941, .6768),3)==0) ||
    !all(s.db==c(3,2,3,3,3,0,3,3,0)) ||
    !all(round(e.db$coef-c(1.0498, -.3001, 1.0217),3)==0)) {
 print(list(s.d=s.d, e.d=e.d, s.db=s.db, e.db=e.db))
 stop("Failed b2degree term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed b2degree term test")
}




#b2factor, bipartite, undirected
num.tests=num.tests+1
s.a <- summary(bipnw~b2factor("Letter"))
e.a <- ergm(bipnw~b2factor("Letter"), MPLEonly=TRUE)
s.ab <- summary(bipnw~b2factor("Letter", base=3))
e.ab <- ergm(bipnw~b2factor("Letter", base=2), MPLEonly=TRUE)
if (!all(s.a==c(19,16)) ||
    !all(round(e.a$coef+c(3.944, 4.119),3)==0) ||
    !all(s.ab==c(25,19)) ||
    !all(round(e.ab$coef+c(3.555, 4.119),3)==0)) {
 print(list(s.a=s.a, e.a=e.a, s.ab=s.ab, e.ab=e.ab))
 stop("Failed b2factor term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed b2factor term test")
}



#b2star, bipartite, undirected
num.tests=num.tests+1
s.k <- summary(bipnw~b2star(1:2))
e.k <- ergm(bipnw~b2star(1:2), MPLEonly=TRUE)
s.ka <- summary(bipnw~b2star(2:3, "Letter"))
e.ka <- ergm(bipnw~b2star(2:2, "Letter"), MPLEonly=TRUE)
if (!all(s.k==c(60,51)) ||
    !all(round(e.k$coef+c(3.3457, .2724),3)==0) ||
    !all(s.ka==c(3,0)) ||
    round(e.ka$coef+ 4.464, 3) != 0) {
 print(list(s.k=s.k, e.k=e.k, s.ka=s.ka, e.ka=e.ka))
 stop("Failed b2star term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed b2star term test")
}




#b2starmix, bipartite, undirected
num.tests=num.tests+1
s.ka <- summary(bipnw2~b2starmix(2, "Letter"))
e.ka <- ergm(bipnw2~b2starmix(1, "Letter"), MPLEonly=TRUE)
s.kab <- summary(bipnw2~b2starmix(1, "Letter", base=2))
e.kab <- ergm(bipnw2~b2starmix(1, "Letter", base=2:3), MPLEonly=TRUE)
s.kad <- summary(bipnw2~b2starmix(1, "Letter", diff=FALSE))
e.kad <- ergm(bipnw2~b2starmix(1, "Letter", diff=FALSE), MPLEonly=TRUE)
s.kabd <- summary(bipnw2~b2starmix(1, "Letter", base=2, diff=FALSE))
e.kabd <- ergm(bipnw2~b2starmix(1, "Letter", base=2:3, diff=FALSE), MPLEonly=TRUE)
if (!all(s.ka==c(1,4,0,2)) ||
    !all(round(e.ka$coef+c(6.275, 5.869, 6.016, 6.016),3)==0) ||
    !all(s.kab==c(14,18,18)) ||
    !all(round(e.kab$coef+c(6.275, 6.016),3)==0) ||
    !all(s.kad==c(71,79)) ||
    !all(round(e.kad$coef+c(6.051, 6.016),3)==0) ||
    s.kabd != 71 || round(e.kabd$coef+ 6.051,3)!=0)  {
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
e.a <- ergm(bipnw2~b2twostar("Letter"), MPLEonly=TRUE)
s.aa <- summary(bipnw2~b2twostar("Letter", "Color"))
e.aa <- ergm(bipnw2~b2twostar("Letter", "Color"), MPLEonly=TRUE)
s.ab <- summary(bipnw2~b2twostar("Letter", base=2:4))
e.ab <- ergm(bipnw2~b2twostar("Letter", base=c(1,3,5)), MPLEonly=TRUE)
s.aab <- summary(bipnw2~b2twostar("Letter", "Color", base=2:4))
e.aab <- ergm(bipnw2~b2twostar("Letter", "Color", base=c(1,3,5)), MPLEonly=TRUE)
if (!all(s.a==c(6,3,16,16,8,6)) ||
    !all(round(e.a$coef+c(5.185, 6.158, 4.892, 4.844, 4.142, 4.901),3)==0) ||
    !all(s.aa==c(10,6,16,11,3,9)) ||
    !all(round(e.aa$coef+c(5.25, 5.698, 4.632, 5.331, 4.931, 4.123),3)==0) ||
    !all(s.ab==c(6,8,6)) ||
    !all(round(e.ab$coef+c(6.158, 4.844, 4.901),3)==0) ||
    !all(s.aab==c(10,3,9)) ||
    !all(round(e.aab$coef+c(5.698, 5.331, 4.123),3)==0)) {
 print(list(s.a=s.a, e.a=e.a, s.aa=s.aa, e.aa=e.aa, s.ab=s.ab, e.ab=e.ab,
            s.aab=s.aab, e.aab=e.aab))
 stop("Failed b2twostar term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed b2twostar term test")
}


# gwb1degree, bipartite
num.tests=num.tests+1
s.d <- summary(bipnw~gwb1degree(.3))
e.d <- ergm(bipnw~gwb1degree(.4), MPLEonly=TRUE)
s.df <- summary(bipnw~gwb1degree(.3, fixed=TRUE))
e.df <- ergm(bipnw~gwb1degree(.2, fixed=TRUE), MPLEonly=TRUE)
s.da <- summary(bipnw~gwb1degree(.1, attrname="Letter"))
e.da <- ergm(bipnw~gwb1degree(.1, attrname="Letter"), MPLEonly=TRUE)
s.dfa <- summary(bipnw~gwb1degree(.1, TRUE, "Letter"))
e.dfa <- ergm(bipnw~gwb1degree(.1, TRUE, "Letter"), MPLEonly=TRUE)
if (round(e.d$coef + 6.979, 3) != 0 ||
    round(s.df -45.4137, 3) != 0 ||
    round(e.df$coef + 8.057, 3) != 0 ||
    !all(round(e.da$coef+c(6.729, 6.762, 6.418),3)==0) ||
    !all(round(s.dfa-c(13.39962, 13.49479, 16.28549),3)==0) ||
    !all(round(e.dfa$coef+c(6.729, 6.762, 6.418),3)==0)) {
 print(list(e.d=e.d, s.df=s.df, e.df=e.df, e.da=e.da, s.dfa=s.dfa, e.dfa=e.dfa))
 stop("Failed gwb1degree term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed gwb1degree term test")
}


# gwb2degree, bipartite
num.tests=num.tests+1
s.d <- summary(bipnw~gwb2degree(.3))
e.d <- ergm(bipnw~gwb2degree(.4), MPLEonly=TRUE)
s.df <- summary(bipnw~gwb2degree(.3, fixed=TRUE))
e.df <- ergm(bipnw~gwb2degree(.2, fixed=TRUE), MPLEonly=TRUE)
s.da <- summary(bipnw~gwb2degree(.1, attrname="Letter"))
e.da <- ergm(bipnw~gwb2degree(.1, attrname="Letter"), MPLEonly=TRUE)
s.dfa <- summary(bipnw~gwb2degree(.1, TRUE, "Letter"))
e.dfa <- ergm(bipnw~gwb2degree(.1, TRUE, "Letter"), MPLEonly=TRUE)
if (round(e.d$coef + 25.99385, 3) != 0 ||
    round(s.df -31.97479, 3) != 0 ||
    round(e.df$coef + 32.78813, 3) != 0 ||
    !all(round(e.da$coef+c(33.82191, 24.76756, 34.28992),3)==0) ||
    !all(round(s.dfa-c(9.809166, 10.598143, 7.598143),3)==0) ||
    !all(round(e.dfa$coef+c(33.82191, 24.76756, 34.28992),3)==0)) {
 print(list(e.d=e.d, s.df=s.df, e.df=e.df, e.da=e.da, s.dfa=s.dfa, e.dfa=e.dfa))
 stop("Failed gwb2degree term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed gwb2degree term test")
}



if(num.passed.tests==num.tests)
  print("Passed all bipartite term tests")
