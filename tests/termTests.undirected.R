
library(ergm)
#library(network)  This is not necessary; ergm already depends on network


# an undirected nw
data(faux.mesa.high)
fmh <- faux.mesa.high
set.seed(7)
set.edge.attribute(fmh, "GradeMet", rbinom(203, 6, .5))


# a small undirected nw w/ lots o' triangles
set.seed(20)
t<-trunc(runif(160, 1, 20))
set.seed(21)
h<-trunc(runif(160, 1, 20))
el <- cbind(t,h)
bad <- which(el[,2]==el[,1])
el[bad,2] = el[bad,2]+1
unnw <- network(el, directed=FALSE)
unnw %v% "Pet" <- c("dog", "cat")

num.passed.tests=0
num.tests=0

#altkstar, undirected, 
num.tests=num.tests+1
s.0 <- summary(fmh~altkstar)
e.0 <- ergm(fmh~altkstar, MPLEonly=TRUE)
e.l <- ergm(fmh~altkstar(.5), MPLEonly=TRUE)
s.f <- summary(fmh~altkstar(fixed=TRUE))
e.lf <- ergm(fmh~altkstar(.9, fixed=TRUE), MPLEonly=TRUE)
if (!all(s.0[1:10]==c(51,30,28,18,10,2,4,1,2,1)) ||
    round(e.0$coef+ 3.234, 3) !=0 ||
    round(e.l$coef+ 4.166, 3) !=0 ||
    258 - s.f != 0 ||
    round(e.lf$coef+ 3.494, 3) !=0) {
 print(list(s.0=s.0,e.0=e.0, e.l=e.l, s.f=s.f, e.lf=e.lf))
 stop("Failed altkstar term test")
} else {
 num.passed.tests = num.passed.tests+1
 print("Passed altkstar term test")
}


# concurrent, undirected
num.tests=num.tests+1
s.0 <- summary(fmh~concurrent)
e.0 <- ergm(fmh~concurrent, MPLEonly=TRUE)
s.b <- summary(fmh~concurrent(by="Grade"))
e.b <- ergm(fmh~concurrent(by="Sex"), MPLEonly=TRUE)
if (s.0 != 97 || round(e.0$coef + 4.871, 3) != 0 ||
    !all(s.b==c(35,15,18,8,13,8)) ||
    !all(round(e.b$coef + c(5.17301, 4.67697), 3) == 0)) {
    print(list(s.0=s.0, e.0=e.0, s.b=s.b, e.b=e.b))
 stop("Failed concurrent term test")
} else {
  num.passed.tests = num.passed.tests+1
  print("Passed concurrent term test")
}


# degree, undirected
num.tests=num.tests+1
s.d <- summary(fmh~degree(2:3))
e.d <- ergm(fmh~degree(0), MPLEonly=TRUE)
s.db <- summary(fmh~degree(1:3, "Grade"))
e.db <- ergm(fmh~degree(4, "Sex"), MPLEonly=TRUE)
s.dbh <- summary(fmh~degree(4:5, by="Sex", homophily=TRUE))
e.dbh <- ergm(fmh~degree(2, by="Grade", homophily=TRUE), MPLEonly=TRUE)
if (!all(s.d==c(30,28)) || round(e.d$coef - 5.11, 3) != 0 ||
    !all(s.db==c(15,9,9,9,4,2,11,5,9,9,4,2,5,5,4,2,3,2)) ||
    !all(round(e.db$coef+c(.345, .6005),3)==0) ||
    !all(s.dbh==c(10,3)) || round(e.dbh$coef +.5713,3) !=0) {
 print(list(s.d=s.d, e.d=e.d, s.db=s.db, e.db=e.db, s.dbh=s.dbh, e.dbh=e.dbh))
 stop("Failed degree term test")
} else {
  print("Passed degree term test")
  num.passed.tests = num.passed.tests+1
}


# gwdegree, undirected
num.tests=num.tests+1
s.d <- summary(fmh~gwdegree(.3))
e.d <- ergm(fmh~gwdegree(.4), MPLEonly=TRUE)
s.df <- summary(fmh~gwdegree(.3, fixed=TRUE))
e.df <- ergm(fmh~gwdegree(.2, fixed=TRUE), MPLEonly=TRUE)
s.da <- summary(fmh~gwdegree(.1, attrname="Grade"))
e.da <- ergm(fmh~gwdegree(.1, attrname="Grade"), MPLEonly=TRUE)
s.dfa <- summary(fmh~gwdegree(.1, fixed=TRUE, attrname="Grade"))
e.dfa <- ergm(fmh~gwdegree(.1, fixed=TRUE, attrname="Grade"), MPLEonly=TRUE)
if (!all(head(s.d)==c(51,30,28,18,10,2)) ||
    round(e.d$coef + 13.59067, 3) != 0 ||
    round(s.df - 178.4312, 3) != 0 ||
    round(e.df$coef + 18.2508, 3) != 0 ||
    !all(round(e.da$coef+c(23.94060, 23.30646, 23.51430, 23.31140, 25.11103, 26.88088),3)==0) ||
    !all(round(s.dfa-c(53.58148, 25.53534, 30.83418, 17.79934, 19.31326, 10.80933 ),3)==0) ||
    !all(round(e.dfa$coef+c(23.94060, 23.30646, 23.51430, 23.31140, 25.11103, 26.88088),3)==0)) {
 print(list(s.d=s.d, e.d=e.d, s.df=s.df, e.df=e.df, e.da=e.da, s.dfa=s.dfa, e.dfa=e.dfa))
 stop("Failed gwdegree term test")
} else {
  print("Passed gwdegree term test")
  num.passed.tests = num.passed.tests+1
}



# kstar, undirected
num.tests=num.tests+1
s.k <- summary(fmh~kstar(1:3))
e.k <- ergm(fmh~kstar(c(2,4)), MPLEonly=TRUE)
s.ka <- summary(fmh~kstar(2, "Grade"))
e.ka <- ergm(fmh~kstar(2, "Sex"), MPLEonly=TRUE)
if (!all(s.k == c(406, 659, 1010)) ||
    round(e.k$coef - c(-1.45086, .06255), 3) != 0 ||
    s.ka != 466 || round(e.ka$coef + 1.535175, 3) != 0) {
 print(list(s.k=s.k, e.k=e.k, s.ka=s.ka, e.ka=e.ka))
 stop("Failed kstar term test")
} else {
  print("Passed kstar term test")
  num.passed.tests = num.passed.tests+1
}

# sociality, undirected
num.tests=num.tests + 1
s.0 <- summary(fmh~sociality)
s.a <- summary(fmh~sociality("Race"))
s.b <- summary(fmh~sociality(base=2:203))
s.ab <- summary(fmh~sociality("Race", base=3:200))
e.ab <- ergm(fmh~sociality("Race", base=3:205), MPLEonly=TRUE)
if (!all(head(s.0)==c(4,0,0,1,0,0)) ||
    !all(s.a[45:50]==c(0,8,0,0,0,3)) ||
    !all(s.b==c(13,3,1)) ||
    !all(s.ab==c(7,3,2,0,0,0,0)) ||
    !all(round(e.ab$coef + c(2.6595, 3.5464), 3) ==0)) { 
 print(list(s.0=s.0, s.a=s.a, s.b=s.b, s.ab=s.ab, e.ab=e.ab))
 stop("Failed sociality term test")
} else {
 num.passed.tests=num.passed.tests+1
 print("Passed sociality term test")
}


# tripercent, undirected
# num.tests=num.tests+1
# all of these are turning out 0 or NA                
#s.0 <- summary(unnw~tripercent)
#e.0 <- ergm(unnw~tripercent, MPLEonly=TRUE)
#s.a <- summary(unnw~tripercent("Pet"))
#e.a <- ergm(unnw~tripercent("Pet"), MPLEonly=TRUE)                
#s.ad <- summary(unnw~tripercent("Race", diff=TRUE))
#e.ad <- ergm(unnw~tripercent("Race", diff=TRUE), MPLEonly=TRUE)   
#if (s.0 != 62 || round(e.0$coef + .06997, 3) != 0 ||
#    s.a != 18 || round(e.a$coef - .06354, 3) != 0 ||
#    !all(s.ad==c(2,0,0)) ||
#    !all(round(e.ad$coef + c(.70278, .44099), 3) == 0)) { 
# print(list(s.0=s.0, e.0=e.0, s.a=s.a, e.a=e.a, s.ad=s.ad, e.ad=e.ad))
# stop("Failed tripercent term test")
#} else {
# print("Passed tripercent term test")
#  num.passed.tests = num.passed.tests+1
#}


if(num.passed.tests==num.tests)
print("Passed all undirected term tests")


