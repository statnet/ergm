#  File tests/termTests.flexible.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
#######################################################################

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
bipnw2 <- as.network(exbip.el, matrix.type="edgelist", bipartite=100, directed=FALSE)
bipnw2 %v% "Letter" <- letters[1:2]
color <- rbinom(400, 1, .4)
color[color ==1] <- "Purple"
color[color ==0] <- "Gold"
bipnw2 %v% "Color" <- color


# a directed nw
data(sampson)
set.seed(42)
set.edge.attribute(samplike, "YearsTrusted", rbinom(88, 4, .5))
set.seed(296)
set.vertex.attribute(samplike, "YearsServed", rbinom(18, 10, .5))
samplike %v% "Trinity" <- c("F", "S", "H")


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

# absdiff, no required type, independent
num.tests=num.tests+1
s.a <- summary(fmh ~ absdiff("Grade"))
e.a <- ergm(fmh ~ absdiff("Grade"))
s.ap <- summary(fmh ~ absdiff("Grade", pow=2))
e.ap <- ergm(fmh ~ absdiff("Grade", pow=2))

if (s.a-79 != 0 || round(e.a$coef + 4.354,3) != 0 ||
    s.ap-195 != 0 || round(e.ap$coef + 3.41,3) != 0) {
 print(list(s.a=s.a,e.a=e.a, s.ap=s.ap, e.ap=e.ap))
 stop("Failed absdiff term test")
}else{
 num.passed.tests=num.passed.tests+1
 print("Passed absdiff term test")
}


# absdiffcat, no required type, independent
num.tests=num.tests+1
s.a <- summary(fmh ~ absdiffcat("Grade"))
e.a <- ergm(fmh ~ absdiffcat("Grade"))
s.ab <- summary(fmh ~ absdiffcat("Grade", base=4:5))
e.ab <- ergm(fmh ~ absdiffcat("Grade", base=4:5))
if (!all(s.a==c(15,15,7,2,1)) ||
    !all(round(e.a$coef+c(6.005,5.788,6.063,6.891,6.611),3)==0) ||
    !all(s.ab==c(15,15,7)) ||
    !all(round(e.ab$coef+c(6.005,5.788,6.063),3)==0)) {
 print(list(s.a=s.a,e.a=e.a, s.ab=s.ab, e.ab=e.ab))
 stop("Failed absdiffcat term test")
} else {
 num.passed.tests=num.passed.tests+1
 print("Passed absdiffcat term test")
}



# balance, dir or undir
num.tests=num.tests+1
s.0 <- summary(fmh~balance)
e.0 <- ergm(fmh~balance, estimate="MPLE")
if (s.0 != 40139 || round(e.0$coef + .02376, 3) != 0) {
 print(list(s.0=s.0, e.0=e.0))
 stop("Failed balance term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed balance term test")
}



# cycle, either
num.tests=num.tests+1
s.k <- summary(fmh~cycle(3:6))
e.k <- ergm(fmh~cycle(c(4,6)), estimate="MPLE")
if(!all(s.k==c(62,80,138,270)) ||
    !all(round(e.k$coef+c(-.1615, .2083),3)==0)) {
 print(list(s.k=s.k,e.k=e.k))
 stop("Failed cycle test")
} else {
 num.passed.tests=num.passed.tests+1
 print("Passed cycle test")
}



# density, either
num.tests=num.tests+1
s.0 <- summary(fmh~density)
e.0 <- ergm(samplike~density, estimate="MPLE")
if (round(s.0 - .009708274,3) != 0 ||
    round(e.0$coef + 277.5904, 3) != 0) {
 print(list(s.0=s.0, e.0=e.0))
 stop("Failed density term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed density term test")
}



# dsp, either
num.tests=num.tests+1
s.d <- summary(fmh~dsp(2:3))
e.d <- ergm(samplike~dsp(4), estimate="MPLE")
if (!all(s.d==c(75,23)) ||
    round(e.d$coef + .04275, 3) != 0) {
   print(list(s.d=s.d, e.d=e.d))
    stop("Failed dsp term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed dsp term test")
}



# dyadcov, either
num.tests=num.tests+1
set.seed(120)
cov <- matrix(rbinom(324, 1, .5),18,18)
cov <- cov+t(cov)
s.x <- summary(samplike~dyadcov(cov))
e.x <- ergm(samplike ~ dyadcov(cov))
s.xa <- summary(fmh~dyadcov(fmh, "GradeMet"))
e.xa <- ergm(fmh ~ dyadcov(fmh, "GradeMet"))
if (!all(s.x==c(31,21,14)) ||
    !all(round(e.x$coef+c(.8546, 1.0732, 1.3467),3)==0) ||
    s.xa!=641 || round(e.xa$coef - 12.31787,3)!=0) {
 print(list(s.x=s.x,e.x=e.x,s.xa=s.xa,e.xa=e.xa))
 stop("Failed dyadcov test")
}else{
 num.passed.tests=num.passed.tests+1
 print("Passed dyadcov test")
}

  

# edgecov, either
num.tests=num.tests+1
set.seed(64)
cov <- matrix(rbinom(324, 3, .5),18,18)
s.x <- summary(samplike~edgecov(cov))
e.x <- ergm(samplike ~ edgecov(cov))
s.xa <- summary(samplike~edgecov(samplike, "YearsTrusted"))
e.xa <- ergm(samplike ~ edgecov(samplike, "YearsTrusted"))
n.x <- try(summary(samplike~edgecov('dummy')),silent=TRUE)
set.network.attribute(samplike,'dummy',cov)
n2.x <- summary(samplike~edgecov('dummy'))
if (s.x!=134 || round(e.x$coef + .5022,3)!=0  ||
    s.xa!=183 || e.xa$coef!=+Inf ||
    !is(n.x,'try-error') || n2.x!=134) {
 print(list(s.x=s.x,e.x=e.x,s.xa=s.xa,e.xa=e.xa))
 stop("Failed edgecov test")
}else{
 num.passed.tests=num.passed.tests+1
 print("Passed edgecov test")
}



# edges, either
num.tests=num.tests+1
s.0 <- summary(fmh~edges)
e.0 <- ergm(samplike~edges, estimate="MPLE")
if (s.0 != 203 || round(e.0$coef + .9072, 3) != 0) {
 print(list(s.0=s.0, e.0=e.0))
 stop("Failed edges term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed edges term test")
}



# esp, either
num.tests=num.tests+1
s.d <- summary(fmh~esp(2:3))
e.d <- ergm(samplike~esp(4), estimate="MPLE")
if (!all(s.d==c(36,13)) ||
    round(e.d$coef - .3093, 3) != 0) {
   print(list(s.d=s.d, e.d=e.d))
    stop("Failed esp term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed esp term test")
}




# gwdsp, either
num.tests=num.tests+1
s.0 <- summary(fmh~gwdsp)
e.0 <- ergm(samplike~gwdsp(fixed=TRUE), estimate="MPLE")
e.a <- ergm(samplike~gwdsp(.8, fixed=TRUE), estimate="MPLE")
s.f <- summary(fmh~gwdsp(fixed=TRUE))
s.af <- summary(fmh~gwdsp(.3, fixed=TRUE))
e.af <- ergm(samplike~gwdsp(.2, fixed=TRUE), estimate="MPLE")
if (!all(head(s.0)==c(431, 75, 23, 1, 1, 0)) ||
    round(e.0$coef + .3309974, 3) != 0 ||
    round(e.a$coef + .1875983, 3) != 0 ||
    s.f!=531 ||
    round(s.af - 558.6369, 3) != 0 ||
    round(e.af$coef + .2829672, 3) != 0) {
  print(list(s.0=head(s.0), e.0=e.0, e.a=e.a, s.f=s.f, s.af=s.af, e.af=e.af))
 stop("Failed gwdsp term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed gwdsp term test")
}



# gwesp, either
num.tests=num.tests+1
s.0 <- summary(fmh~gwesp)
e.0 <- ergm(samplike~gwesp(fixed=TRUE), estimate="MPLE")
e.a <- ergm(samplike~gwesp(.8, fixed=TRUE), estimate="MPLE")
s.f <- summary(fmh~gwesp(fixed=TRUE))
s.af <- summary(fmh~gwesp(.3, fixed=TRUE))
e.af <- ergm(samplike~gwesp(.2, fixed=TRUE), estimate="MPLE")
if (!all(head(s.0)==c(70,36,13,0,1,0)) ||
    round(e.0$coef + .4115515, 3) != 0 ||
    round(e.a$coef + .1898684, 3) != 0 ||
    s.f!=120 ||
    round(s.af - 133.9215, 3) != 0 ||
    round(e.af$coef + .3371385, 3) != 0) {
  print(list(s.0=head(s.0), e.0=e.0, e.a=e.a, s.f=s.f, s.af=s.af, e.af=e.af))
 stop("Failed gwesp term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed gwesp term test")
}



# gwnsp, either
num.tests=num.tests+1
s.0 <- summary(fmh~gwnsp)
e.0 <- ergm(samplike~gwnsp(fixed=TRUE), estimate="MPLE")
e.a <- ergm(samplike~gwnsp(.8, fixed=TRUE), estimate="MPLE")
s.f <- summary(fmh~gwnsp(fixed=TRUE))
s.af <- summary(fmh~gwnsp(.3, fixed=TRUE))
e.af <- ergm(samplike~gwnsp(.2, fixed=TRUE), estimate="MPLE")
if (!all(head(s.0)==c(361,39,10,1,0,0)) ||
    round(e.0$coef + .4189, 3) != 0 ||
    round(e.a$coef + .3123, 3) != 0 ||
    s.f!=411 ||
    round(s.af - 424.7154, 3) != 0 ||
    round(e.af$coef + .3934841, 3) != 0) {
  print(list(s.0=head(s.0), e.0=e.0, e.a=e.a, s.f=s.f, s.af=s.af, e.af=e.af))
 stop("Failed gwnsp term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed gwnsp term test")
}



# hamming, any
num.tests=num.tests+1
mat.d <- matrix(0,18,18)
mat.u <- matrix(0, 205, 205)
set.seed(456)
# Using a covariate matrix that matches the edges exactly is too easy.
cov.d <- cbind(as.edgelist(samplike)[,2:1], rbinom(88, 3, .5))
set.seed(145)
cov.u <- cbind(as.edgelist(fmh), rbinom(203, 3, .5))

# although there are 4 non-required inputs, giving
# 16 combinations of inputs, I've exlcuded most that
# don't involve 'x' because w/o 'x', the results are
# 0 or largely negative, as the hamming distance is 
# compared between identical networks
s.0 <- 0# COMMENTED OUT FOR NOW BECAUSE IT'S BROKEN:  summary(samplike~hamming)
s.x <- summary(samplike~hamming(mat.d))
# and everything commented below is broke.

# should this really be NA
#e.x <- ergm(fmh~hamming(mat.u), estimate="MPLE")
## OK
s.xc <- summary(samplike~hamming(mat.d, cov=cov.d))
# NA
#e.xc <- ergm(fmh~hamming(mat.u, cov=cov.u), estimate="MPLE")
# OK
s.xd <- summary(samplike~hamming(mat.d, defaultweight=.3))
# NA value
#e.xd <- ergm(samplike~hamming(mat.d, defaultweight=.3), estimate="MPLE")
# OK
s.xca <- summary(samplike~hamming(mat.d, cov=samplike, attrname="YearsTrusted"))
# NA
#e.xca <- ergm(fmh~hamming(mat.u, cov=fmh, attrname="Grade"), estimate="MPLE")
# OK
s.xcd<- summary(samplike~hamming(mat.d, cov=cov.d, defaultweight=.5))
# NA
#e.xcd<- ergm(samplike~hamming(mat.d, cov=cov.d, defaultweight=.5), estimate="MPLE")
# 0 & NA
#s.xcad<- summary(samplike~hamming(mat.d, samplike, "YearsTrusted", .5))
#e.xcad<- ergm(samplike~hamming(mat.d, samplike, "YearsTrusted", .5), estimate="MPLE")
if (FALSE && !all.equal(as.vector(c(s.0, s.x, s.xc, s.xd, s.xca, s.xcd)),
               as.vector(c(  0,  88,   84, 26.4,   183,   100)))) {
 print(list(s.0=s.0, s.x=s.x, s.xc=s.xc, s.xd=s.xd, s.xca=s.xca, s.xcd=s.xcd))
 stop("Failed hamming term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed hamming term test")
}


                     
# isolates, either
num.tests=num.tests+1
s.0 <- summary(samplike~isolates)
e.0 <- ergm(fmh~isolates, estimate="MPLE")
if (s.0 != 0 || round(e.0$coef - 5.10979, 3) != 0) {
 print(list(s.0=s.0, e.0=e.0))
 stop("Failed isolates term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed isolates term test")
}


                     
## localtriangle, either
#num.tests=num.tests+1
#set.seed(85)
#x <- matrix(rbinom(324, 2, .5),18,18)
#s.x <- summary(samplike~localtriangle(x))
#e.x <- ergm(samplike~localtriangle(x), estimate="MPLE")
#s.xa <- summary(fmh~localtriangle(fmh, "GradeMet"))
#if (s.x != 56 || round(e.x$coef + .1553, 3) != 0 ||
#    s.xa != 61) {
# print(list(s.x=s.x, e.x=e.x, s.xa=s.xa))
# stop("Failed localtriangle term test")
#} else {
#  num.passed.tests=num.passed.tests+1
#  print("Passed localtriangle term test")
#}


                
# meandeg, either
num.tests=num.tests+1
s.0 <- summary(samplike~meandeg)
e.0 <- ergm(fmh~meandeg, estimate="MPLE")
if (round(s.0 - 9.77777, 3) != 0 ||
    round(e.0$coef + 474.0647, 3) != 0) {
 print(list(s.0=s.0, e.0=e.0))
 stop("Failed meandeg term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed meandeg term test")
}


                
#nodecov, either
num.tests=num.tests+1
s.a <- summary(samplike~nodecov("YearsServed"))
e.a <- ergm(fmh~nodecov("Grade"), estimate="MPLE")
s.at <- summary(samplike~nodecov("YearsServed", function(x)x^2))
e.at <- ergm(fmh~nodecov("Grade", function(x)x^2), estimate="MPLE")
s.att <- summary(samplike~nodecov("YearsServed", function(x)x^2, "squared"))
if (s.a != 906 || round(e.a$coef + .271, 3) != 0 ||
    s.at != 5036 || round(e.at$coef + .03199, 3) != 0) {
 print(list(s.a=s.a, e.a=e.a, s.at=s.at, e.at=e.at))
 stop("Failed nodecov term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed nodecov term test")
}



#nodefactor, either
num.tests=num.tests+1
s.a <- summary(fmh~nodefactor("Grade"))
e.a <- ergm(samplike~nodefactor("group"), estimate="MPLE")
s.ab <- summary(fmh~nodefactor("Sex", base=4:5))
e.ab <- ergm(samplike~nodefactor("Trinity", base=0), estimate="MPLE")
if (!all(s.a==c(75, 65, 36, 49, 28)) ||
    !all(round(e.a$coef+c(.9480, .3273),3)==0) ||
    !all(s.ab==c(235,171)) ||
    !all(round(e.ab$coef+c(.4451, .4451, .4706),3)==0)) {
  print(list(s.a=s.a,e.a=e.a, s.ab=s.ab, e.ab=e.ab))
  stop("Failed nodefactor term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed nodefactor term test")
}



#nodematch, either
num.tests=num.tests+1
s.a <- summary(fmh~nodematch("Race"))
e.a <- ergm(samplike~nodematch("Trinity"), estimate="MPLE")
s.ad <- summary(samplike~nodematch("group", diff=TRUE))
e.ad <- ergm(fmh~nodematch("Sex", diff=TRUE), estimate="MPLE")
s.ak <- summary(fmh~nodematch("Grade", keep=3:4))
e.ak <- ergm(samplike~nodematch("group", keep=2), estimate="MPLE")
s.adk <- summary(samplike~nodematch("Trinity", TRUE, 1:2))
e.adk <- ergm(fmh~nodematch("Race", TRUE, 2), estimate="MPLE")
if (s.a != 103 || round(e.a$coef + 1.45725,3)!=0  ||
    !all(s.ad==c(23,10,30)) ||
    !all(round(e.ad$coef+c(4.06317, 4.7032),3)==0) ||
    s.ak!=32 || !all(round(e.ad$coef+c(4.063173, 4.703204),3)!=0 ||
    !all(s.adk==c(8,4)) ||
    round(e.adk$coef+ 4.700995,3)==0)) {
 print(list(s.0=s.0, e.0=e.0, s.a=s.a, e.a=e.a, s.ad=s.ad, e.ad=e.ad, s.ak=s.ak,
            e.ak=e.ak, s.adk=e.adk))
 stop("Failed nodematch term test")
} else {
  num.passed.tests=num.passed.tests+1
 print("Passed nodematch term test")
}

                
# nodemix, any
num.tests=num.tests+1
s.a <- summary(fmh ~ nodemix("Grade"))
e.a <- ergm(samplike ~ nodemix("group"), estimate="MPLE")
s.ab <- summary(bipnw ~ nodemix("Letter"), base=0)
e.ab <- ergm(bipnw ~ nodemix("Letter", base=2:6))
s.ab2 <- summary(fmh ~ nodemix("Race", base=1))
e.ab2 <- ergm(samplike ~ nodemix("Trinity", base=3:9))                
if (!all(s.a == c(75, 0, 33, 0, 2, 23, 1, 4, 7, 9, 1,
                  2, 6, 1, 17, 1, 1, 4, 5, 5, 6)) ||
    !all(round(e.a$coef - c(0.1910552, -3.2958369, -2.1747517, -2.5649494,
                           1.6094379, -3.2958369, -1.4916549, -1.0986123,
                            0.9162907), 3) == 0) ||
    !all(s.ab==c(9,8,8,7,7,5,4,6,6)) ||
    !all(round(e.ab$coef+c(3.497, 4.431, 3.989, 3.989),3)==0) ||
    !all(s.ab2==c(8,53,13,41,46,0,1,0,0,5,22,10,0,4)) ||
    !all(round(e.ab2$coef+c(1.0116, .82098),3)==0)) {
  print(list(s.a=s.a, e.a=e.a, s.ab=s.ab, e.ab=e.ab, s.ab2=s.ab2, e.ab2=e.ab2))
  stop("Failed nodemix term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed nodemix term test")
}



# nsp, either
num.tests=num.tests+1
s.d <- summary(fmh~nsp(2:3))
e.d <- ergm(samplike~nsp(4), estimate="MPLE")
if (!all(s.d==c(39, 10)) ||
   round(e.d$coef + 1.1096, 3) != 0) {
   print(list(s.d=s.d, e.d=e.d))
   stop("Failed nsp term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed nsp term test")
}

                
                
# smalldiff
num.tests=num.tests+1
s.ac.d <- summary(samplike~smalldiff("YearsServed", 3))
s.ac.u <- summary(fmh~smalldiff("Grade", 2))
s.ac.b <- summary(bipnw~smalldiff("Cost", 1))                
e.ac.d <- ergm(samplike~smalldiff("YearsServed", 3), estimate="MPLE")
e.ac.u <- ergm(fmh~smalldiff("Grade", 2), estimate="MPLE")
e.ac.b <- ergm(bipnw~smalldiff("Cost", 1), estimate="MPLE")                
if (s.ac.d != 78 || s.ac.u != 193 || s.ac.b != 48 ||
    round(e.ac.d$coef + .86903, 3) != 0 ||
    round(e.ac.u$coef + 4.3525, 3) != 0 ||
    round(e.ac.b$coef + 3.8318, 3) != 0 ) {
print(list(s.ac.d=s.ac.d, s.ac.u=s.ac.u, s.ac.b=s.ac.b,
           e.ac.d=e.ac.d, e.ac.u=e.ac.u, e.ac.b=e.ac.b))
 stop("Failed smalldiff term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed smalldiff term test")
}


                
# threepath, either
num.tests=num.tests+1
s.0 <- summary(samplike~threepath)
e.0 <- ergm(fmh~threepath, estimate="MPLE")
s.k <- summary(samplike~threepath(keep=2))
e.k <- ergm(samplike~threepath(keep=1:2), estimate="MPLE")
if (!all(s.0==c(2103, 2326, 1749, 1897)) ||
    round(e.0$coef + .2842, 3) != 0 ||
    s.k!=2326 ||
    !all(round(e.k$coef+c(.0188, -.0077),3)==0)) {
 print(list(s.0=s.0, e.0=e.0, s.k=s.k, e.k=e.k))
 stop("Failed threepath term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed threepath term test")
}

                

# triangles, either
num.tests=num.tests+1
s.0 <- summary(fmh~triangles)
e.0 <- ergm(samplike~triangles, estimate="MPLE")
s.a <- summary(fmh~triangles("Race"))
e.a <- ergm(samplike~triangle("group"), estimate="MPLE")                
s.ad <- summary(samplike~triangles("Trinity", diff=TRUE))
e.ad <- ergm(fmh~triangle("Sex", diff=TRUE), estimate="MPLE")   
if (s.0 != 62 || round(e.0$coef + .06997, 3) != 0 ||
    s.a != 18 || round(e.a$coef - .06354, 3) != 0 ||
    !all(s.ad==c(2,0,0)) ||
    !all(round(e.ad$coef + c(.70278, .44099), 3) == 0)) { 
 print(list(s.0=s.0, e.0=e.0, s.a=s.a, e.a=e.a, s.ad=s.ad, e.ad=e.ad))
 stop("Failed triangles term test")
} else {
  num.passed.tests=num.passed.tests+1
 print("Passed triangles term test")
}



                
# triadcensus, either
num.tests=num.tests+1
s.0 <- summary(samplike~triadcensus)
e.0 <- ergm(fmh~triadcensus, estimate="MPLE")
s.d <- summary(samplike~triadcensus(3))
e.d <- ergm(fmh~triadcensus(2:3), estimate="MPLE")
if (!all(s.0==c(205, 190, 12, 24, 24, 68, 34, 5, 0, 35, 15, 6, 5, 18, 8)) ||
    !all(round(e.0$coef+c(.02559, .06254, -2.61531),3)==0) ||
    s.d != 12 || !all(round(e.d$coef - c(-1.749635, 2.228183), 3) ==0)) {
 print(list(s.0=s.0, e.0=e.0, s.d=s.d, e.d=e.d))
 stop("Failed triadcensus term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed triadcensus term test")
}


                
# twopath, either
num.tests=num.tests+1
s.0 <- summary(samplike~twopath)
e.0 <- ergm(fmh~twopath, estimate="MPLE")
if (s.0 != 378 || round(e.0$coef + 1.297362, 3)){
 print(list(s.0=s.0, e.0=e.0))
 stop("Failed twopath term test")
} else {
  num.passed.tests=num.passed.tests+1
  print("Passed twopath term test")
}


if (num.passed.tests==num.tests)
  print("Passed all flexible term tests")
