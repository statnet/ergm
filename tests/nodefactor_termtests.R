library(ergm)
data(faux.mesa.high)
data(sampson)

#nodefactor
s <- summary(faux.mesa.high~nodefactor("Grade", base=0))
fmh <- faux.mesa.high
fmh[146,192]<-0 # eliminate only "Other"-race tie to test drop=TRUE
e<-ergm(fmh ~ nodefactor("Race",base=0)+nodefactor("Sex")+nodefactor("Grade",base=3:4))
if (!all(s==c(153, 75, 65, 36, 49, 28)) ||
    !all(round(e$coef[!is.infinite(e$coef)]-c(-1.500, -2.613, -2.203,  -2.163, 
         -0.386, 0.488, 0.181, 0.244, 0.381),3)==0) ||
    e$coef[4]!=-Inf) {
  print(list(s=s,e=e))
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
if (!all(s==c(0, 6, 9, 15, 0, 0)) ||
    !all(round(e$coef[!is.infinite(e$coef)]-c(-1.299, -1.492),3)==0) ||
    e$coef[3]!=-Inf || e$coef[4]!=-Inf) {
  print(list(s=s,e=e))
  stop("Failed nodeifactor and nodeofactor test")
} else {
  print("Passed nodeifactor and nodeofactor test")
}


