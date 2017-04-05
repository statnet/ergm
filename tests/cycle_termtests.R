library(ergm)
data(faux.mesa.high)
data(sampson)

fmh <- faux.mesa.high
s1 <- summary(fmh~cycle(3:6))-c(62,80,138,270)
s2 <- summary(samplike~cycle(2:6))-c(28,39,111,260,651)
if(!all(s1==0) || !all(s2==0)) {
 print(list(s1=s1,s2=s2))
 stop("Failed cycle test")
} else {
 print("Passed cycle test")
}


