library(ergm)
data(florentine)
f.miss<-network.copy(flomarriage)
f.miss[2,1] <- NA
f.miss[3,1] <- NA
fit <- ergm(f.miss ~ edges)
OK.2.1<-OK.3.1<-FALSE
for(s in 1:1000){
  sim <- simulate(fit, constraints=~observed)
  d<-(sim[,] != flomarriage[,])
  if(d[2,1]) OK.2.1<-TRUE
  if(d[3,1]) OK.3.1<-TRUE

  d[2,1]<-d[3,1]<-d[1,2]<-d[1,3]<-FALSE
  if(any(d)) stop("Observed dyad toggled!")

  if(OK.2.1 && OK.3.1) break
}

if(!OK.2.1 || !OK.3.1) stop("No toggles of a missing dyad in 1000 trials.")



