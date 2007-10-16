n<-100
base.net<-matrix(0,n,n)
base.net<-as.network.matrix(base.net,dir=FALSE,ma="a")
ergm.fit<-ergm(base.net~degree(c(0,1,2)),meanstats=n*c(.2,.4,.2))
