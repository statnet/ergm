fulldistdist<-function(dat,geodist.precomp=NULL, directed=FALSE){
   #Get the counts matrix
   if(is.null(geodist.precomp))
      cnt<-ergm.geodistR(dat)$gdist
   else
      cnt<-geodist.precomp$gdist
   if(directed){
     allcnt <- cnt[row(cnt)!=col(cnt)]
   }else{
     allcnt <- cnt[row(cnt)<col(cnt)]
   }
#  tabulate(allcnt,nbins=ncol(cnt))      
   allcnt[allcnt==Inf] <- ncol(cnt)-1
   tabulate(allcnt,nbins=ncol(cnt))
}
