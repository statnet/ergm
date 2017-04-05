#  File ergm/R/ergm.geodist.R
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
# Copyright 2003 Mark S. Handcock, University of Washington
#                David R. Hunter, Penn State University
#                Carter T. Butts, University of California - Irvine
#                Steven M. Goodreau, University of Washington
#                Martina Morris, University of Washington
# Copyright 2007 The statnet Development Team
######################################################################
#  File ergm/R/ergm.geodist.R
#
# Carter's code from SNA
#
#geodist - Find the numbers and lengths of geodesics among nodes in a graph 
#using a BFS, a la Brandes (2000).  (Thanks, Ulrik!)
ergm.geodist<-function(dat,inf.replace=Inf){
   n<-dim(dat)[2]
   #Initialize the matrices
   sigma<-matrix(0,nrow=n,ncol=n)
   gd<-matrix(Inf,nrow=n,ncol=n)
   #Perform the calculation
   geo<-.C("geodist_R",as.double(dat),as.double(n),gd=as.double(gd), sigma=as.double(sigma),NAOK=TRUE, PACKAGE='ergm')
   #Return the results
   o<-list()
   o$counts<-matrix(geo$sigma,n,n)
   o$gdist<-matrix(geo$gd,n,n)
   o$gdist[o$gdist==Inf]<-inf.replace  #Patch Infs, if desired
   o
}

fullgcount<-function(dat,geodist.precomp=NULL, directed=FALSE){
   #Get the counts matrix
   if(is.null(geodist.precomp))
      cnt<-ergm.geodist(dat)$counts
   else
      cnt<-geodist.precomp$counts
   if(directed){
     allcnt <- cnt[row(cnt)!=col(cnt)]
   }else{
     allcnt <- cnt[row(cnt)<col(cnt)]
   }
   tabulate(allcnt,nbins=ncol(cnt))
}

ergm.geodistdist<-function(nw, directed=is.directed(nw)){
 ergm.geodistn(edgelist=as.matrix.network(nw,matrix.type="edgelist"),
               n=nw$gal$n, directed=directed)/(2-is.directed(nw))
}

ergm.geodesicmatrix<-function(nw, directed=is.directed(nw)){
 ergm.geodesicmatrix.edgelist(edgelist=as.matrix.network(nw,matrix.type="edgelist"),
               n=nw$gal$n, directed=directed)
}

ergm.geodesicmatrix.old<-function(nw, directed=is.directed(nw), n=nw$gal$n){
  ans<-matrix(0,n,n)
  for(i in 1:n){
    ans[i,] <- ergm.nodegeodesics(
     edgelist=as.matrix.network(nw,matrix.type="edgelist"),i,n,directed)
  }
  ans[ans==n]<-Inf
  ans
}
