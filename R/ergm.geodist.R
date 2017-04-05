#=======================================================================================
# This file contains the following 5 functions from Carter's SNA code
#          <ergm.geodist>            <ergm.geodesicmatrix>
#          <fullgcount>              <fulldistdist>
#          <ergm.geodesicmatrix>
#=======================================================================================




##############################################################################
# The <ergm.geodist> function finds the geodesic distance between all pairs of
# nodes via <geodist_R.C> using a BFS, a la Brandes (2000). (Thanks, Ulrik!)
# 
# --PARAMETERS--
#   dat        :  an adjacency matrix
#   inf.replace:  a replacement value for infinite values; default=Inf
#
# --RETURNED--
#   o: a list of two matrices:
#      sigma: the matrix of ?? between all pairs of nodes
#      gd   : the matrix of geodesic distances between all pairs of nodes
#
###############################################################################

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



ergm.geodesicmatrix<-function(nw, directed=is.directed(nw)){
 ergm.geodesicmatrix.edgelist(edgelist=as.matrix.network(nw,matrix.type="edgelist"),
               n=network.size(nw), directed=directed)
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



fulldistdist<-function(dat,geodist.precomp=NULL, directed=FALSE){
   #Get the counts matrix
   if(is.null(geodist.precomp))
      cnt<-ergm.geodist(dat)$gdist
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
