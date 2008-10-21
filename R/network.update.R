#  File ergm/R/network.update.R
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
network.update<-function(nw,newmatrix,matrix.type=NULL)
{
#  print(paste("old:",network.edgecount(nw)," new:", nrow(newmatrix),collapse=" "))
  unw <- network.copy(nw)
  if(is.null(matrix.type)){
    warning("Don't leave matrix type to chance! Pass matrix.type to network.update!")
    matrix.type <- which.matrix.type(newmatrix)
    if(nrow(newmatrix)==0){matrix.type <- "edgelist"}
  }
  if(matrix.type=="adjacency" 
     && max(abs(newmatrix))==1 && max(abs(newmatrix-as.integer(newmatrix)))==0){
    unw[,] <- newmatrix
  }else if(matrix.type=="edgelist"){
#  cnw <- as.matrix.network(nw,matrix.type="edgelist")
#  unw[cnw[,2],cnw[,1]] <- 0
#  unw[,] <- 0
#  eid<-vector()
#  for(i in 1:network.size(nw)){  
#    eid <- c(eid,get.edgeIDs(unw,i))
#  }
   eid <- c(unlist(unw$iel),unlist(unw$oel))
   delete.edges(unw,eid)
   if(!is.null(newmatrix) && nrow(newmatrix)>0){
#   unw[newmatrix] <- 1
    add.edges(unw,head=newmatrix[,2],tail=newmatrix[,1])
   }
  }
  unw
}
