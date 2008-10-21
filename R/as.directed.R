#  File ergm/R/as.directed.R
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
as.directed<-function(x){
  if(!is.network(x))
    stop("as.directed requires an argument of class network.\n")
  else{
    if(get.network.attribute(x,"directed")){
     newmatrix <- as.matrix.network(x,matrix.type="edgelist")
     n1 <- network.size(x)+1
     eid <- cbind(pmax(newmatrix[,1],newmatrix[,2]),
                  pmin(newmatrix[,1],newmatrix[,2]))
     newmatrix <- eid[,1]+n1*eid[,2]
     newmatrix <- sort(unique(newmatrix))   
     eid <- trunc(newmatrix/n1)
     newmatrix <- cbind(newmatrix-eid*n1,eid)
     unw <- network.copy(x)
     eid<-vector()
     for(i in 1:network.size(x)){  
       eid <- c(eid,get.edgeIDs(unw,i))
     }
     delete.edges(unw,eid)
     unw %n% "directed" <- FALSE
     if(nrow(newmatrix)>0){
      add.edges(unw,head=newmatrix[,2],tail=newmatrix[,1])
     }
     unw
    }else{
     return(x)
    }
  }
}
