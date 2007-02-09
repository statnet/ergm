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
