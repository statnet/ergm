network.update<-function(nw,newmatrix,matrix.type=NULL,output="network")
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
  if(output=="as.edgelist.compressed") unw <- as.edgelist.compressed(unw)
  unw
}
#Force the input into edgelist form.  Network size, directedness, and vertex
#names are stored as attributes, since they cannot otherwise be included
#A copy of as.edgelist.sna
as.edgelist.compressed<-function(x, attrname=NULL, force.bipartite=FALSE){
  #In case of lists, process independently
  if(is.list(x)&&(!(class(x)%in%c("network"))))
    return(lapply(x,as.edgelist.compressed, attrname=attrname, force.bipartite=force.bipartite))
  #Begin with network objects
  if(class(x)=="network"){
    require("network")  #Must have network library to process network objects
    out<-as.matrix.network.edgelist(x,attrname=attrname)
#   if(!is.directed(x)){
#    out <- out[1:(nrow(x)/2),]
#   }
    if(NCOL(out)==2)                        #If needed, add edge values
      out<-cbind(out,rep(1,NROW(out)))
    attr(out,"n")<-network.size(x)
    attr(out,"directed")<-is.directed(x)
    attr(out,"vnames")<-network.vertex.names(x)
    van<-list.vertex.attributes(x)
    if(length(van)>0){
     va <- vector(mode = "list", length(van))
     for (i in (1:length(van))){ 
      va[[i]]<-get.vertex.attribute(x,van[i],unlist=TRUE)
     }
     names(va)<-van
     attr(out,"vertex.attributes")<-va
    }
    if(is.bipartite(x))
      attr(out,"bipartite")<-get.network.attribute(x,"bipartite")
    else if(force.bipartite)
      out<-as.edgelist.compressed(out,attrname=attrname,force.bipartite=force.bipartite)
  }else{
    warning("as.edgelist.compressed input must be network, or list thereof.\n Returning the original object.\n")
  }
  #Return the result
  out
}
as.network.uncompressed<-function(x, 
        na.rm=FALSE, edge.check=FALSE){
  #Initialize the network object
  if(class(x)=="network"){return(x)}
  if(is.null(attr(x,"vnames"))){
   warning("as.network.uncompressed input must be a compressed network, or a network.\n Returning the original object.\n")
   return(x)
  }
  n<-attr(x,"n")
  directed<-attr(x,"directed")
  g<-network.initialize(n,directed=directed)
  #Call the specific coercion routine, depending on matrix type
# g<-network.edgelist(x,g,na.rm=na.rm,edge.check=edge.check)
  g<-add.edges(g,as.list(x[,1]),as.list(x[,2]),edge.check=edge.check)
  va <- attr(x,"vertex.attributes")
  if(length(va)>0){
   for (i in (1:length(va))){ 
    g <- set.vertex.attribute(g,names(va)[i], va[[i]])
   }
  }
  #Return the result
  g
}
