#  File R/plot.network.ergm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2015 Statnet Commons
#######################################################################
#########################################################################
# The <plot.network.ergm> function produces a two-dimensional network
# visualization based on <plot.network.default>; a variety of options are 
# available to control vertex placement, display details, color, etc; the 
# function is based on the plotting capabilities of the network package 
# with additional pre-processing of arguments; some of the capabilites 
# require the latentnet package; see <plot.network> in the network package
# for details.
#
# --PARAMETERS--
#   x: a network
#    : (see the man page for descriptions of the other 52 input params)
#
# --RETURNED--
#   the plot as an invisible list containing:
#    x        : the x-coordinates used in the plot
#    y        : the y-coordinates used in the plot
#    latentfit: the latent fit as a list containing:
#       Z.mle : ??
#       beta  : ??
#
#########################################################################

"plot.network.ergm"<-function(x,
    attrname=NULL,
    label=network.vertex.names(x),
    coord=NULL,
    jitter=TRUE,
    thresh=0,
    usearrows=TRUE,
    mode="fruchtermanreingold",
    displayisolates=TRUE,
    interactive=FALSE,
    xlab=NULL,                                      
    ylab=NULL,
    xlim=NULL,
    ylim=NULL,
    pad=0.2,
    label.pad=0.5,
    displaylabels=FALSE,
    boxed.labels=TRUE,
    label.pos=0,
    label.bg="white",
    vertex.sides=8,
    vertex.rot=0,
    arrowhead.cex=1,
    label.cex=1,
    loop.cex=1,
    vertex.cex=1,
    edge.col=1,
    label.col=1,
    vertex.col=2,
    label.border=1,
    vertex.border=1,
    edge.lty=1,
    label.lty=NULL,
    vertex.lty=1,
    edge.lwd=0,
    label.lwd=par("lwd"),
    edge.len=0.5,
    edge.curve=0.1,
    edge.steps=50,
    loop.steps=20,
    object.scale=0.01,
    uselen=FALSE,
    usecurve=FALSE,
    suppress.axes=TRUE,
    vertices.last=TRUE,
    new=TRUE,
    layout.par=NULL,
    cex.main=par("cex.main"), 
    cex.sub=par("cex.sub"),
    seed=NULL,
    latent.control=list(maxit=500,trace=0,dyadsample=10000,
               penalty.sigma=c(5,0.5), nsubsample=200),
    colornames="rainbow",
    verbose=FALSE, latent=FALSE, ...){
#
   #Extract the network to be displayed
   if(is.hyper(x)){    #Is this a hypergraph?  If so, use two-mode form.
     d<-as.matrix.network(x,matrix.type="incidence",attrname=attrname)
     n<-sum(dim(d))
     temp<-matrix(0,nrow=n,ncol=n)
     if(is.directed(x)){  #If directed, depict as such.
       temp[1:dim(d)[1],(dim(d)[1]+1):n]<-abs(pmin(d,0))     #Tail set
       temp[(dim(d)[1]+1):n,1:dim(d)[1]]<-t(abs(pmax(d,0)))  #Head set
       d<-temp
     }else{
       temp[1:dim(d)[1],(dim(d)[1]+1):n]<-d
       temp[lower.tri(temp)]<-t(temp)[lower.tri(temp)]
       d<-temp
       usearrows<-FALSE   #Don't use labels for undirected networks
     }
     if(length(label)==network.size(x))  #Fix labels, if needed
       label<-c(label,paste("e",1:(n-network.size(x)),sep=""))
   }else{
     n<-network.size(x)
     d<-as.matrix.network(x,matrix.type="adjacency",attrname=attrname)
     if(!is.directed(x))
       usearrows<-FALSE
   }
   diag<-has.loops(x)         #Check for existence of loops
   #Replace NAs with 0s
   d[is.na(d)]<-0
#
#  Start ergm specific
#
# Get missingness matrix
#
   Ydesign <- is.na(x)  
   if(network.edgecount(Ydesign)!=0){
    if(is.network(Ydesign)){Ydesign <- as.sociomatrix(Ydesign)}
    Ydesign <- Ydesign==0
   }   
   current.warn <- options()$warn
   if(is.null(current.warn)){current.warn <- 0}
#
   if(!is.null(seed)) set.seed(seed)
#
#  "ergm" reads
#
#  xlims <- xlim; ylims <- ylim; rm("xlim","ylim")
   par("cex.main"=cex.main,"cex.sub"=cex.sub)
   if((length(label)==1) && is.logical(label) && label){
     displaylabels=TRUE
   }
   if((length(label.bg)==1) && is.character(label.bg)){
     if(is.na(match(label.bg,colors()))){
       label.bg <- as.numeric(as.factor(get.vertex.attribute(x,label.bg)))
     }
   }
   if(!missing(vertex.col) && length(vertex.col)==1 && is.character(vertex.col)){
     trycol <- as.numeric(as.factor(get.vertex.attribute(x,vertex.col)))
     if(!all(is.na(trycol))){
       trycol <- as.vector(unclass(factor(trycol)))
       if(is.numeric(trycol)){
        trycol[is.na(trycol)] <- max(trycol,na.rm=TRUE)+1
       }else{
        trycol[is.na(trycol)] <- "purple"
       }
       if(missing(colornames)){
         vertex.col <- rainbow(length(unique(as.vector(trycol))))[trycol]
       }else{
         vertex.col <- colornames[trycol]
       }
     }
   }
   if(missing(vertex.col) && is.bipartite(x)){
    nb1 <- get.network.attribute(x,"bipartite")
    nb2 <- network.size(x) - nb1
    if(missing(colornames)){
      vertex.col <- rep(rainbow(2),c(nb1, nb2))
    }else{
      vertex.col <- rep(colornames[1:2],c(nb1, nb2))
    }
   }
   if(length(vertex.col)==1){
     vertex.col <- rep(vertex.col,n)
   }
   if(!missing(vertex.sides) && length(vertex.sides)==1 && is.character(vertex.sides)){
     trycol <- unlist(get.vertex.attribute(x,vertex.sides))
     if(!any(is.na(trycol))){
       trycol <- c(8,4,2,6,10,3)[as.vector(unclass(factor(trycol)))]
       vertex.sides <- trycol
     }
   }
   if(missing(vertex.sides) && is.bipartite(x)){
    nb1 <- get.network.attribute(x,"bipartite")
    nb2 <- network.size(x) - nb1
    vertex.sides <- rep(c(8,3),c(nb1, nb2))
   }
   if(length(vertex.sides)==1){
     vertex.sides <- rep(vertex.sides,n)
   }
   if(!missing(vertex.cex) && length(vertex.cex)==1 && is.character(vertex.cex)){
     trycex <- unlist(get.vertex.attribute(x,vertex.cex))
     if(!any(is.na(trycex))){
       trycex <- trycex-min(trycex)
       vertex.cex <- 1.5*trycex/median(trycex)+0.5
     }
   }
   if(length(vertex.cex)==1){
     vertex.cex <- rep(vertex.cex,n)
   }
#  if(!missing(vertex.text) && length(vertex.text)==1 && is.character(vertex.text)){
#    trycol <- unlist(get.vertex.attribute(x,vertex.text))
#    if(!any(is.na(trycol))){
#      vertex.text <- paste(trycol)
#    }
#  }
#  if((!missing(label)|!missing(vertex.label)) 
   if((!missing(label)) 
       && length(label)==1 && is.character(label)){
     trycol <- unlist(get.vertex.attribute(x,label))
     if(!all(is.na(trycol))){
       label <- paste(trycol)
     }
   }else{
     trycol <- unlist(get.vertex.attribute(x,"vertex.names"))
     if(!all(is.na(trycol))){
       label <- paste(trycol)
     }
   }
#  if(!missing(edge.lwd) && length(edge.lwd)==1 && is.character(edge.lwd)){
#    trycol <- as.sociomatrix(x, attrname=edge.lwd)
#    if(!any(is.na(trycol))){
#      edge.lwd <- trycol
#    }
#  }
#  if(length(edge.lwd)==1){
#    edge.lwd <- matrix(edge.lwd,n,n)
#  }
   if(!missing(edge.col) && length(edge.col)==1 && is.character(edge.col)){
     trycol <- is.null(get.edge.attribute(x$mel,edge.col))
     if(!trycol){
       trycol <- as.sociomatrix(x, attrname=edge.col)
       ucols <- sort(unique(as.vector(trycol)))
       edge.col <- matrix(0,n,n)
       if(missing(colornames)){
         edgecol <- rainbow(length(ucols))
       }else{
         edgecol <- colornames
       }
       for(i in seq(along=ucols)){
        edge.col[trycol==ucols[i]] <- 
         edgecol[i-length(edgecol)*trunc(i/length(edgecol))]
       }
       edge.col[trycol==0] <- "white"
     }
   }
   if(length(edge.col)==1){
     edge.col <- matrix(edge.col,n,n)
   }
#
#  End "ergm" reads
#
   #Which nodes should we use?
   if(is.bipartite(x)){
     use<-displayisolates|(apply(d,1,sum)>0)
   }else{
     use<-displayisolates|((apply(d,1,sum)+apply(d,2,sum))>0)
   }
#
#  Which nodes should we use?
#
   if(sum(use)<3){
    if(all(!use)){warning("The network is empty and was not plotted.")}
    else{warning("The network has only one tie and was not plotted.")}
    options(warn=current.warn)
    return(invisible(x))
   }
   if(is.bipartite(x)){
     d <- d[use,]
   }else{
     d <- d[use,use]
   }
#
   if(!is.null(coord) | mode=="fruchtermanreingold" | interactive){
      latent <- FALSE
   }
   if(latent){
    latent <- try(requireNamespace('latentnet'))
   }
   if(is.list(coord)){
     coord <- cbind(coord$x, coord$y)
   }

   latentfit <- NULL
   if(latent){ #Place by latent space MLE
    requireNamespace('latentnet')
    plotfile <- paste("gplottmp",Sys.getpid(),sep="")
    curplot <- dev.cur()
    pictex(file = plotfile)
    options(warn=-1)
    coord <- plot.network.default(x,
         attrname=attrname,
         label=label,
         coord=coord,
         jitter=jitter,
         thresh=thresh,
         usearrows=usearrows,
         mode=mode,
         displayisolates=displayisolates,
         interactive=interactive,
         xlab=xlab,
         ylab=ylab,
         xlim=xlim,
         ylim=ylim,
         pad=pad,
         label.pad=label.pad,
         displaylabels=displaylabels,
         boxed.labels=boxed.labels,
         label.pos=label.pos,
         label.bg=label.bg,
         vertex.sides=vertex.sides,
         vertex.rot=vertex.rot,
         arrowhead.cex=arrowhead.cex,
         label.cex=label.cex,
         loop.cex=loop.cex,
         vertex.cex=vertex.cex,
         edge.col=edge.col,
         label.col=label.col,
         vertex.col=vertex.col,
         label.border=label.border,
         vertex.border=vertex.border,
         edge.lty=edge.lty,
         label.lty=label.lty,
         vertex.lty=vertex.lty,
         edge.lwd=edge.lwd,
         label.lwd=label.lwd,
         edge.len=edge.len,
         edge.curve=edge.curve,
         edge.steps=edge.steps,
         loop.steps=loop.steps,
         object.scale=object.scale,
         uselen=uselen,
         usecurve=usecurve,
         suppress.axes=suppress.axes,
         vertices.last=vertices.last,
         new=new,
         layout.par=layout.par,
         ...
     )
     dev.off()
     unlink(plotfile)
     dev.set(curplot)
     options(warn=current.warn)
#
#    Adjust for the components
#
     dimSpace <- 2
     if(verbose) cat("Calling geodesic distances\n")
     D <- ergm.geodesicmatrix(network(d,directed=is.directed(x)))
     D[D==Inf] <- max(D[D!=Inf])+1
     reach <- D!=max(D)
#
#    sample
#
     if(!is.null(Ydesign)){   
       nsample <- sum(reach[Ydesign]) 
     }else{   
       nsample <- sum(reach)
     }   
     if(!is.null(latent.control$dyadsample) &&
         nsample > latent.control$dyadsample){
       dyadsample <- sample(nsample,size=latent.control$dyadsample,replace=FALSE)
       samreach <- reach & !reach     
       if(!is.null(Ydesign)){   
         samreach[reach[Ydesign]][dyadsample] <- TRUE  
       }else{   
         samreach[reach][dyadsample] <- TRUE
       }   
     }else{
       samreach <- reach
       if(!is.null(Ydesign)){   
         samreach[!Ydesign] <- FALSE
       }   	
     }

     Y <- d[use,use]
     g <- nrow(Y)
     if(is.null(latent.control$MCMLE.maxit)){latent.control$MCMLE.maxit <- 40}
     if(is.null(latent.control$trace)){latent.control$trace <-6}
     if(is.null(latent.control$dyadsample)){latent.control$dyadsample <- 1000}
     if(is.null(latent.control$penalty.sigma)){latent.control$penalty.sigma <- c(10,0.5)}
     if(g > 100){
       maxit <- round(200/sqrt(g))
       trace <- 6
     }

     Z <- coord
     beta <- 0
     ## now find the mle
     abvZ <- c(beta,Z)

     if(latent.control$penalty.sigma[1]>0){
      penalty.factor <- c(
       1/(latent.control$penalty.sigma[1]*latent.control$penalty.sigma[1]),
       latent.control$penalty.sigma[2])
     }else{
      penalty.factor <- c(0,latent.control$penalty.sigma[2])
     }
     #BFGS  use previous found MDS fit to plug into quasi newton raphson
     abz.list <- list(Y=Y,nnodes=g,dimSpace=dimSpace,
                      penalty.factor=penalty.factor,
                      reach=samreach[use,use],directed=is.directed(x),
		      Ydesign=Ydesign[use,use])  
     if(verbose) cat("Calling latent MLE fit\n")
# Commented out by DH because mlpY.plot and mlpY.grad.plot do not exist.
#     MLE.fit <- try(
#                optim(par=abvZ,fn=mlpY.plot,gr=mlpY.grad.plot,
#                 method="BFGS",
#                 control=list(fnscale=-1, maxit=latent.control$MCMLE.maxit, 
#                              trace=latent.control$trace),
#                 abz.list=abz.list)
#                 )
#     if(inherits(MLE.fit,"try-error")){
      stop("MLE could not be found.")
#     }else{
#      abvZ <- MLE.fit$par
#     }
     Z.mle <- matrix(abvZ[-1],nrow=g,ncol=dimSpace)
     latentfit <- list(Z.mle=Z.mle, beta=abvZ[1])
     coord[use,1]<-Z.mle[,1]
     coord[!use,1]<-min(coord[use,1])
     coord[use,2]<-Z.mle[,2]
     coord[!use,2]<-min(coord[use,2])
#   End of latent
   }
   options(warn=-1)
    coord <- plot.network.default(x,
         attrname=attrname,
         label=label,
         coord=coord,
         jitter=jitter,
         thresh=thresh,
         usearrows=usearrows,
         mode=mode,
         displayisolates=displayisolates,
         interactive=interactive,
         xlab=xlab,
         ylab=ylab,
         xlim=xlim,
         ylim=ylim,
         pad=pad,
         label.pad=label.pad,
         displaylabels=displaylabels,
         boxed.labels=boxed.labels,
         label.pos=label.pos,
         label.bg=label.bg,
         vertex.sides=vertex.sides,
         vertex.rot=vertex.rot,
         arrowhead.cex=arrowhead.cex,
         label.cex=label.cex,
         loop.cex=loop.cex,
         vertex.cex=vertex.cex,
         edge.col=edge.col,
         label.col=label.col,
         vertex.col=vertex.col,
         label.border=label.border,
         vertex.border=vertex.border,
         edge.lty=edge.lty,
         label.lty=label.lty,
         vertex.lty=vertex.lty,
         edge.lwd=edge.lwd,
         label.lwd=label.lwd,
         edge.len=edge.len,
         edge.curve=edge.curve,
         edge.steps=edge.steps,
         loop.steps=loop.steps,
         object.scale=object.scale,
         uselen=uselen,
         usecurve=usecurve,
         suppress.axes=suppress.axes,
         vertices.last=vertices.last,
         new=new,
         layout.par=layout.par,
         ...
     )
   options(warn=current.warn)
#
#  Back to the original code
#
   invisible(list(x=coord[,1],y=coord[,2], latentfit=latentfit))
}
