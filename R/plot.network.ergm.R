#  File R/plot.network.ergm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################


#' Two-Dimensional Visualization of Networks
#'
#' \code{\link{plot.network.ergm}} produces a simple two-dimensional plot of
#' the network object \code{x}. A variety of options are available to control
#' vertex placement, display details, color, etc.  The function is based on the
#' plotting capabilities of the \code{\link[network]{network}} package with
#' additional pre-processing of arguments.  Some of the capabilites require the
#' \code{\link[latentnet]{latentnet}} package.  See
#' \code{\link[network]{plot.network}} in the \code{\link[network]{network}}
#' package for details.
#'
#' \code{\link[network]{plot.network}} is a version of the standard network
#' visualization tool within the \code{sna} package.  By means of clever
#' selection of display parameters, a fair amount of display flexibility can be
#' obtained.  Network layout -- if not specified directly using \code{coord} --
#' is determined via one of the various available algorithms.  These are
#' (briefly) as follows: \enumerate{
#' 
#' \item \code{latentPrior}: Use a two-dimensional latent space model based on
#' a Bayesian minimum Kullback-Leibler fit.  See documentation for
#' \code{latent()} in \code{\link{ergm}}.
#' 
#' \item \code{random}: Vertices are placed (uniformly) randomly within a
#' square region about the origin.
#' 
#' \item \code{circle}: Vertices are placed evenly about the unit circle.
#' 
#' \item \code{circrand}: Vertices are placed in a "Gaussian donut," with
#' distance from the origin following a normal distribution and angle relative
#' to the X axis chosen (uniformly) randomly.
#' 
#' \item \code{eigen}, \code{princoord}: Vertices are placed via (the real
#' components of) the first two eigenvectors of: \enumerate{ \item
#' \code{eigen}: the matrix of correlations among (concatenated) rows/columns
#' of the adjacency matrix
#' 
#' \item \code{princoord}: the raw adjacency matrix.  }
#' 
#' \item \code{mds}, \code{rmds}, \code{geodist}, \code{adj}, \code{seham}:
#' Vertices are placed by a metric MDS.  The distance matrix used is given by:
#' \enumerate{ \item \code{mds}: absolute row/column differences within the
#' adjacency matrix
#' 
#' \item \code{rmds}: Euclidean distances between rows of the adjacency matrix
#' 
#' \item \code{geodist}: geodesic distances between vertices within the network
#' 
#' \item \code{adj}: \eqn{(\max A)-A}{(max A)-A}, where \eqn{A}{A} is the raw
#' adjacency matrix
#' 
#' \item \code{seham}: structural (dis)equivalence distances (i.e., as per
#' \code{sedist} in the package \code{sna}) based on the Hamming metric }
#' 
#' \item \code{spring}, \code{springrepulse}: Vertices are placed using a
#' simple spring embedder.  Parameters for the embedding model are given by
#' \code{embedder.params}, in the following order: vertex mass; equilibrium
#' extension; spring coefficient; repulsion equilibrium distance; and base
#' coefficient of friction.  Initial vertex positions are in random order
#' around a circle, and simulation proceeds -- increasing the coefficient of
#' friction by the specified base value per unit time -- until "motion"
#' within the system ceases.  If \code{springrepulse} is specified, then an
#' inverse-cube repulsion force between vertices is also simulated; this force
#' is calibrated so as to be exactly equal to the force of a unit spring
#' extension at a distance specified by the repulsion equilibrium distance.  }
#' 
#' @param x an object of class \code{\link[network]{network}}.
#' @param attrname an optional edge attribute, to be used to set edge values.
#' @param label a vector of vertex labels, if desired; defaults to the vertex
#' labels returned by \code{\link[network]{network.vertex.names}}.
#' @param coord user-specified vertex coordinates, in an NCOL(dat)x2 matrix.
#' Where this is specified, it will override the \code{mode} setting.
#' @param jitter boolean; should the output be jittered?
#' @param thresh real number indicating the lower threshold for tie values.
#' Only ties of value >\code{thresh} are displayed.  By default,
#' \code{thresh}=0.
#' @param usearrows boolean; should arrows (rather than line segments) be used
#' to indicate edges?
#' @param mode the vertex placement algorithm; this must correspond to a
#' \code{network.layout} function.  These include \code{"latent"},
#' \code{"latentPrior"}, and \code{"fruchtermanreingold"}.
#' @param displayisolates boolean; should isolates be displayed?
#' @param interactive boolean; should interactive adjustment of vertex
#' placement be attempted?
#' @param xlab x axis label.
#' @param ylab y axis label.
#' @param xlim the x limits (min, max) of the plot.
#' @param ylim the y limits of the plot.
#' @param pad amount to pad the plotting range; useful if labels are being
#' clipped.
#' @param label.pad amount to pad label boxes (if \code{boxed.labels==TRUE}),
#' in character size units.
#' @param displaylabels boolean; should vertex labels be displayed?
#' @param boxed.labels boolean; place vertex labels within boxes?
#' @param label.pos position at which labels should be placed, relative to
#' vertices.  \code{0} results in labels which are placed away from the center
#' of the plotting region; \code{1}, \code{2}, \code{3}, and \code{4} result in
#' labels being placed below, to the left of, above, and to the right of
#' vertices (respectively); and \code{label.pos>=5} results in labels which are
#' plotted with no offset (i.e., at the vertex positions).
#' @param label.bg background color for label boxes (if
#' \code{boxed.labels==TRUE}); may be a vector, if boxes are to be of different
#' colors.
#' @param vertex.sides number of polygon sides for vertices; may be given as a
#' vector or a vertex attribute name, if vertices are to be of different types.
#' @param vertex.rot angle of rotation for vertices (in degrees); may be given
#' as a vector or a vertex attribute name, if vertices are to be rotated
#' differently.
#' @param arrowhead.cex expansion factor for edge arrowheads.
#' @param label.cex character expansion factor for label text.
#' @param loop.cex expansion factor for loops; may be given as a vector or a
#' vertex attribute name, if loops are to be of different sizes.
#' @param vertex.cex expansion factor for vertices; may be given as a vector or
#' a vertex attribute name, if vertices are to be of different sizes.
#' @param edge.col color for edges; may be given as a vector, adjacency matrix,
#' or edge attribute name, if edges are to be of different colors.
#' @param label.col color for vertex labels; may be given as a vector or a
#' vertex attribute name, if labels are to be of different colors.
#' @param vertex.col color for vertices; may be given as a vector or a vertex
#' attribute name, if vertices are to be of different colors.
#' @param label.border label border colors (if \code{boxed.labels==TRUE}); may
#' be given as a vector, if label boxes are to have different colors.
#' @param vertex.border border color for vertices; may be given as a vector or
#' a vertex attribute name, if vertex borders are to be of different colors.
#' @param edge.lty line type for edge borders; may be given as a vector,
#' adjacency matrix, or edge attribute name, if edge borders are to have
#' different line types.
#' @param label.lty line type for label boxes (if \code{boxed.labels==TRUE});
#' may be given as a vector, if label boxes are to have different line types.
#' @param vertex.lty line type for vertex borders; may be given as a vector or
#' a vertex attribute name, if vertex borders are to have different line types.
#' @param edge.lwd line width scale for edges; if set greater than 0, edge
#' widths are scaled by \code{edge.lwd*dat}.  May be given as a vector,
#' adjacency matrix, or edge attribute name, if edges are to have different
#' line widths.
#' @param label.lwd line width for label boxes (if \code{boxed.labels==TRUE});
#' may be given as a vector, if label boxes are to have different line widths.
#' @param edge.len if \code{uselen==TRUE}, curved edge lengths are scaled by
#' \code{edge.len}.
#' @param edge.curve if \code{usecurve==TRUE}, the extent of edge curvature is
#' controlled by \code{edge.curv}.  May be given as a fixed value, vector,
#' adjacency matrix, or edge attribute name, if edges are to have different
#' levels of curvature.
#' @param edge.steps for curved edges (excluding loops), the number of line
#' segments to use for the curve approximation.
#' @param loop.steps for loops, the number of line segments to use for the
#' curve approximation.
#' @param object.scale base length for plotting objects, as a fraction of the
#' linear scale of the plotting region. Defaults to 0.01.
#' @param uselen boolean; should we use \code{edge.len} to rescale edge
#' lengths?
#' @param usecurve boolean; should we use \code{edge.curve}?
#' @param suppress.axes boolean; suppress plotting of axes?
#' @param vertices.last boolean; plot vertices after plotting edges?
#' @param new boolean; create a new plot?  If \code{new==FALSE}, vertices and
#' edges will be added to the existing plot.
#' @param layout.par parameters to the \code{network.layout} function specified
#' in \code{mode}.
#' @param cex.main Character expansion for the plot title.
#' @param cex.sub Character expansion for the plot sub-title.
#' @param seed Integer for seeding random number generator.  See
#' \code{\link[base]{set.seed}}.
#' @param latent.control A list of parameters to control the \code{latent} and
#' \code{latentPrior} models, \code{dyadsample} determines the size above which
#' to sample the latent dyads; see \code{\link{ergm}} and
#' \code{\link[stats]{optim}} for details.
#' @param colornames A vector of color names that can be selected by index for
#' the plot. By default it is \code{colors()}.
#' @param verbose logical; if this is \code{TRUE}, we will print out more
#' information as we run the function.
#' @param latent logical; use a two-dimensional latent space model based on the
#' MLE fit.  See documentation for \code{ergmm()} in
#' \code{\link[latentnet]{latentnet}}.
#' @param \dots additional arguments to \code{\link{plot}}.
#' @return None.
#' @section Requires: \code{mva}
#' @author Carter T. Butts \email{buttsc@@uci.edu}
#' @seealso \code{\link{plot}}
#' @references Wasserman, S., and Faust, K.  (1994).  "Social Network
#' Analysis: Methods and Applications." Cambridge: Cambridge University Press.
#' @keywords hplot graphs
#' @examples
#' 
#' data(florentine)
#' plot(flomarriage)  #Plot the Florentine Marriage data
#' plot(network(10))  #Plot a random network
#' \dontrun{plot(flomarriage,interactive="points")}
#' 
#' @export plot.network.ergm
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
     trycol <- is.null(get.edge.value(x,edge.col))
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
