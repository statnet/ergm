#################################################################################
# The <plot.ergm> function does it plotting via the <mcmc.diagnostics> function.
# This function basically serves as a wrapper
#
# --PARAMETERS--
#   x: an ergm object
#   *: a host of parameters, all of which are ignored; for details see the
#      R documentation for <plot.ergm>
#
# --RETURNED--
#   NULL
# 
###############################################################################

"plot.ergm" <- function (x, ..., mle=FALSE, comp.mat = NULL,
            label = NULL, label.col = "black",
            xlab, ylab, main, label.cex = 0.8, edge.lwd = 1,
            edge.col=1, al = 0.1,
            contours=0, density=FALSE, only.subdens = FALSE, 
            drawarrows=FALSE,
            contour.color=1, plotnetwork=FALSE, pie = FALSE, piesize=0.07,
            vertex.col=1,vertex.pch=19,vertex.cex=2,
            mycol=c("black","red","green","blue","cyan",
              "magenta","orange","yellow","purple"),
            mypch=15:19,mycex=2:10)
{
  if(!missing(vertex.col) && length(vertex.col)==1 && is.character(vertex.col)){
    trycol <- as.numeric(as.factor(unlist(get.vertex.attribute(x$network,vertex.col))))
    if(!all(is.na(trycol))){
      vertex.col <- mycol[trycol]
    }
  }
  if((length(vertex.col)==1) && (vertex.col==1)) vertex.col <- mycol[x$class]

  if(!missing(label) && length(label)==1 && is.character(label)){
    trycol <- unlist(get.vertex.attribute(x$network,label))
    if(!all(is.na(trycol))){
      label <- trycol
    }
  }
  
  if(!missing(vertex.pch) && length(vertex.pch)==1 && is.character(vertex.pch)){
    trycol <- as.numeric(as.factor(unlist(get.vertex.attribute(x$network,vertex.pch))))
    if(!all(is.na(trycol))){
      vertex.pch <- mypch[trycol]
    }
  }
  
  if(!missing(vertex.cex) && length(vertex.cex)==1 && is.character(vertex.cex)){
    trycol <- as.numeric(as.factor(unlist(get.vertex.attribute(x$network,vertex.cex))))
    if(!all(is.na(trycol))){
      vertex.cex <- mycex[trycol]
    }
  }

  newnetwork <- as.sociomatrix(x$newnetwork)
#  if (!is.latent(x)) {    
#
#   So regular (non-latent) call
#
#   if(is.null(x$mcmc.loglikelihood))
#     cat("No plotting method available for non-latent models with no mcmc.loglikelihood\n")
#   else
#     plot(x$mcmc.loglikelihood[, 1], type = "l",
#          main = "MCMC log Likelihood Values", 
#          xlab = "sample", ylab = "log likelihood",...)
    if(is.null(x$sample))
      cat("No plotting method available for non-latent models with MCMC sample\n")
    else
      mcmc.diagnostics(x)
#  }
#  else {
##
## So latent call
##
#    if (missing(xlab)) 
#      xlab <- ""
#    if (missing(ylab)) 
#      ylab <- ""
#    if(mle)
#      {
#        z.pos <- x$Z.mle
#        if (missing(main)) 
#          main <- paste("MLEs of Latent Positions of", 
#                        deparse(substitute(x)))
#      }
#    else
#      {
#        z.pos <- x$Z.mkl
#        if (missing(main)) 
#          main <- paste("KL Latent Positions of", 
#                        deparse(substitute(x)))
#      }
#    if(density[1])
#    {
##
##   Plot densities
##
#      require(KernSmooth,quietly=TRUE)
#      if(density[1]>1)
#      {
##
##   Plot 2-dimensional densities
##
#        if(length(density)>1)
#          opar <- par(mfrow=c(density[1],density[2]),mar=c(2.5,2.5,1,1))
#        else
#          opar <- par(mfrow=c(density,density),mar=c(2.5,2.5,1,1))
#      }
#      if(!only.subdens)
#      {
#        Z.use <- cbind(as.vector(x$Z[,1,]),as.vector(x$Z[,2,]))
#        plot(Z.use,type='n',xlab=xlab,ylab=ylab,xlim=range(Z.use),
#             ylim=range(Z.use), ...)
#        title(main=paste("Posterior density of",main), cex.main=0.7, ...)
#        temp <- bkde2D(Z.use,0.2,c(201,201))
#        image(temp$x1,temp$x2,temp$fhat,col=grey(seq(1,0,length=255)),add=TRUE)
#        if(drawarrows)
#          for (i in 1:network.size(x$network))
#            for (j in 1:network.size(x$network))
#              if (newnetwork[i,j])
#                midarrow(x0 = z.pos[i, 1], y0 = z.pos[i, 2], x1 = z.pos[j,1],
#                         y1 = z.pos[j, 2], length = 0.1,col="yellow")
#        box()
#      }
#
#      if(density[1]>1)
#      {
#        Z.Ki.v <- as.vector(t(x$Ki))
#        Z.proc.mean <- cbind(as.vector(x$Z[,1,]),as.vector(x$Z[,2,]))
#        for(i in 1:x$ngroups)
#          {
#            plot(Z.proc.mean,
#                 xlim=range(Z.proc.mean),ylim=range(Z.proc.mean),
#                 main=paste("Class",i),
#                 type='n', ...)
#            temp <- bkde2D(Z.proc.mean[Z.Ki.v==i,],0.2,c(101,101))
#            use.col <- as.vector(col2rgb(mycol[i])/255)
#            image(temp$x1,temp$x2,temp$fhat,add=TRUE,
#                  col=rgb(seq(1,use.col[1],length=255),
#                    seq(1,use.col[2],length=255),
#                    seq(1,use.col[3],length=255)))
#            contour(temp$x1,temp$x2,temp$fhat,add=TRUE, nlevels=4,
#                    drawlabels=FALSE,
#                  col="white")
#            box()
#          }
#        if(plotnetwork)
#          plot(x,label=label)
#        par(opar)
#      }
#    }
#    else if(contours>0)
#    {
##do a contours x contours array of contour plots of posterior densities
#      opar <- par(mfrow=c(contours,contours),mar=c(0,0,1,0), omi=c(.5,.5,1,.5),
#                  mgp=c(1.5,.25,0))
#      temp.x.pos <- range(x$Z[,1,])
#      temp.y.pos <- range(x$Z[,2,])
#      require(KernSmooth,quietly=TRUE)
#      for(k in 1:nrow(x$Z.mle)){
#        plot(x=x$Z[k,1,],y=x$Z[k,2,],
#             xlab="",ylab="",type="n", asp=1,
#             xlim=temp.x.pos,ylim=temp.y.pos,
#             xaxt='n',yaxt='n')
##             main=paste(flo.names[k],sum(flo[k,])),xaxt='n',yaxt='n')
#        if((k-(contours^2)*trunc((k-1)/(contours^2)))/contours>(contours-1))
#          axis(side=1)
#        if(nrow(x$Z.mle)- k < contours)
#          axis(side=1)
#        if(contours*trunc(k/contours)+1==k){
#          axis(side=2)
#        }
#        est <- bkde2D(x=t(x$Z[k,,]), bandwidth=c(0.2,0.2))
#        est$fhat <- est$fhat/max(est$fhat,na.rm=TRUE)
#        contour(est$x1, est$x2, est$fhat, add=TRUE, nlevels=5,
#                drawlabels=FALSE,col=contour.color)
#        for (i in 1:network.size(x$network))
#          for (j in 1:network.size(x$network))
#            if (newnetwork[i,j])
#            {
#              midarrow(x0 = z.pos[i, 1], y0 = z.pos[i, 2], x1 = z.pos[j,1],
#                       y1 = z.pos[j, 2], length = al/10,lwd=edge.lwd,col=edge.col)
#            }
#        if(!is.null(label) && length(label)==1 && is.character(label)){
#          trycol <- unlist(get.vertex.attribute(x,label))
#          if(!any(is.na(trycol))){
#            label <- trycol
#          }
#        }
#        if(is.null(label) | length(label)==1){
#          label <- 1:network.size(x$network)
#        }
#        text(z.pos[, 1], z.pos[, 2], label, col = label.col, cex = label.cex)
#      }
#      par(opar)
#    }
#    else
#    {
##
##    If not densities or contours then
##    just plot the points with their MKL classes
##
#      if (!is.null(comp.mat)) {
#        procr <- function(Z, Zo) {
#          A <- t(Z) %*% Zo %*% t(Zo) %*% Z
#          A.eig <- eigen(A, symmetric = TRUE)
#          A.sqrt <- A.eig$vec %*% diag(1/sqrt(A.eig$val)) %*% 
#            t(A.eig$vec)
#          Z %*% A.sqrt %*% t(Z) %*% Zo
#        }
#        z.pos <- procr(z.pos, comp.mat)
#      }
#      plot(z.pos, xlab = xlab, ylab = ylab,
#           type = "n", asp = 1, main = main, ...)
#      aaa <- x$formula   
#      xformula <- paste(aaa[2],aaa[1],aaa[-c(1:2)],collapse=" ")
#      title(main = xformula, line = 1, cex.main = 0.7)
#      for (i in 1:network.size(x$network))
#        for (j in 1:network.size(x$network))
#          if (newnetwork[i,j])
#          {
##           midarrow(x0 = z.pos[i, 1], y0 = z.pos[i, 2], x1 = z.pos[j,1],
##                    y1 = z.pos[j, 2], length = al,lwd=edge.lwd,col=edge.col)
#            segments(x0 = z.pos[i, 1], y0 = z.pos[i, 2], x1 = z.pos[j,1],
#                     y1 = z.pos[j, 2], length = al,lwd=edge.lwd,col=edge.col)
#          }
#      if (is.null(label)) 
#        label <- 1:network.size(x$network)
#      if((!is.null(x$cluster)) & !pie)
#        {
#          points(z.pos,cex=vertex.cex,col=vertex.col,pch=vertex.pch)
#          if(label.col=="black") label.col <- "white"
#        }
#      else if(!pie)
#        {
#          points(z.pos,pch=vertex.pch,cex=vertex.cex,col=vertex.col)
##         if(label.col=="black") label.col <- "white"
#        }
#
#      if(pie)
#        for(i in 1:nrow(x$Z.mkl))
#        {
#          drawpie(z.pos[i,],piesize,x$qig[i,],n=50,cols=mycol)
#          text(z.pos[, 1], z.pos[, 2], label, col = "white", cex = label.cex)
#        }
#      else      
#        text(z.pos[, 1], z.pos[, 2], label, col = label.col, cex = label.cex)
#    }
#  }
  invisible(NULL)
}
