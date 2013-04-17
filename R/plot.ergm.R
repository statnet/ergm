#  File R/plot.ergm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
#######################################################################
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
    if(is.null(x$sample))
      cat("No plotting method available for non-latent models with MCMC sample\n")
    else
      mcmc.diagnostics(x)
  invisible(NULL)
}
