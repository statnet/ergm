#  File R/midarrow.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
#######################################################################
#=====================================================================
# This file contains the following 3 plotting functions:
#           <midarrow>
#           <drawcircle>
#           <drawpie>
#=====================================================================




##############################################################################
# The <midarrow> function adds a line, with an arrow placed in the middle of
# the line to an active plot
#
# --PARAMTERS--
#   x0,y0 : the starting x and y coordinates of the line
#   x1,y1 : the ending x and y coordinates of the line
#   length: the length of the edges of the arrow head, in inches; default=.25
#   angle : the angle from the shaft of the arrow to the edge of the arrow head;
#           default=30
#   code  : a code to indicate which way the arrow should point:
#      1 -- towards the x0,y0 point
#      2 -- towards the x1,y1 point
#      3 -- towards both points
#   col, lty, lwd, xpd have their typical par interpretations
#
##############################################################################

midarrow <- function(x0, y0, x1, y1, length = 0.25, angle = 30, code = 2, 
                     col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)
{
  ymid <- (y0+y1)/2
  xmid <- (x0+x1)/2
  arrows(x0, y0, xmid, ymid, length = length, angle = angle, code = code, 
         col = col, lty = lty, lwd = lwd, xpd = xpd)
  if(length(col)<length(x1))
    col <- rep(col[1],length(x1))
  if(length(lwd)<length(x1))
    lwd <- rep(lwd[1],length(x1))
  if(!is.null(lty))
    if(length(lty)<length(x1))
      lty <- rep(lty[1],length(x1))
  if(!is.null(xpd))
    if(length(xpd)<length(x1))
      xpd <- rep(xpd[1],length(x1))
  for(i in 1:length(x1))
    lines(c(xmid[i],x1[i]),c(ymid[i], y1[i]), col = col[i], lty = lty[i], lwd = lwd[i], xpd = xpd[i])
}



drawcircle <- function(center,radius,length=50,...)
{
  x0 <- seq(-radius,radius,length=length)
  x1 <- seq(radius,-radius,length=length)
  x <- c(x0,x1)
  y <- c(sqrt(radius^2 - x0^2),-sqrt(radius^2 - x1^2))
  lines(x+center[1],y+center[2],...)
}



drawpie <- function(center,radius,probs,n=50,cols=1:length(probs),...)
{
  x <- c(0,cumsum(probs)/sum(probs))
  dx <- diff(x)
  np <- length(probs)
  for (i in 1:np)
  {
    t2p <- 2 * pi * seq(x[i], x[i + 1], length = n)
    xc <- center[1] + c(cos(t2p), 0) * radius
    yc <- center[2] + c(sin(t2p), 0) * radius
    polygon(xc, yc, border = FALSE, col = cols[i])
  }
}
