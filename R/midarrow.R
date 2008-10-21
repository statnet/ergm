#  File ergm/R/midarrow.R
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
