#  File R/midarrow.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
#=====================================================================
# This file contains the following 3 plotting functions:
#           <drawpie>
#=====================================================================







drawpie <- function(center,radius,probs,n=50,cols=1:length(probs),...)
{
  .Deprecated("'ergmm.drawpie' in the latentnet package")
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
