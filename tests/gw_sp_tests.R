#  File tests/gw_sp_tests.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################

niter <- 20

library(ergm)

data("faux.dixon.high")

test.approx = function(a, b, message, tol=1e-5) {
  if (abs(a-b) > tol) stop(message)
}

run.tests <- function(cache.sp){
  message("Testing with shared partner cache ",if(cache.sp)"enabled"else"disabled",".")
  
  ctrl <- control.ergm(term.options=list(cache.sp=cache.sp))

# testing on a larger network (comparing with functions in the ergm package)
net <- faux.dixon.high
if (!all(summary(cache.sp=cache.sp,net ~ dgwdsp(fixed=F, type='OTP'))[1:5] == c(4871, 1077, 228, 53, 17)))
  stop("DSP OTP count on large network incorrect")
test.approx(summary(cache.sp=cache.sp,net ~ dgwdsp(fixed=T, type='OTP', decay=0.1)), 6379.609, tol=1e-3, "GWDSP summary")

test.approx(coef(ergm(control=ctrl,net ~ edges + dgwdsp(fixed=T, decay=0.1, type="OTP"), estimate = "MPLE"))[2], 0.006467, "MPLE estimate on larger nw")

# OTP and ITP are the same for dyadwise SP
if (summary(cache.sp=cache.sp,net ~ ddsp(type="OTP", d=2)) != 1077) stop("ddsp OTP count error")
if (summary(cache.sp=cache.sp,net ~ ddsp(type="ITP", d=2)) != 1077) stop("ddsp ITP count error")
if (summary(cache.sp=cache.sp,net ~ ddsp(type="OSP", d=2)) != 948) stop("ddsp OSP count error")
if (summary(cache.sp=cache.sp,net ~ ddsp(type="ISP", d=2)) != 1132) stop("ddsp ISP count error")
if (any(summary(cache.sp=cache.sp,net ~ ddsp(type="RTP", d=1:2)) != c(1096,96))) stop("ddsp RTP count error")

# GWDSP OTP test

net <- network.initialize(10, directed = T)
net[1,2] <- net[2,4] <- net[1,3] <- net[3,4] <- 1
net[5,6] <- net[6,8] <- net[8,5] <- 1
net[10,7] <- net[10,9] <- net[9,7] <- 1
#plot(net, edge.lwd=2, arrowhead.cex=2, label=1:network.size(net))

for (i in 1:niter) {
  (espcounts <- summary(cache.sp=cache.sp,net ~ dgwdsp(fixed=F, type="OTP")) )
  if (!all(espcounts[1:4]==c(4,1,0,0))) stop("DSP OTP mis-count")
}
espcounts

test.approx(summary(cache.sp=cache.sp,net~dgwdsp(fixed=T, decay=0.1, type="OTP")), 5.095163, "GWDSP_OTP wrong stat")

test.approx(coef(ergm(control=ctrl,net ~ edges + dgwdsp(fixed=T, decay=0.1, type="OTP"), estimate = "MPLE"))[2], -1.2014723, "GWDSP_OTP ergm MPLE wrong estimate")

# larger network


# GWDSP ITP test

net <- network.initialize(10, directed = T)

net[1,2] <- net[2,4] <- net[1,3] <- net[3,4] <- 1
net[5,6] <- net[6,8] <- net[8,5] <- 1
net[10,7] <- net[10,9] <- net[9,7] <- 1

#plot(net, edge.lwd=2, arrowhead.cex=2, label=1:network.size(net))

for (i in 1:niter) {
  (espcounts <- summary(cache.sp=cache.sp,net ~ dgwdsp(fixed=F, type="ITP")) )
  if (!all(espcounts[1:4]==c(4,1,0,0))) stop("DSP ITP mis-count")
}
espcounts

test.approx(summary(cache.sp=cache.sp,net~dgwdsp(fixed=T, decay=0.1, type="ITP")), 5.095163, "GWDSP_ITP wrong stat")

test.approx(coef(ergm(control=ctrl,net ~ edges + dgwdsp(fixed=T, decay=0.1, type="ITP"), estimate = "MPLE"))[2], -1.2014723, "GWDSP_ITP ergm MPLE wrong estimate")



# GWDSP OSP test

net <- network.initialize(10, directed = T)

net[2,1] <- net[2,4] <- net[3,1] <- net[3,4] <- 1
net[5,6] <- net[6,8] <- net[5,8] <- 1
net[10,7] <-  net[9,7] <- 1

#plot(net, edge.lwd=2, arrowhead.cex=2, label=1:network.size(net))

for (i in 1:niter) {
  (espcounts <- summary(cache.sp=cache.sp,net ~ dgwdsp(fixed=F, type="OSP")) )
  if (!all(espcounts[1:4]==c(4,2,0,0))) stop("DSP_OSP mis-count")
}
espcounts

test.approx(summary(cache.sp=cache.sp,net~dgwdsp(fixed=T, decay=0.1, type="OSP")), 6.190325, "GWDSP_OSP wrong stat")

test.approx(coef(ergm(control=ctrl,net ~ edges + dgwdsp(fixed=T, decay=0.1, type="OSP"), estimate = "MPLE"))[2], -0.2104664, "GWDSP_OSP ergm MPLE wrong estimate")


# GWDSP ISP test

net <- network.initialize(10, directed = T)

net[2,1] <- net[2,4] <- net[3,1] <- net[3,4] <- 1
net[5,6] <- net[6,8] <- net[5,8] <- 1
net[7,10] <-  net[7,9] <- 1

#plot(net, edge.lwd=2, arrowhead.cex=2, label=1:network.size(net))

for (i in 1:niter) {
  (espcounts <- summary(cache.sp=cache.sp,net ~ dgwdsp(fixed=F, type="ISP")) )
  if (!all(espcounts[1:4]==c(4,2,0,0))) stop("DSP ISP mis-count")
}
espcounts

test.approx(summary(cache.sp=cache.sp,net~dgwdsp(fixed=T, decay=0.1, type="ISP")), 6.190325, "GWDSP_ISP wrong stat")

test.approx(coef(ergm(control=ctrl,net ~ edges + dgwdsp(fixed=T, decay=0.1, type="ISP"), estimate = "MPLE"))[2], -0.2104664, "GWDSP_ISP ergm MPLE wrong estimate")


# ============


library(ergm)

data("faux.dixon.high")

test.approx = function(a, b, message, tol=1e-6) {
  if (abs(a-b) > tol) stop(message)
}

# testing on a larger network (comparing with functions in the ergm package)
net <- faux.dixon.high
if (!all(summary(cache.sp=cache.sp,net ~ dgwesp(fixed=F, type='OTP'))[1:5] == c(516, 222, 65, 21, 3)))
  stop("ESP OTP count on large network incorrect")
test.approx(summary(cache.sp=cache.sp,net ~ dgwesp(fixed=T, type='OTP', decay=0.1)), 857.4225, tol=1e-3, "GWESP summary")
test.approx(coef(ergm(control=ctrl,net ~ edges + dgwesp(fixed=T, decay=0.1, type="OTP"), estimate = "MPLE"))[2], 1.807995, "MPLE estimate on larger nw")

if (summary(cache.sp=cache.sp,net ~ desp(type="OTP", d=2)) != 222) stop("desp OTP count error")
if (summary(cache.sp=cache.sp,net ~ desp(type="ITP", d=2)) != 185) stop("desp ITP count error")
if (summary(cache.sp=cache.sp,net ~ desp(type="OSP", d=2)) != 215) stop("desp OSP count error")
if (summary(cache.sp=cache.sp,net ~ desp(type="ISP", d=2)) != 228) stop("desp ISP count error")


# GWESP OTP test

net <- network.initialize(10, directed = T)
net[1,2] <- net[2,4] <- net[1,3] <- net[3,4] <- net[1,4] <- 1
net[5,6] <- net[6,8] <- net[8,5] <- 1
net[10,7] <- net[10,9] <- net[9,7] <- 1
#plot(net, edge.lwd=2, arrowhead.cex=2, label=1:network.size(net))

for (i in 1:niter) {
  (espcounts <- summary(cache.sp=cache.sp,net ~ dgwesp(fixed=F, type="OTP")) )
  if (!all(espcounts[1:4]==c(1,1,0,0))) stop("ESP OTP mis-count")
}
espcounts

test.approx(summary(cache.sp=cache.sp,net~dgwesp(fixed=T, decay=0.1, type="OTP")), 2.095163, "GWESP_OTP wrong stat")

test.approx(coef(ergm(control=ctrl,net ~ edges + dgwesp(fixed=T, decay=0.1, type="OTP"), estimate = "MPLE"))[2], 0.9157069, "GWESP_OTP ergm MPLE wrong estimate")


# GWESP ITP test

net <- network.initialize(10, directed = T)
net[1,2] <- net[2,4] <- net[1,3] <- net[3,4] <- net[4,1] <- 1
net[5,6] <- net[6,8] <- net[8,5] <- 1
net[10,7] <- net[10,9] <- net[9,7] <- 1

#plot(net, edge.lwd=2, arrowhead.cex=2, label=1:network.size(net))

for (i in 1:niter) {
  (espcounts <- summary(cache.sp=cache.sp,net ~ dgwesp(fixed=F, type="ITP")) )
  if (!all(espcounts[1:4]==c(7,1,0,0))) stop("ESP ITP mis-count")
}
espcounts

test.approx(summary(cache.sp=cache.sp,net~dgwesp(fixed=T, decay=0.1, type="ITP")), 8.095163, "GWESP_ITP wrong stat")

test.approx(coef(ergm(control=ctrl,net ~ edges + dgwesp(fixed=T, decay=0.1, type="ITP"), estimate = "MPLE"))[2], 1.97197, "GWESP_ITP ergm MPLE wrong estimate")



# GWESP OSP test

net <- network.initialize(10, directed = T)

net[2,1] <- net[2,4] <- net[3,1] <- net[3,4] <- net[3,2] <- 1
net[5,6] <- net[6,8] <- net[5,8] <- 1
net[10,7] <-  net[9,7] <- 1

#plot(net, edge.lwd=2, arrowhead.cex=2, label=1:network.size(net))

for (i in 1:niter) {
  (espcounts <- summary(cache.sp=cache.sp,net ~ dgwesp(fixed=F, type="OSP")) )
  if (!all(espcounts[1:4]==c(1,1,0,0))) stop("ESP_OSP mis-count")
}
espcounts

test.approx(summary(cache.sp=cache.sp,net~dgwesp(fixed=T, decay=0.1, type="OSP")), 2.095163 , "GWESP_OSP wrong stat")

test.approx(coef(ergm(control=ctrl,net ~ edges + dgwesp(fixed=T, decay=0.1, type="OSP"), estimate = "MPLE"))[2], 1.137139, "GWESP_OSP ergm MPLE wrong estimate")


# GWESP ISP test

net <- network.initialize(10, directed = T)

net[2,1] <- net[2,4] <- net[3,1] <- net[3,4] <- net[1,4] <- 1
net[5,6] <- net[6,8] <- net[5,8] <- 1
net[7,10] <-  net[7,9] <- 1

#plot(net, edge.lwd=2, arrowhead.cex=2, label=1:network.size(net))

for (i in 1:niter) {
  (espcounts <- summary(cache.sp=cache.sp,net ~ dgwesp(fixed=F, type="ISP")) )
  if (!all(espcounts[1:4]==c(1,1,0,0))) stop("ESP ISP mis-count")
}
espcounts

test.approx(summary(cache.sp=cache.sp,net~dgwesp(fixed=T, decay=0.1, type="ISP")), 2.095163, "GWESP_ISP wrong stat")

test.approx(coef(ergm(control=ctrl,net ~ edges + dgwesp(fixed=T, decay=0.1, type="ISP"), estimate = "MPLE"))[2], 1.137139, "GWESP_ISP ergm MPLE wrong estimate")

# ========


library(ergm)

data("faux.dixon.high")

test.approx = function(a, b, message, tol=1e-5) {
  if (abs(a-b) > tol) stop(message)
}

# testing on a larger network (comparing with functions in the ergm package)
net <- faux.dixon.high
if (!all(summary(cache.sp=cache.sp,net ~ dgwnsp(fixed=F, type='OTP'))[1:5] == c(4355, 855,163,32,14)))
  stop("NSP OTP count on large network incorrect")
test.approx(summary(cache.sp=cache.sp,net ~ dgwnsp(fixed=T, type='OTP', decay=0.1)), 5522.186 , "GWNSP summary", tol=1e-3)
test.approx(coef(ergm(control=ctrl,net ~ edges + dgwnsp(fixed=T, decay=0.1, type="OTP"), estimate = "MPLE"))[2], -0.07421213 , "MPLE estimate on larger nw")

test.approx(summary(cache.sp=cache.sp,net ~ dgwdsp(fixed=T, type='OTP', decay=0.1)), summary(cache.sp=cache.sp,net ~ dgwesp(fixed=T, type='OTP', decay=0.1)) + summary(cache.sp=cache.sp,net ~ dgwnsp(fixed=T, type='OTP', decay=0.1)),
            "GWDSP should equal GWESP + GWNSP")

if (summary(cache.sp=cache.sp,net ~ dnsp(type="OTP", d=2)) != 855) stop("dnsp OTP count error")
if (summary(cache.sp=cache.sp,net ~ dnsp(type="ITP", d=2)) != 892) stop("dnsp ITP count error")
if (summary(cache.sp=cache.sp,net ~ dnsp(type="OSP", d=2)) != 733) stop("dnsp OSP count error")
if (summary(cache.sp=cache.sp,net ~ dnsp(type="ISP", d=2)) != 904) stop("dnsp ISP count error")

# GWNSP OTP test

net <- network.initialize(10, directed = T)
net[1,2] <- net[2,4] <- net[1,3] <- net[3,4] <- 1
net[5,6] <- net[6,8] <- net[8,5] <- 1
net[10,7] <- net[10,9] <- net[9,7] <- 1
#plot(net, edge.lwd=2, arrowhead.cex=2, label=1:network.size(net))

for (i in 1:niter) {
  (espcounts <- summary(cache.sp=cache.sp,net ~ dgwnsp(fixed=F, type="OTP")) )
  if (!all(espcounts[1:4]==c(3,1,0,0))) stop("NSP OTP mis-count")
}
espcounts

test.approx(summary(cache.sp=cache.sp,net~dgwnsp(fixed=T, decay=0.1, type="OTP")), 4.095163, "GWNSP_OTP wrong stat")

test.approx(coef(ergm(control=ctrl,net ~ edges + dgwnsp(fixed=T, decay=0.1, type="OTP"), estimate = "MPLE"))[2], -0.8563381 , "GWNSP_OTP ergm MPLE wrong estimate")



# GWNSP ITP test

net <- network.initialize(10, directed = T)

net[1,2] <- net[2,4] <- net[1,3] <- net[3,4] <- 1
net[5,6] <- net[6,8] <- net[8,5] <- 1
net[10,7] <- net[10,9] <- net[9,7] <- 1

#plot(net, edge.lwd=2, arrowhead.cex=2, label=1:network.size(net))

for (i in 1:niter) {
  (espcounts <- summary(cache.sp=cache.sp,net ~ dgwnsp(fixed=F, type="ITP")) )
  if (!all(espcounts[1:4]==c(1,1,0,0))) stop("NSP ITP mis-count")
}
espcounts

test.approx(summary(cache.sp=cache.sp,net~dgwnsp(fixed=T, decay=0.1, type="ITP")), 2.095163, "GWNSP_ITP wrong stat")

test.approx(coef(ergm(control=ctrl,net ~ edges + dgwnsp(fixed=T, decay=0.1, type="ITP"), estimate = "MPLE"))[2], -1.760084 , "GWNSP_ITP ergm MPLE wrong estimate")



# GWNSP OSP test

net <- network.initialize(10, directed = T)

net[2,1] <- net[2,4] <- net[3,1] <- net[3,4] <- 1
net[5,6] <- net[6,8] <- net[5,8] <- 1
net[10,7] <-  net[9,7] <- 1

#plot(net, edge.lwd=2, arrowhead.cex=2, label=1:network.size(net))

for (i in 1:niter) {
  (espcounts <- summary(cache.sp=cache.sp,net ~ dgwnsp(fixed=F, type="OSP")) )
  if (!all(espcounts[1:4]==c(3,2,0,0))) stop("NSP_OSP mis-count")
}
espcounts

test.approx(summary(cache.sp=cache.sp,net~dgwnsp(fixed=T, decay=0.1, type="OSP")), 5.190325 , "GWNSP_OSP wrong stat")

test.approx(coef(ergm(control=ctrl,net ~ edges + dgwnsp(fixed=T, decay=0.1, type="OSP"), estimate = "MPLE"))[2], -0.2865438, "GWNSP_OSP ergm MPLE wrong estimate")


# GWNSP ISP test

net <- network.initialize(10, directed = T)

net[2,1] <- net[2,4] <- net[3,1] <- net[3,4] <- 1
net[5,6] <- net[6,8] <- net[5,8] <- 1
net[7,10] <-  net[7,9] <- 1

#plot(net, edge.lwd=2, arrowhead.cex=2, label=1:network.size(net))

for (i in 1:niter) {
  (espcounts <- summary(cache.sp=cache.sp,net ~ dgwnsp(fixed=F, type="ISP")) )
  if (!all(espcounts[1:4]==c(3,2,0,0))) stop("NSP ISP mis-count")
}
espcounts

test.approx(summary(cache.sp=cache.sp,net~dgwnsp(fixed=T, decay=0.1, type="ISP")), 5.190325, "GWNSP_ISP wrong stat")

test.approx(coef(ergm(control=ctrl,net ~ edges + dgwnsp(fixed=T, decay=0.1, type="ISP"), estimate = "MPLE"))[2], -0.2865438, "GWNSP_ISP ergm MPLE wrong estimate")


# undirected tests
data(faux.mesa.high)
net <- faux.mesa.high

# dsp
if (!all(summary(cache.sp=cache.sp,net ~ dgwdsp(fixed=F)) == summary(cache.sp=cache.sp,net ~ gwdsp(fixed=F))))
  stop("DSP OTP count on large network incorrect")
test.approx(summary(cache.sp=cache.sp,net ~ dgwdsp(fixed=T, decay=0.1)), summary(cache.sp=cache.sp,net ~ gwdsp(fixed=T, decay=0.1)), tol=1e-3, "GWDSP summary")

test.approx(coef(ergm(control=ctrl,net ~ edges + dgwdsp(fixed=T, decay=0.1), estimate = "MPLE"))[2], coef(ergm(control=ctrl,net ~ edges + gwdsp(fixed=T, decay=0.1), estimate = "MPLE"))[2], "MPLE estimate on larger nw")

if (summary(cache.sp=cache.sp,net ~ ddsp(d=2)) != summary(cache.sp=cache.sp,net ~ dsp(d=2))) stop("ddsp UTP count error")

# esp
if (!all(summary(cache.sp=cache.sp,net ~ dgwesp(fixed=F)) == summary(cache.sp=cache.sp,net ~ gwesp(fixed=F))))
  stop("ESP UTP count on large network incorrect")
test.approx(summary(cache.sp=cache.sp,net ~ dgwesp(fixed=T, decay=0.1)), summary(cache.sp=cache.sp,net ~ gwesp(fixed=T, decay=0.1)), tol=1e-3, "GWESP summary")

test.approx(coef(ergm(control=ctrl,net ~ edges + dgwesp(fixed=T, decay=0.1), estimate = "MPLE"))[2], coef(ergm(control=ctrl,net ~ edges + gwesp(fixed=T, decay=0.1), estimate = "MPLE"))[2], "MPLE estimate on larger nw")

if (summary(cache.sp=cache.sp,net ~ desp(d=2)) != summary(cache.sp=cache.sp,net ~ esp(d=2))) stop("desp UTP count error")


# nsp
if (!all(summary(cache.sp=cache.sp,net ~ dgwnsp(fixed=F)) == summary(cache.sp=cache.sp,net ~ gwnsp(fixed=F))))
  stop("NSP UTP count on large network incorrect")
test.approx(summary(cache.sp=cache.sp,net ~ dgwnsp(fixed=T, decay=0.1)), summary(cache.sp=cache.sp,net ~ gwnsp(fixed=T, decay=0.1)), tol=1e-3, "GWNSP summary")

test.approx(coef(ergm(control=ctrl,net ~ edges + dgwnsp(fixed=T, decay=0.1), estimate = "MPLE"))[2], coef(ergm(control=ctrl,net ~ edges + gwnsp(fixed=T, decay=0.1), estimate = "MPLE"))[2], "MPLE estimate on larger nw")

if (summary(cache.sp=cache.sp,net ~ dnsp(d=2)) != summary(cache.sp=cache.sp,net ~ nsp(d=2))) stop("dnsp UTP count error")
}

run.tests(FALSE)
run.tests(TRUE)
