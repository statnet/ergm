#  File tests/testthat/test-term-gw-sp.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################


niter <- 20

data(faux.dixon.high)
data(faux.mesa.high)

test.approx = function(a, b, tol=1e-6) {
  expect_lte(max(abs(a - b)), tol)
  if(!is.null(names(b))) expect_named(a, names(b))
}

test_that(paste0("DSP"), {
                                        # testing on a larger network (comparing with functions in the ergm package)
  net <- faux.dixon.high

  expect_equal(summary(net ~ dgwdsp(fixed=F, type='OTP', cutoff=5)), setNames(c(4871, 1077, 228, 53, 17), paste0("dsp.OTP#", 1:5))) # DSP OTP count
  expect_error(summary(net ~ dgwdsp(fixed=F, type='OTP', cutoff=3))) # Cutoff too small
  test.approx(summary(net ~ dgwdsp(fixed=T, type='OTP', decay=0.1)), c(gwdsp.OTP.fixed.0.1=6379.609), tol=1e-3) # GWDSP summary
  test.approx(coef(ergm(net ~ edges + dgwdsp(fixed=T, decay=0.1, type="OTP"), estimate = "MPLE"))[2], c(gwdsp.OTP.fixed.0.1=0.006467)) # MPLE estimate on larger nw

                                        # OTP and ITP are the same for dyadwise SP
  expect_equal(summary(net ~ ddsp(type="OTP", d=2)), c(dsp.OTP2=1077)) # ddsp OTP count error
  expect_equal(summary(net ~ ddsp(type="ITP", d=2)), c(dsp.ITP2=1077)) # ddsp ITP count error
  expect_equal(summary(net ~ ddsp(type="OSP", d=2)), c(dsp.OSP2=948)) # ddsp OSP count error
  expect_equal(summary(net ~ ddsp(type="ISP", d=2)), c(dsp.ISP2=1132)) # ddsp ISP count error
  expect_equal(summary(net ~ ddsp(type="RTP", d=1:2)), c(dsp.RTP1=1096,dsp.RTP2=96)) # ddsp RTP count error
})

test_that(paste0("GWDSP OTP test"), {
  net <- network.initialize(10, directed = T)
  net[1,2] <- net[2,4] <- net[1,3] <- net[3,4] <- 1
  net[5,6] <- net[6,8] <- net[8,5] <- 1
  net[10,7] <- net[10,9] <- net[9,7] <- 1
                                        #plot(net, edge.lwd=2, arrowhead.cex=2, label=1:network.size(net))

  for (i in 1:niter) {
    (espcounts <- summary(net ~ dgwdsp(fixed=F, type="OTP")) )
    expect_equal(espcounts[1:4], c(4,1,0,0), ignore_attr=TRUE) # DSP OTP mis-count
  }
                                        # espcounts

  test.approx(summary(net~dgwdsp(fixed=T, decay=0.1, type="OTP")), 5.095163) # GWDSP_OTP wrong stat

  test.approx(coef(ergm(net ~ edges + dgwdsp(fixed=T, decay=0.1, type="OTP"), estimate = "MPLE"))[2], -1.2014723) # GWDSP_OTP ergm MPLE wrong estimate
})

test_that(paste0("GWDSP ITP test"), {
  net <- network.initialize(10, directed = T)
  net[1,2] <- net[2,4] <- net[1,3] <- net[3,4] <- 1
  net[5,6] <- net[6,8] <- net[8,5] <- 1
  net[10,7] <- net[10,9] <- net[9,7] <- 1
                                        #plot(net, edge.lwd=2, arrowhead.cex=2, label=1:network.size(net))

  for (i in 1:niter) {
    (espcounts <- summary(net ~ dgwdsp(fixed=F, type="ITP")) )
    expect_equal(espcounts[1:4], c(4,1,0,0), ignore_attr=TRUE) # DSP ITP mis-count
  }
  espcounts

  test.approx(summary(net~dgwdsp(fixed=T, decay=0.1, type="ITP")), 5.095163) # GWDSP_ITP wrong stat

  test.approx(coef(ergm(net ~ edges + dgwdsp(fixed=T, decay=0.1, type="ITP"), estimate = "MPLE"))[2], -1.2014723) # GWDSP_ITP ergm MPLE wrong estimate
})

test_that(paste0("GWDSP OSP test"), {
  net <- network.initialize(10, directed = T)
  net[2,1] <- net[2,4] <- net[3,1] <- net[3,4] <- 1
  net[5,6] <- net[6,8] <- net[5,8] <- 1
  net[10,7] <-  net[9,7] <- 1
                                        #plot(net, edge.lwd=2, arrowhead.cex=2, label=1:network.size(net))

  for (i in 1:niter) {
    (espcounts <- summary(net ~ dgwdsp(fixed=F, type="OSP")) )
    expect_equal(espcounts[1:4], c(4,2,0,0), ignore_attr=TRUE) # DSP_OSP mis-count
  }
  espcounts

  test.approx(summary(net~dgwdsp(fixed=T, decay=0.1, type="OSP")), 6.190325) # GWDSP_OSP wrong stat

  test.approx(coef(ergm(net ~ edges + dgwdsp(fixed=T, decay=0.1, type="OSP"), estimate = "MPLE"))[2], -0.2104664) # GWDSP_OSP ergm MPLE wrong estimate
})

test_that(paste0("GWDSP ISP test"), {
  net <- network.initialize(10, directed = T)
  net[2,1] <- net[2,4] <- net[3,1] <- net[3,4] <- 1
  net[5,6] <- net[6,8] <- net[5,8] <- 1
  net[7,10] <-  net[7,9] <- 1
                                        #plot(net, edge.lwd=2, arrowhead.cex=2, label=1:network.size(net))

  for (i in 1:niter) {
    (espcounts <- summary(net ~ dgwdsp(fixed=F, type="ISP")) )
    expect_equal(espcounts[1:4], c(4,2,0,0), ignore_attr=TRUE) # DSP ISP mis-count
  }
  espcounts

  test.approx(summary(net~dgwdsp(fixed=T, decay=0.1, type="ISP")), 6.190325) # GWDSP_ISP wrong stat

  test.approx(coef(ergm(net ~ edges + dgwdsp(fixed=T, decay=0.1, type="ISP"), estimate = "MPLE"))[2], -0.2104664) # GWDSP_ISP ergm MPLE wrong estimate
})

                                        # ============

test_that(paste0("ESP"), {
                                        # testing on a larger network (comparing with functions in the ergm package)
  net <- faux.dixon.high
  expect_equal(summary(net ~ dgwesp(fixed=F, type='OTP', cutoff=5)), setNames(c(516, 222, 65, 21, 3), paste0("esp.OTP#", 1:5))) # ESP OTP count
  expect_error(summary(net ~ dgwesp(fixed=F, type='OTP', cutoff=4))) # Cutoff too small
  test.approx(summary(net ~ dgwesp(fixed=T, type='OTP', decay=0.1)), c(gwesp.OTP.fixed.0.1=857.4225), tol=1e-3) # GWESP summary
  test.approx(coef(ergm(net ~ edges + dgwesp(fixed=T, decay=0.1, type="OTP"), estimate = "MPLE"))[2], c(gwesp.OTP.fixed.0.1=1.807995)) # MPLE estimate on larger nw

  expect_equal(summary(net ~ desp(type="OTP", d=2)), c(esp.OTP2=222)) # desp OTP count error
  expect_equal(summary(net ~ desp(type="ITP", d=2)), c(esp.ITP2=185)) # desp ITP count error
  expect_equal(summary(net ~ desp(type="OSP", d=2)), c(esp.OSP2=215)) # desp OSP count error
  expect_equal(summary(net ~ desp(type="ISP", d=2)), c(esp.ISP2=228)) # desp ISP count error
})

test_that(paste0("GWESP OTP test"), {
  net <- network.initialize(10, directed = T)
  net[1,2] <- net[2,4] <- net[1,3] <- net[3,4] <- net[1,4] <- 1
  net[5,6] <- net[6,8] <- net[8,5] <- 1
  net[10,7] <- net[10,9] <- net[9,7] <- 1
                                        #plot(net, edge.lwd=2, arrowhead.cex=2, label=1:network.size(net))

  for (i in 1:niter) {
    (espcounts <- summary(net ~ dgwesp(fixed=F, type="OTP")) )
    expect_equal(espcounts[1:4], c(1,1,0,0), ignore_attr=TRUE) # ESP OTP mis-count
  }
  espcounts

  test.approx(summary(net~dgwesp(fixed=T, decay=0.1, type="OTP")), 2.095163) # GWESP_OTP wrong stat

  test.approx(coef(ergm(net ~ edges + dgwesp(fixed=T, decay=0.1, type="OTP"), estimate = "MPLE"))[2], 0.9157069) # GWESP_OTP ergm MPLE wrong estimate
})

test_that(paste0("GWESP ITP test"), {
  net <- network.initialize(10, directed = T)
  net[1,2] <- net[2,4] <- net[1,3] <- net[3,4] <- net[4,1] <- 1
  net[5,6] <- net[6,8] <- net[8,5] <- 1
  net[10,7] <- net[10,9] <- net[9,7] <- 1
                                        #plot(net, edge.lwd=2, arrowhead.cex=2, label=1:network.size(net))

  for (i in 1:niter) {
    (espcounts <- summary(net ~ dgwesp(fixed=F, type="ITP")) )
    expect_equal(espcounts[1:4], c(7,1,0,0), ignore_attr=TRUE) # ESP ITP mis-count
  }
  espcounts

  test.approx(summary(net~dgwesp(fixed=T, decay=0.1, type="ITP")), 8.095163) # GWESP_ITP wrong stat

  test.approx(coef(ergm(net ~ edges + dgwesp(fixed=T, decay=0.1, type="ITP"), estimate = "MPLE"))[2], 1.97197) # GWESP_ITP ergm MPLE wrong estimate
})

test_that(paste0("GWESP OSP test"), {
  net <- network.initialize(10, directed = T)
  net[2,1] <- net[2,4] <- net[3,1] <- net[3,4] <- net[3,2] <- 1
  net[5,6] <- net[6,8] <- net[5,8] <- 1
  net[10,7] <-  net[9,7] <- 1
                                        #plot(net, edge.lwd=2, arrowhead.cex=2, label=1:network.size(net))

  for (i in 1:niter) {
    (espcounts <- summary(net ~ dgwesp(fixed=F, type="OSP")) )
    expect_equal(espcounts[1:4], c(1,1,0,0), ignore_attr=TRUE) # ESP_OSP mis-count
  }
  espcounts

  test.approx(summary(net~dgwesp(fixed=T, decay=0.1, type="OSP")), 2.095163 ) # GWESP_OSP wrong stat

  test.approx(coef(ergm(net ~ edges + dgwesp(fixed=T, decay=0.1, type="OSP"), estimate = "MPLE"))[2], 1.137139) # GWESP_OSP ergm MPLE wrong estimate
})

test_that(paste0("GWESP ISP test"), {
  net <- network.initialize(10, directed = T)
  net[2,1] <- net[2,4] <- net[3,1] <- net[3,4] <- net[1,4] <- 1
  net[5,6] <- net[6,8] <- net[5,8] <- 1
  net[7,10] <-  net[7,9] <- 1
                                        #plot(net, edge.lwd=2, arrowhead.cex=2, label=1:network.size(net))

  for (i in 1:niter) {
    (espcounts <- summary(net ~ dgwesp(fixed=F, type="ISP")) )
    expect_equal(espcounts[1:4], c(1,1,0,0), ignore_attr=TRUE) # ESP ISP mis-count
  }
  espcounts

  test.approx(summary(net~dgwesp(fixed=T, decay=0.1, type="ISP")), 2.095163) # GWESP_ISP wrong stat

  test.approx(coef(ergm(net ~ edges + dgwesp(fixed=T, decay=0.1, type="ISP"), estimate = "MPLE"))[2], 1.137139) # GWESP_ISP ergm MPLE wrong estimate
})

                                        # ========

test_that(paste0("NSP"), {
                                        # testing on a larger network (comparing with functions in the ergm package)
  net <- faux.dixon.high

  expect_equal(summary(net ~ dgwnsp(fixed=F, type='OTP', cutoff=5)), setNames(c(4355, 855,163,32,14), paste0("nsp.OTP#", 1:5))) # NSP OTP count
  expect_error(summary(net ~ dgwnsp(fixed=F, type='OTP', cutoff=2))) # Cutoff too small
  test.approx(summary(net ~ dgwnsp(fixed=T, type='OTP', decay=0.1)), c(gwnsp.OTP.fixed.0.1=5522.186), tol=1e-3) # GWNSP summary
  test.approx(coef(ergm(net ~ edges + dgwnsp(fixed=T, decay=0.1, type="OTP"), estimate = "MPLE"))[2], c(gwnsp.OTP.fixed.0.1=-0.07421213)) # MPLE estimate on larger nw

  test.approx(summary(net ~ dgwdsp(fixed=T, type='OTP', decay=0.1)), unname(summary(net ~ dgwesp(fixed=T, type='OTP', decay=0.1)) + summary(net ~ dgwnsp(fixed=T, type='OTP', decay=0.1)))) # GWDSP should equal GWESP + GWNSP

  expect_equal(summary(net ~ dnsp(type="OTP", d=2)), c(nsp.OTP2=855)) # dnsp OTP count error
  expect_equal(summary(net ~ dnsp(type="ITP", d=2)), c(nsp.ITP2=892)) # dnsp ITP count error
  expect_equal(summary(net ~ dnsp(type="OSP", d=2)), c(nsp.OSP2=733)) # dnsp OSP count error
  expect_equal(summary(net ~ dnsp(type="ISP", d=2)), c(nsp.ISP2=904)) # dnsp ISP count error
})

test_that(paste0("GWNSP OTP test"), {
  net <- network.initialize(10, directed = T)
  net[1,2] <- net[2,4] <- net[1,3] <- net[3,4] <- 1
  net[5,6] <- net[6,8] <- net[8,5] <- 1
  net[10,7] <- net[10,9] <- net[9,7] <- 1
                                        #plot(net, edge.lwd=2, arrowhead.cex=2, label=1:network.size(net))

  for (i in 1:niter) {
    (espcounts <- summary(net ~ dgwnsp(fixed=F, type="OTP")) )
    expect_equal(espcounts[1:4], c(3,1,0,0), ignore_attr=TRUE) # NSP OTP mis-count
  }
  espcounts

  test.approx(summary(net~dgwnsp(fixed=T, decay=0.1, type="OTP")), 4.095163) # GWNSP_OTP wrong stat

  test.approx(coef(ergm(net ~ edges + dgwnsp(fixed=T, decay=0.1, type="OTP"), estimate = "MPLE"))[2], -0.8563381 ) # GWNSP_OTP ergm MPLE wrong estimate
})

test_that(paste0("GWNSP ITP test"), {
  net <- network.initialize(10, directed = T)
  net[1,2] <- net[2,4] <- net[1,3] <- net[3,4] <- 1
  net[5,6] <- net[6,8] <- net[8,5] <- 1
  net[10,7] <- net[10,9] <- net[9,7] <- 1
                                        #plot(net, edge.lwd=2, arrowhead.cex=2, label=1:network.size(net))

  for (i in 1:niter) {
    (espcounts <- summary(net ~ dgwnsp(fixed=F, type="ITP")) )
    expect_equal(espcounts[1:4], c(1,1,0,0), ignore_attr=TRUE) # NSP ITP mis-count
  }
  espcounts

  test.approx(summary(net~dgwnsp(fixed=T, decay=0.1, type="ITP")), 2.095163) # GWNSP_ITP wrong stat

  test.approx(coef(ergm(net ~ edges + dgwnsp(fixed=T, decay=0.1, type="ITP"), estimate = "MPLE"))[2], -1.760084 ) # GWNSP_ITP ergm MPLE wrong estimate
})

test_that(paste0("GWNSP OTP test"), {
  net <- network.initialize(10, directed = T)
  net[2,1] <- net[2,4] <- net[3,1] <- net[3,4] <- 1
  net[5,6] <- net[6,8] <- net[5,8] <- 1
  net[10,7] <-  net[9,7] <- 1
                                        #plot(net, edge.lwd=2, arrowhead.cex=2, label=1:network.size(net))

  for (i in 1:niter) {
    (espcounts <- summary(net ~ dgwnsp(fixed=F, type="OSP")) )
    expect_equal(espcounts[1:4], c(3,2,0,0), ignore_attr=TRUE) # NSP_OSP mis-count
  }
  espcounts

  test.approx(summary(net~dgwnsp(fixed=T, decay=0.1, type="OSP")), 5.190325 ) # GWNSP_OSP wrong stat

  test.approx(coef(ergm(net ~ edges + dgwnsp(fixed=T, decay=0.1, type="OSP"), estimate = "MPLE"))[2], -0.2865438) # GWNSP_OSP ergm MPLE wrong estimate
})

test_that(paste0("GWNSP ITP test"), {
  net <- network.initialize(10, directed = T)
  net[2,1] <- net[2,4] <- net[3,1] <- net[3,4] <- 1
  net[5,6] <- net[6,8] <- net[5,8] <- 1
  net[7,10] <-  net[7,9] <- 1
                                        #plot(net, edge.lwd=2, arrowhead.cex=2, label=1:network.size(net))

  for (i in 1:niter) {
    (espcounts <- summary(net ~ dgwnsp(fixed=F, type="ISP")) )
    expect_equal(espcounts[1:4], c(3,2,0,0), ignore_attr=TRUE) # NSP ISP mis-count
  }
  espcounts

  test.approx(summary(net~dgwnsp(fixed=T, decay=0.1, type="ISP")), c(gwnsp.ISP.fixed.0.1=5.190325)) # GWNSP_ISP wrong stat

  test.approx(coef(ergm(net ~ edges + dgwnsp(fixed=T, decay=0.1, type="ISP"), estimate = "MPLE")), c(edges = -1.90541331842123, gwnsp.ISP.fixed.0.1 = -0.286543786068899 )) # GWNSP_ISP ergm MPLE wrong estimate
})

test_that(paste0("DSP Undirect tests"), {
  net <- faux.mesa.high

  expect_equal(summary(net ~ dgwdsp(fixed=F, cutoff=6)), setNames(c(431, 75, 23, 1, 1, 0), paste0("dsp#", 1:6))) # DSP OTP count
  test.approx(summary(net ~ dgwdsp(fixed=T, decay=0.1)), c(gwdsp.fixed.0.1=540.7445), tol=1e-3) # GWDSP summary

  test.approx(coef(ergm(net ~ edges + dgwdsp(fixed=T, decay=0.1), estimate = "MPLE")), c(edges = -4.89690639, gwdsp.fixed.0.1 = 0.06665661)) # MPLE estimate on larger nw

  expect_equal(summary(net ~ ddsp(d=2)), c(dsp2=75)) # ddsp UTP count error
})

test_that(paste0("ESP Undirect tests"), {
  net <- faux.mesa.high

  expect_equal(summary(net ~ dgwesp(fixed=F, cutoff=6)), setNames(c(70, 36, 13, 0, 1, 0), paste0("esp#", 1:6))) # ESP UTP count
  test.approx(summary(net ~ dgwesp(fixed=T, decay=0.1)), c(gwesp.fixed.0.1=124.8859), tol=1e-3) # GWESP summary

  test.approx(coef(ergm(net ~ edges + dgwesp(fixed=T, decay=0.1), estimate = "MPLE")), c(edges=-5.295790, gwesp.fixed.0.1=1.730727)) # MPLE estimate on larger nw

  expect_equal(summary(net ~ desp(d=2)), c(esp2=36)) # desp UTP count
})

test_that(paste0("NSP Undirect tests"), {
  net <- faux.mesa.high
  expect_equal(summary(net ~ dgwnsp(fixed=F, cutoff=5)), setNames(c(361, 39, 10, 1, 0), paste0("nsp#", 1:5))) # NSP UTP count
  test.approx(summary(net ~ dgwnsp(fixed=T, decay=0.1)), c(gwnsp.fixed.0.1=415.8586), tol=1e-3) # GWNSP summary

  test.approx(coef(ergm(net ~ edges + dgwnsp(fixed=T, decay=0.1), estimate = "MPLE")), c(edges=-4.2641998, gwnsp.fixed.0.1=-0.1072723)) # MPLE estimate on larger nw

  expect_equal(summary(net ~ dnsp(d=2)), c(nsp2=39)) # dnsp UTP count error
})
