#  File tests/testthat/test-distance.R in package ergm, part of the Statnet
#  suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################
library(ergm)

#Create some coordinates
co<-rbind(
 c(1,2,3),
 c(0,0,0),
 c(3,2,1)
)

#Compute L1 and L2 distances
d1<-as.matrix(dist(co, method="manhattan"))
d2<-as.matrix(dist(co, method="euclidean"))

#Compute stats on a test graph
g<-rbind(
 c(0,1,1),
 c(1,0,0),
 c(0,0,0)
)
net<-network(g)
net%n%"co"<-co         #Stuff coords in a network object
net%v%"x"<-co[,1]      #Also list them as individual vertex coords
net%v%"y"<-co[,2]
net%v%"z"<-co[,3]

#Verify that values are correct
test_that("L2 works", {
  expect_equal(summary(net~distance(co,log=FALSE)), sum(d2[g>0]), ignore_attr=TRUE)
})
test_that("L1 works", {
  expect_equal(summary(net~distance(co,metric=1,log=FALSE)), sum(d1[g>0]), ignore_attr=TRUE)
})
test_that("L2 works", {
  expect_equal(summary(net~distance(co)), sum(log(d2[g>0])), ignore_attr=TRUE)
})
test_that("log L1 works", {
  expect_equal(summary(net~distance(co,metric=1)), sum(log(d1[g>0])), ignore_attr=TRUE)
})
test_that("scaling works", {
  expect_equal(summary(net~distance(co,metric=1,log=FALSE,scale=3)), 3*sum(d1[g>0]), ignore_attr=TRUE)
})
test_that("powers work", {
  expect_equal(summary(net~distance(co,metric=1,log=FALSE,pow=3)), sum(d1[g>0]^3), ignore_attr=TRUE)
})
test_that("thresholding works", {
  expect_equal(summary(net~distance(co,metric=1,log=TRUE,scale=1e-10,mindist=5)), log(5)*sum(g>0), ignore_attr=TRUE)
})

#Verify that alternative specifications work
test_that("character matches matrix", {
  expect_equal(summary(net~distance(co)), summary(net~distance("co")), ignore_attr=TRUE)
})
test_that("character matches data.frame", {
  expect_equal(summary(net~distance(as.data.frame(co))), summary(net~distance("co")), ignore_attr=TRUE)
})
test_that("vertex attributes matches matrix", {
  expect_equal(summary(net~distance(co)), summary(net~distance(c("x","y","z"))), ignore_attr=TRUE)
})
test_that("single attribute defaults to vertex", {
  expect_equal(summary(net~distance(co[,1])), summary(net~distance("x")), ignore_attr=TRUE)
})


#Test for lat/lon calculations
co<-rbind(
  AustinTX=c(30.26694,-97.74278),
  BaltimoreMD=c(39.29028,-76.6125),
  BatonRougeLA=c(30.45056,-91.15444),
  CharlotteNC=c(35.22694,-80.84333),
  ChicagoIL=c(41.85,-87.65),
  ColumbusOH=c(39.96111,-82.99889),
  DetroitMI=c(42.33139,-83.04583),
  DurhamNC=c(35.99389,-78.89889),
  HoustonTX=c(29.76306,-95.36306),
  IndianapolisIN=c(39.76833,-86.15806),
  IrvineCA=c(33.66944,-117.8222),
  IthacaNY=c(42.44056,-76.49694),
  JacksonvilleFL=c(30.33194,-81.65583),
  KonaHI=c(19.64083,-155.9856),
  LasVegasNV=c(36.175,-115.1364),
  LosAngelesCA=c(34.05222,-118.2428),
  MemphisTN=c(35.14944,-90.04889),
  MiamiFL=c(25.77389,-80.19389),
  MilwaukeeWI=c(43.03889,-87.90639),
  MobileAL=c(30.69417,-88.04306),
  NewOrleansLA=c(29.95444,-90.075),
  NewYorkNY=c(40.71417,-74.00639),
  PhiladelphiaPA=c(39.95222,-75.16417),
  PhoenixAZ=c(33.44833,-112.0733),
  PittsburghPA=c(40.39556,-79.83889),
  SanAntonioTX=c(29.42389,-98.49333),
  SanDiegoCA=c(32.71528,-117.1564),
  SanFranciscoCA=c(37.775,-122.4183),
  SanJoseCA=c(37.33944,-121.8939),
  SeattleWA=c(47.60639,-122.3308)
)

#By default, spherical distances give us great circle distances
#on the geosphere, in kilometers
net <- network.initialize(30, directed = FALSE)
net[1,8] <- 1  #Create a tie from Austin, TX to Durham, NC
test_that("first great circle distance is correct", {
  expect_equal(round(summary(net ~ distance(co, sphere = TRUE, log = FALSE)),3), 1862.879, ignore_attr=TRUE)
})
net[1,8] <- 0
net[14,22] <- 1  #Now try Kona, HI to NYC, NY
test_that("second great circle distance is correct", {
  expect_equal(round(summary(net ~ distance(co, sphere = TRUE, log = FALSE)),3), 7940.019, ignore_attr=TRUE)
})
net[1,8] <- 1
test_that("great circle distance sum is correct", {
  expect_equal(round(summary(net ~ distance(co, sphere = TRUE, log = FALSE)),3), 9802.898, ignore_attr=TRUE)
})
test_that("log great circle distance sum is correct", {
  expect_equal(round(summary(net ~ distance(co, sphere = TRUE, log = TRUE)),3), 16.51, ignore_attr=TRUE)
})

