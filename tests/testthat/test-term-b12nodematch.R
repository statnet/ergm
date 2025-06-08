#  File tests/testthat/test-term-b12nodematch.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################


## bipartite network
 el <- cbind( c(17, 20, 22, 26, 19, 24, 16, 22, 18, 23, 28, 20,
               22, 23, 17, 21, 25, 21, 27, 16, 19, 18, 23, 26, 16),
            c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 9, 10,
              10, 11, 11, 4, 2))
 mynw <- network(el, bipartite=15, directed=FALSE)
 mynw %v% "names" <- rep(letters[1:3], c(10,10,8))

# Plot is possible if desired; normally commented out
# plot(mynw, label=mynw %v% "names", vertex.col=rep(2:3,c(15,13)))
# plot(mynw, label=paste(mynw %v% "names",1:28,sep=""), vertex.col=rep(2:3,c(15,13)))

plot(mynw, label=paste(mynw %v% "names",1:28,sep=""), vertex.col=rep(2:3,c(15,13)))

test_that("term b12 node match", {
  expect_equal(summary(mynw~b1nodematch(~names, beta=1)), 12, ignore_attr=TRUE)
  expect_equal(summary(mynw~b1nodematch(function(x) x %v% "names", beta=0)), 9, ignore_attr=TRUE)

  expect_equal(summary(mynw~b2nodematch("names", beta=1)), 7, ignore_attr=TRUE)
  expect_equal(summary(mynw~b2nodematch(function(x) x %v% "names", beta=0)), 7, ignore_attr=TRUE)

  expect_equal(summary(mynw~b1nodematch(~names, diff=TRUE, levels=1, beta=1)), 12, ignore_attr=TRUE)
  expect_equal(summary(mynw~b1nodematch(function(x) x %v% "names", diff=TRUE, keep=1, beta=0)), 9, ignore_attr=TRUE)

  expect_equal(summary(mynw~b2nodematch("names", diff=TRUE, levels=2, beta=1)), 5, ignore_attr=TRUE)
  expect_equal(summary(mynw~b2nodematch(function(x) x %v% "names", diff=TRUE, keep=2, beta=0)), 5, ignore_attr=TRUE)

  expect_equal(summary(mynw~b1nodematch("names", diff=TRUE, alpha=0)), c(10,0), ignore_attr=TRUE)

  expect_equal(summary(mynw~b1nodematch(function(x) x %v% "names", diff=FALSE, alpha=0)), 10, ignore_attr=TRUE)

  expect_equal(summary(mynw~b2nodematch(~names, diff=TRUE, alpha=0)), c(2,4), ignore_attr=TRUE)

  expect_equal(summary(mynw~b2nodematch(function(x) x %v% "names", diff=FALSE, alpha=0)), 6, ignore_attr=TRUE)
})
