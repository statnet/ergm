#  File tests/testthat/test-bipartite-missing-data.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################

test_that("bipartite terms abort for missing data when they should and not when they shouldn't", {
  nw <- network.initialize(10,bip=4,dir=FALSE)
  nw %v% "b1only" <- c(1,2,1,2,NA,NA,NA,NA,NA,NA)
  nw %v% "b2only" <- c(NA,NA,NA,NA,3,3,3,4,4,4)
  nw %v% "both" <- c(1,2,1,2,3,3,3,4,4,4)
  
  expect_error(summary(nw ~ b1degrange(from=1, by="b1only")), NA)
  expect_error(summary(nw ~ b1degrange(from=1, by="b2only")))
  expect_error(summary(nw ~ b1degrange(from=1, by="both")), NA)
               
  expect_error(summary(nw ~ b1degrange(from=1, by="b1only", homophily=TRUE)))
  expect_error(summary(nw ~ b1degrange(from=1, by="b2only", homophily=TRUE)))
  expect_error(summary(nw ~ b1degrange(from=1, by="both", homophily=TRUE)), NA)

  expect_error(summary(nw ~ b2degrange(from=1, by="b1only")))
  expect_error(summary(nw ~ b2degrange(from=1, by="b2only")), NA)
  expect_error(summary(nw ~ b2degrange(from=1, by="both")), NA)
                             
  expect_error(summary(nw ~ b2degrange(from=1, by="b1only", homophily=TRUE)))
  expect_error(summary(nw ~ b2degrange(from=1, by="b2only", homophily=TRUE)))
  expect_error(summary(nw ~ b2degrange(from=1, by="both", homophily=TRUE)), NA)

  expect_error(summary(nw ~ gwb1degree(decay=0, fixed=TRUE, attr="b1only")), NA)
  expect_error(summary(nw ~ gwb1degree(decay=0, fixed=TRUE, attr="b2only")))
  expect_error(summary(nw ~ gwb1degree(decay=0, fixed=TRUE, attr="both")), NA)

  expect_error(summary(nw ~ gwb2degree(decay=0, fixed=TRUE, attr="b1only")))
  expect_error(summary(nw ~ gwb2degree(decay=0, fixed=TRUE, attr="b2only")), NA)
  expect_error(summary(nw ~ gwb2degree(decay=0, fixed=TRUE, attr="both")), NA)
  
  expect_error(summary(nw ~ b1nodematch(attr="b1only")), NA)
  expect_error(summary(nw ~ b1nodematch(attr="b2only")))
  expect_error(summary(nw ~ b1nodematch(attr="both")), NA)
  expect_error(summary(nw ~ b1nodematch(attr="b1only", byb2attr="b1only")))
  expect_error(summary(nw ~ b1nodematch(attr="b1only", byb2attr="b2only")), NA)
  expect_error(summary(nw ~ b1nodematch(attr="b1only", byb2attr="both")), NA)
  expect_error(summary(nw ~ b1nodematch(attr="b2only", byb2attr="b1only")))
  expect_error(summary(nw ~ b1nodematch(attr="b2only", byb2attr="b2only")))
  expect_error(summary(nw ~ b1nodematch(attr="b2only", byb2attr="both")))
  expect_error(summary(nw ~ b1nodematch(attr="both", byb2attr="b1only")))
  expect_error(summary(nw ~ b1nodematch(attr="both", byb2attr="b2only")), NA)
  expect_error(summary(nw ~ b1nodematch(attr="both", byb2attr="both")), NA)

  expect_error(summary(nw ~ b2nodematch(attr="b1only")))
  expect_error(summary(nw ~ b2nodematch(attr="b2only")), NA)
  expect_error(summary(nw ~ b2nodematch(attr="both")), NA)
  expect_error(summary(nw ~ b2nodematch(attr="b1only", byb1attr="b1only")))
  expect_error(summary(nw ~ b2nodematch(attr="b1only", byb1attr="b2only")))
  expect_error(summary(nw ~ b2nodematch(attr="b1only", byb1attr="both")))
  expect_error(summary(nw ~ b2nodematch(attr="b2only", byb1attr="b1only")), NA)
  expect_error(summary(nw ~ b2nodematch(attr="b2only", byb1attr="b2only")))
  expect_error(summary(nw ~ b2nodematch(attr="b2only", byb1attr="both")), NA)
  expect_error(summary(nw ~ b2nodematch(attr="both", byb1attr="b1only")), NA)
  expect_error(summary(nw ~ b2nodematch(attr="both", byb1attr="b2only")))
  expect_error(summary(nw ~ b2nodematch(attr="both", byb1attr="both")), NA)
})
