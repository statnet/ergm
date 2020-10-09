#  File tests/testthat/test-checkpointing.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################

data(florentine)
set.seed(0)
tmpf <- tempfile()

teardown({
  unlink(paste0(tmpf,"_00",1:2,".RData"))
})

test_that("checkpointing works",{
 gest <- ergm(flomarriage ~ kstar(1:2) + absdiff("wealth") + triangle, control=control.ergm(checkpoint=paste0(tmpf,"_%03d.RData")))
})

test_that("resuming works",{
  gest <- ergm(flomarriage ~ kstar(1:2) + absdiff("wealth") + triangle, control=control.ergm(resume=paste0(tmpf,"_002.RData")))
  expect_equal(gest$iterations, 1) # It'll take at least 2 if resume fails.
})
