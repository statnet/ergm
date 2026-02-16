#  File tests/testthat/test-checkpointing.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################

data(florentine)
set.seed(0)

niterations <- NA
tmpf <- tempfile()
checkpoints <- paste0(tmpf,"_00",1:2,".RData")

withr::local_file(checkpoints, {
  test_that("checkpointing works",{
    expect_error(gest <- ergm(flomarriage ~ kstar(1:2) + absdiff("wealth") + triangle, control=control.ergm(MCMLE.termination="Hummel", checkpoint=paste0(tmpf,"_%03d.RData"))), NA)
    niterations <<- gest$iterations # Save for later.
  })

  test_that("resuming works",{
    gest <- ergm(flomarriage ~ kstar(1:2) + absdiff("wealth") + triangle, control=control.ergm(MCMLE.termination="Hummel", resume=sprintf(paste0(tmpf,"_%03d.RData"),niterations)))
    expect_equal(gest$iterations, 1) # It'll take at least 2 if resume fails.
  })
})
