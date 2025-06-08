#  File tests/testthat/test-ergm-proposal-unload.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
test_that("proposal table is purged when package is unloaded", {
  library(ergm.count)
  expect_true("DiscTNT" %in% ergm_proposal_table()$Proposal)
  unloadNamespace("ergm.count")
  expect_false("DiscTNT" %in% ergm_proposal_table()$Proposal)
})

# Ensure that ergm.count is still attached after this test is finished.
library(ergm.count)
