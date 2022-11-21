test_that("proposal table is purged when package is unloaded", {
  library(ergm.count)
  expect_true("DiscTNT" %in% ergm_proposal_table()$Proposal)
  unloadNamespace("ergm.count")
  expect_false("DiscTNT" %in% ergm_proposal_table()$Proposal)
})

# Ensure that ergm.count is still attached after this test is finished.
library(ergm.count)
