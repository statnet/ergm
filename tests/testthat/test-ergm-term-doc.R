#  File tests/testthat/test-ergm-term-doc.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################

library(ergm.count)

test_that("check initialisation functions", {
  for (term in ergm:::ergmTermCache('ergmTerm')) {
    # every term must include either 'valued' or 'binary'
    # check that there is a visable init function defined for the term
    # some terms have both valued an binary forms
    if ('valued'%in%term$concepts){
      expect_true(is.function(eval(locate_prefixed_function(term$name, 'InitWtErgmTerm'))))
    }
    if ('binary'%in%term$concepts){
      expect_true(is.function(eval(locate_prefixed_function(term$name, 'InitErgmTerm'))))
    }
  }
})

test_that("test search ergm term", {
  # crude checks for search.ergmTerms are in the search.ergmTerms man page

  # expect to find at least eight terms mentioning triangles
  expect_equal(length(search.ergmTerms('triangle', packages='ergm')), 9)

  # search using a bipartite net as a template
  myNet<-network.initialize(5,bipartite=3,directed=FALSE)
  expect_equal(length(search.ergmTerms(net=myNet)), 31)

  expect_equal(length(search.ergmTerms(keywords = 'bipartite', packages='ergm')), 36)

  expect_gt(length(search.ergmTerms(name = 'b2factor', packages='ergm')), 0)
  expect_equal(length(search.ergmTerms(name = 'b3factor', packages='ergm')), 0)

  expect_equal(length(search.ergmTerms(keywords = 'bipartite', packages='ergm')), 36)

  ## expect_gt(length(search.ergmTerms(keywords = 'valued')), 44)
  expect_equal(length(search.ergmTerms(keywords = 'valued', packages='ergm')), 44)
  ## expect_gt(length(search.ergmTerms(keywords = 'valued', packages=c('ergm', 'ergm.count'))), 44)
})

test_that("test search ergm reference", {
  expect_equal(length(search.ergmReferences('dyad', packages='ergm')), 4)

  ## expect_equal(length(search.ergmReferences(keywords = 'binary')), 1)
  expect_equal(length(search.ergmReferences(keywords = 'binary', packages='blah')), 0)
  expect_equal(length(search.ergmReferences(keywords = 'binary', packages='ergm')), 1)

  expect_gt(length(search.ergmReferences(name = 'Bernoulli', packages='ergm')), 0)
  expect_equal(length(search.ergmReferences(name = 'Cernoulli', packages='ergm')), 0)
})

test_that("test search ergm constraint", {
  expect_equal(length(search.ergmConstraints('degree', packages='ergm')), 9)

  ## expect_equal(length(search.ergmConstraints(keywords = 'directed')), 16)
  expect_equal(length(search.ergmConstraints(keywords = 'directed', packages='blah')), 0)
  expect_equal(length(search.ergmConstraints(keywords = 'directed', packages='ergm')), 16)

  expect_gt(length(search.ergmConstraints(name = 'b1degrees', packages='ergm')), 0)
  expect_equal(length(search.ergmConstraints(name = 'b3degrees', packages='ergm')), 0)
})

test_that("test search ergm proposal", {
  expect_equal(length(search.ergmProposals('bipartite', packages='ergm')), 2)

  expect_equal(length(search.ergmProposals(constraints='.dyads', packages='ergm')), 5)

  ## expect_equal(length(search.ergmProposals(reference='Bernoulli')), 17)
  expect_equal(length(search.ergmProposals(reference='Bernoulli', packages='ergm.count')), 0)
  expect_equal(length(search.ergmProposals(reference='Bernoulli', packages='ergm')), 18)

  expect_equal(length(search.ergmProposals(name = 'randomtoggle', packages='ergm')), 1)
  expect_equal(length(search.ergmProposals(name = 'mandomtoggle', packages='ergm')), 0)
})

test_that("ergm term cache unloading", {
  library(ergm.count)
  et <- ergm:::ergmTermCache('ergmTerm')
  expect_true("ergm.count" %in% sapply(et, `[[`, "package"))

  unloadNamespace("ergm.count")
  et <- ergm:::ergmTermCache('ergmTerm')
  expect_false("ergm.count" %in% sapply(et, `[[`, "package"))
})

# Ensure that ergm.count is still attached after this test is finished.
library(ergm.count)
