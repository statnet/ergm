#  File tests/ergm_term_doc_test.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################

# tests to veryify the structure of the ergm-terms documentation
# this executes functions similar to what would be run when building the ergm=term=crossRef.Rmd file
# and verifys that the structure of ergm-terms.Rd matches what the parser expects
test_that("check terms keywords", {
  # every term must include at least one of 'directed', 'undirected', or 'operator'
  for (term in ergm:::ergmTermCache('ergmTerm')) {
    expect_true(any(c('directed','undirected','operator')%in%term$concepts))
    expect_true(any(c('binary','valued')%in%term$concepts))
  }
})

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
  expect_equal(length(search.ergmTerms('triangle')), 9)

  # search using a bipartite net as a template
  myNet<-network.initialize(5,bipartite=3,directed=FALSE)
  expect_equal(length(search.ergmTerms(net=myNet)), 38)

  expect_equal(length(search.ergmTerms(keywords = 'bipartite')), 38)

  expect_equal(length(search.ergmTerms(keywords = 'bipartite', packages='ergm')), 38)

  library(ergm.count)
  expect_equal(length(search.ergmTerms(keywords = 'valued')), 85)
  expect_equal(length(search.ergmTerms(keywords = 'valued', packages='ergm')), 83)
  expect_equal(length(search.ergmTerms(keywords = 'valued', packages=c('ergm', 'ergm.count'))), 85)
})
