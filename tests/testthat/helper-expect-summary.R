#  File tests/testthat/helper-expect-summary.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################
expect_summary <- function(s, e, value, coefficients, tolerance=0.001) {
  expect_equal(s, value, tolerance=tolerance, ignore_attr=is.null(names(value)))
  expect_equal(coef(e)[seq_along(coefficients)], coefficients, tolerance=tolerance, ignore_attr=is.null(names(coefficients)))
}

gf_summary <- function(object, ..., extra_toggles = NULL, tolerance = 0.001) {
  nw <- ergm.getnetwork(object)
  # Avoid reinitializing.
  m <- ergm_model(object, ...)
  s1 <- ergm_state(nw, m, stats = numeric(nparam(m, canonical = TRUE)))
  summ1 <- summary(s1)

  if (length(summ1) == 0L) return(summ1)

  el <- as.edgelist(s1, output = "tibble")
  el0 <- el[0L, ]
  s0 <- update(s1, el = el0)
  summ0 <- summary(s0)

  el <- as.matrix(el)
  el2 <- rbind(el, el)

  # Additional random toggles, always cancelling themselves.
  ele <- if(length(extra_toggles) == 1L && extra_toggles > 0L)
           as.network(network.size(nw), directed = is.directed(nw),
                      bipartite = nw %n% "bipartite", loops = has.loops(nw),
                      numedges = extra_toggles) |> as.edgelist()
         else matrix(as.integer(extra_toggles), ncol = 2)
  ele <- rbind(ele, ele)

  scramble <- function(m) m[sample.int(nrow(m)), , drop = FALSE]
  sort_matrix <- function(m) m[order(m[,1], m[,2]), , drop = FALSE]


  gf <- ergm.godfather(s1,
                       changes = list(rbind(ele, el), # Remove all edges in lex order. -> summ0 - summ1
                                      sort_matrix(rbind(ele, el2)), # Add and remove each edge immediately. -> 0
                                      scramble(rbind(ele, el)), # Add all edges back in random order. -> summ1 - summ0
                                      scramble(rbind(ele, el2))), # Add and remove all edges in random order. -> 0
                       changes.only = TRUE) |> as.matrix()

  expect_equal(diff(rbind(0,gf)),
               rbind(summ0 - summ1, 0, summ1 - summ0, 0),
               tolerance = tolerance)

  summ1 # Return the summary.
}
