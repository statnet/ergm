#  File tests/testthat/helper-htests.R in package ergm, part of the Statnet
#  suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################
expect_within_mc_err <- function(object, expected, idx = TRUE, alpha = 0.001) {
  act <- quasi_label(rlang::enquo(object), arg = "object")

  if(!is(act$val, "ergm")) stop(sQuote("object"), " should be an ", sQuote("ergm"), " fit.")

  x <- coef(act$val)[idx]
  v <- vcov(act$val, source = "estimation")[idx, idx, drop = FALSE]

  chi2 <- xTAx_seigen(x - expected, v)
  df <- attr(chi2, "rank")

  expect((pval <- pchisq(chi2, df, lower.tail = FALSE)) > alpha,
         if(isTRUE(idx)) sprintf("%s: χ² = %f, df = %i, p-val. = %f ≤ %f = α.", act$lab, chi2, df, pval, alpha)
         else sprintf("%s[%s]: χ² = %f, df = %i, p-val. = %f ≤ %f = α.", act$lab, deparse1(idx), chi2, df, pval, alpha))

  invisible(act$val)
}

expect_within_mc_err2 <- function(object1, object2, idx1 = TRUE, idx2 = idx1, alpha = 0.001) {
  obj1 <- quasi_label(rlang::enquo(object1), arg = "object1")
  obj2 <- quasi_label(rlang::enquo(object2), arg = "object2")

  if(!is(obj1$val, "ergm")) stop(sQuote("object1"), " should be an ", sQuote("ergm"), " fit.")
  if(!is(obj2$val, "ergm")) stop(sQuote("object2"), " should be an ", sQuote("ergm"), " fit.")

  x1 <- coef(obj1$val)[idx1]
  v1 <- vcov(obj1$val, source = "estimation")[idx1, idx1, drop = FALSE]
  x2 <- coef(obj2$val)[idx2]
  v2 <- vcov(obj2$val, source = "estimation")[idx2, idx2, drop = FALSE]

  chi2 <- xTAx_seigen(x1 - x2, v1 + v2)
  df <- attr(chi2, "rank")

  expect((pval <- pchisq(chi2, df, lower.tail = FALSE)) > alpha,
         if(isTRUE(idx1) && isTRUE(idx2)) sprintf("%s vs. %s: χ² = %f, df = %i, p-val. = %f ≤ %f = α.", obj1$lab, obj2$lab, chi2, df, pval, alpha)
         else sprintf("%s[%s] vs. %s[%s]: χ² = %f, df = %i, p-val. = %f ≤ %f = α.", obj1$lab, deparse1(idx1), obj2$lab, deparse1(idx2), chi2, df, pval, alpha))

  invisible(obj1$val)
}

expect_t_test <- function(object, expected, alpha = 0.001) {
  act <- quasi_label(rlang::enquo(object), arg = "object")
  exp <- quasi_label(rlang::enquo(expected), arg = "expected")

  tst <- t.test(act$val, mu = exp$val)

  expect(tst$p.value > alpha,
         sprintf("H₀: E(%s) = %s: t = (%f - %f) / %f  = %f, df = %i, p-val. = %f ≤ %f = α.", act$lab, exp$lab, tst$estimate, tst$null.value, tst$stderr, tst$statistic, tst$parameter, tst$p.value, alpha))

  invisible(act$val)
}
