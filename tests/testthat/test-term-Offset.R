test_that("Estimation with Offset() operator works", {
  local_edition(3)
  nw <- network(8, dir = F)
  offset_RE <- ".*offset.* decorator used on term .* with no free parameters is meaningless.*"

  expect_message(
    ergm_model(nw ~ edges + offset(Offset(~triangle, which = 1, coef = 0.1))),
    offset_RE)

  expect_failure(expect_message(
    off <- ergm(nw ~ edges + offset(triangle), offset.coef = 0.1),
    offset_RE), ".* did not throw the expected message.*")

  expect_failure(expect_message(
    Off <- ergm(nw ~ edges + Offset(~triangle, which = 1, coef = 0.1)),
    offset_RE), ".* did not throw the expected message.*")

  expect_lt(
  (coef(off)[1] - coef(Off)) /
  sqrt(vcov(off, sources = "estimation")[1,1] + vcov(Off, sources = "estimation")),
  qnorm(.9999))
})