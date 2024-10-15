geodist <- function(M){
  n <- nrow(M)
  R <- M
  D <- matrix(Inf, n, n)
  diag(D) <- 0
  for(d in seq_len(n-1)){
    D[c(R != 0)] <- pmin(D[c(R != 0)], d)
    R <- R %*% M
  }

  D
}


geodistn <- function(M){
  c(table(factor(c(geodist(M)), levels = c(seq_len(nrow(M)-1), Inf))))
}


test_that("geodesic distance calculation (undirected)", {
  data(florentine)
  expect_equal(
    ergm.geodistdist(flomarriage),
    geodistn(as.matrix(flomarriage))/2
  )
})


test_that("geodesic distance calculation (undirected)", {
  data(sampson)
  expect_equal(
    ergm.geodistdist(samplike),
    geodistn(as.matrix(samplike))
  )
})
