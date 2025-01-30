add_eattr <- function(nw, attr, f = runif){
  nwn <- substitute(nw)
  nw %e% attr <- f(network.edgecount(nw, na.omit = FALSE))

  eval.parent(call("<-", nwn, nw), 1)
}
