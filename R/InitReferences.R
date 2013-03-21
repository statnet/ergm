InitReference.Bernoulli <- function(lhs.nw, ...){
  list(name="Bernoulli")  
}

InitReference.Unif <- function(lhs.nw, a, b, ...){
  list(name="Unif", a=a, b=b)  
}

InitReference.DiscUnif <- function(lhs.nw, a, b, ...){
  list(name="DiscUnif", a=a, b=b)
}
