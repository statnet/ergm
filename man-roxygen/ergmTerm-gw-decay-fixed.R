#' @param decay nonnegative decay parameter for the <%= multiplicand %>; required if `fixed=TRUE` and ignored with a warning otherwise.
#' @param fixed optional argument indicating
#'   whether the `decay` parameter is fixed at the given value, or is to be fit as a curved
#'   exponential-family model (see Hunter and Handcock, 2006). The
#'   default is `FALSE` , which means the scale parameter is not
#'   fixed and thus the model is a curved exponential family. 
