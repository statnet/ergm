
print.stergm <- function (x, digits = max(3, getOption("digits") - 3), ...) {
  cat("Formation Coefficients:\n")
  print.default(format(x$formation.fit$coef, digits = digits), print.gap = 2, quote = FALSE)
  cat("Dissolution Coefficients:\n")
  print.default(format(x$dissolution.fit$coef, digits = digits), print.gap = 2, quote = FALSE)
}
