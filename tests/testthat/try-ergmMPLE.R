
data("faux.mesa.high")
summary(faux.mesa.high ~ edges + offset(nodematch("Sex", diff=FALSE)))
##              edges offset(nodematch.Sex) 
##                203                   132 
fit <- ergm(faux.mesa.high ~ edges +  offset(nodematch("Sex", diff=FALSE)), offset.coef = 2)
predict(fit)

as.list(body(ergmMPLE))

trace(
  "ergmMPLE",
  quote(print(str(pl))),
  at = 14
)

trace(
  "ergmMPLE",
  quote(assign("lastpl", pl, envir=.GlobalEnv)),
  at = 14
)


ergmMPLE(
  faux.mesa.high ~ edges + offset(nodematch("Sex", diff=FALSE)),
  output = "fit",
  offset.coef = 2
)

# trace(
#   "ergmMPLE",
#   tracer = quote({
#     cat("Value of pl$theta.offset:\n")
#     print(pl$theta.offset)
#     }),
#   at = 14
# )
# untrace("ergmMPLE")

debugonce(ergmMPLE)

frm <- faux.mesa.high ~ edges + offset(nodematch("Sex", diff=FALSE)) + offset(nodefactor("Sex"))
m <- ergmMPLE(
  frm,
  theta.offset = c(F,T,T), # which parameters are offsets
  output = "matrix"
)
fit <- with(m,
            glm(
              response ~ . - 1,
              data = as.data.frame(predictor),
              weights = weights,
              family = "binomial"
            )
)
