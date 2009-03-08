install.packages("statnet")
library("statnet")

data("faux.magnolia.high")
fmh <- faux.magnolia.high
fmh
plot(fmh, displayisolates = FALSE)
table(component.dist(fmh)$csize)
summary(fmh)
plot(fmh, displayisolates = FALSE, vertex.col = "Grade")
fmh.degreedist <- table(degree(fmh, cmode = "indegree"))
fmh.degreedist
summary(fmh ~ degree(0:8))
help(package = "sna")
summary(fmh ~ degree(0:8, "Sex"))
summary(fmh ~ triangle)
summary(fmh ~ edges + triangle)
mixingmatrix(fmh, "Grade")
gr <- fmh %v% "Grade"
table(gr)
save.image()

model1 <- ergm(fmh ~ edges)
summary(model1)
names(model1)
model1$coef
model1$mle.lik
model2 <- ergm(fmh ~ edges + nodematch("Grade") + nodematch("Race") +
   nodematch("Sex"))
summary(model2)
model2$mle.lik
sim2 <- simulate(model2, burnin = 1e+6, verbose = TRUE, seed = 9)
mixingmatrix(sim2, "Race")
mixingmatrix(fmh, "Race")
plot(summary(fmh ~ degree(0:10)), type = "l", lty = 1, lwd = 2,
   xlab = "Degree", ylab = "Count")
lines(summary(sim2 ~ degree(0:10)), lty = 2, lwd = 3)
legend("topright", legend = c("Observed", "Simulated"), lwd = 3, 
   lty = 1:2)
c(fmh = summary(fmh ~ triangle), sim2 = summary(sim2 ~ triangle))

model3 <- ergm(fmh ~ edges + triangle, seed = 99)
pdf("model3diagnostics.pdf")
mcmc.diagnostics(model3)
dev.off()
model3 <- ergm(fmh ~ edges + triangle, verbose = TRUE, seed = 99)
model3.take2 <- ergm(fmh ~ edges + triangle, MCMCsamplesize = 1e+5,
   interval = 1000, verbose = TRUE, seed = 88)
model3.take3 <- ergm(fmh ~ edges + triangle, maxit = 25, seed = 888,
   control = control.ergm(steplength = 0.25), verbose = TRUE)

model4.take1 <- ergm(fmh ~ edges + nodematch("Grade") +
     nodematch("Race") + nodematch("Sex") + gwesp(0, fixed = TRUE),
     MCMCsamplesize = 1e+5, maxit = 15, verbose = TRUE,
     control = control.ergm(steplength = 0.25), seed = 123)
model4.take1$coef
model4.take2 <- ergm(fmh ~ edges + nodematch("Grade") +
     nodematch("Race") + nodematch("Sex") + gwesp(0.1, fixed = TRUE),
     MCMCsamplesize = 1e+5, maxit = 15, verbose = TRUE,
     control = control.ergm(steplength = 0.25), seed = 123)
model4.take3 <- ergm(fmh ~ edges + nodematch("Grade") +
     nodematch("Race") + nodematch("Sex") + gwesp(0.2, fixed = TRUE),
     MCMCsamplesize = 1e+5, maxit = 15, verbose = TRUE,
     control = control.ergm(steplength = 0.25), seed = 123)
c(model4.take1$mle.lik, model4.take2$mle.lik, model4.take3$mle.lik)
model4 <- model4.take3
model4$coef

sim4 <- simulate(model4, burnin = 1e+5, interval = 1e+5,
     nsim = 100, verbose = TRUE, seed = 321)
class(sim4)
names(sim4)
class(sim4$networks[[1]])
model4.tridist <- sapply(sim4$networks, function(x) summary(x ~ triangle))
hist(model4.tridist, xlab = "Triangles")
fmh.tri <- summary(fmh ~ triangle)
arrows(fmh.tri, 20, fmh.tri, 5, col = "red", lwd = 3)
sum(fmh.tri <= model4.tridist)
gof4.deg <- gof(model4 ~ degree, verbose = TRUE, burnin = 1e+5,
     interval = 1e+5, seed = 246)
plot(gof4.deg)
gof4.deg
gof4.esp.dist <- gof(model4 ~ espartners + distance, verbose = TRUE,
     burnin = 1e+5, interval = 1e+5, seed = 642)
get(getOption("device"))(width = 8, height = 4)
par(mfrow = c(1, 2))
plot(gof4.esp.dist, plotlogodds = TRUE)


