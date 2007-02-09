"anova.ergmlist" <-
  function (object, ..., scale = 0, test = "F") 
{
  objects <- list(object, ...)
  responses <- as.character(lapply(objects, function(x) deparse(x$formula[[2]])))
  sameresp <- responses == responses[1]
  if (!all(sameresp)) {
    objects <- objects[sameresp]
    warning("Models with response ", deparse(responses[!sameresp]), 
            " removed because response differs from ", "model 1")
  }
  nmodels <- length(objects)
  if (nmodels == 1) 
    return(anova.ergm(object))
  n <- network.size(object$newnetwork)
  logl <- df <- Rdf <- rep(0, nmodels)
  pmodemean2 <- rep(NA, nmodels)
  for (i in 1:nmodels) {
    if (!is.null(objects[[i]]$Z.mkl)) {
      Z.use <- objects[[i]]$Z.mkl
      if(!is.null(objects[[i]]$cluster))
        Z.use <- Z.use * apply(objects[[i]]$Beta,2,mean)[ncol(objects[[i]]$Beta)]
      z.pmode.mean <- apply(Z.use, 2, mean)
      Z.pmode <- sweep(Z.use, 2, z.pmode.mean)
      normnorm <- function(x) sum(x^2)
      pmodemean2[i] <- mean(apply(Z.pmode, 1, normnorm))
      p <- ncol(objects[[i]]$Z.mkl)
    }else {
      p <- 0
    }
    nodes<- network.size(objects[[i]]$newnetwork)
    n <- network.dyadcount(objects[[i]]$newnetwork)
    df[i] <- length(objects[[i]]$coef) + (nodes -  (p + 1)/2)*p
    if(!is.null(objects[[i]]$cluster))
      df[i] <- length(objects[[i]]$coef) + objects[[i]]$ngroups*(p+2) - 1 # ng-1 + ng *p + ng
    Rdf[i] <- n - df[i]
    logl[i] <- objects[[i]]$mle.lik
  }
  k <- nmodels
# k <- 1 + length(objects[[i]]$glm$coef)
#
# if (k >= 2) {
#    k <- k+1
#    if(length(object$glm$coef) > 3)
#      varlist <- attr(object$terms, "variables")
#    x <- if (n <- match("x", names(object$glm), 0))
#      object$glm[[n]]
#    else
#      model.matrix(object$glm)
#    varseq <- attr(x, "assign")
#    nvars <- max(0, varseq)
#    resdev <- resdf <- NULL
#    if(nvars>1)
#      for(i in 1:(nvars - 1))
#      {
#        fit <- glm.fit(x = x[, varseq <= i, drop = FALSE], 
#                     y = object$glm$y, weights = object$prior.weights, 
#                     start = object$glm$start, offset = object$glm$offset, 
#                     family = object$glm$family, control = object$glm$control)
#        resdev <- c(resdev, fit$deviance)
#        resdf <- c(resdf, fit$df.residual)
#      }
#
#    df <- c(0, object$glm$df.null - object$glm$df.residual, df)
#    Rdf <- c(object$glm$df.null, resdf,object$glm$df.residual, Rdf)
#    df <- n - Rdf
#    if(length(resdev>0))
#      logl <- c(-object$glm$null.deviance/2, -resdev/2,-object$glm$deviance/2, logl)
#    else logl <- c(-object$glm$null.deviance/2, -object$glm$deviance/2, logl)
#  } else {
    df <- c(0, df)
#   Rdf <- c(object$glm$df.null, Rdf)
#   logl <- c(-object$glm$null.deviance/2, logl)
    Rdf <- c(n, Rdf)
    logl <- c(-n*log(2), logl)
#  }
  var.ex <- 1 - pmodemean2[-1]/pmodemean2[-length(pmodemean2)]
  pv <- pchisq(abs(2 * diff(logl)), abs(diff(df)), lower.tail = FALSE)

  table <- data.frame(c(NA, pmodemean2), c(rep(NA, 2), var.ex),
                      c(NA, -diff(Rdf)), c(NA, diff(2 * logl)), 
                      Rdf, -2 * logl, c(NA, pv))
  variables <- lapply(objects, function(x) paste(deparse(formula(x)), 
                                                 collapse = "\n"))
  colnames(table) <- c("RSS", "R^2","Df","Deviance", "Resid. Df",
                              "Resid. Dev", "Pr(>|Chisq|)")
  if (k > 2) 
    rownames(table) <- c("NULL", object$glm.names,1:nmodels)
  else
    rownames(table) <- c("NULL", 1:nmodels)

  title <- "Analysis of Variance Table\n"
  topnote <- paste("Model ", format(1:nmodels), ": ", variables, 
                   sep = "", collapse = "\n")
  structure(table, heading = c(title, topnote), class = c("anova", 
                                                  "data.frame"))
}
