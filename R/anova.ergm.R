"anova.ergm" <- function (object, ...) 
{
  if (length(list(object, ...)) > 1) 
    return(anova.ergmlist(object, ...))
  logl <- object$mle.lik
  nodes<- network.size(object$newnetwork)
# if(is.directed(object$newnetwork)){
#  n <- nodes * (nodes-1)
# }else{
#  n <- nodes * (nodes-1)/2
# }
  n<- network.dyadcount(object$network)
  if(!is.null(object$Z.mkl)){
    p <- ncol(object$Z.mkl)
    Z.use <- object$Z.mkl
    if(!is.null(object$cluster))
      Z.use <- Z.use * apply(object$Beta,2,mean)[ncol(object$Beta)]
    z.pmode.mean <- apply(Z.use, 2, mean)
    Z.pmode <- sweep(Z.use, 2, z.pmode.mean)
    normnorm <- function(x) sum(x^2)
    pmodemean2 <- mean(apply(Z.pmode, 1, normnorm))
  }else{
    pmodemean2 <- NA
    p <- 0
  }
  df <- length(object$coef) + (nodes - (p + 1)/2) * p
  if(!is.null(object$cluster))
    df <- length(object$coef) + object$ngroups*(p+2) - 1 # ng-1 + ng *p + ng
  Rdf <- n - df
#
  k <- 1 + (length(object$mplefit$glm$coef) >= 2)
#  if (k > 1 & is.latent(object)) {
#    if(length(object$mplefit$glm$coef) > 3)
#      varlist <- attr(object$terms, "variables")
#    x <- if (n <- match("x", names(object$mplefit$glm), 0))
#      object$mplefit$glm[[n]]
#    else
#      model.matrix(object$mplefit$glm)
#    varseq <- attr(x, "assign")
#    nvars <- max(0, varseq)
#    resdev <- resdf <- NULL
#    for(i in 1:(nvars - 1))
#    {
#      fit <- glm.fit(x = x[, varseq <= i, drop = FALSE], 
#                     y = object$mplefit$glm$y, weights = object$prior.weights, 
#                     start = object$mplefit$glm$start, offset = object$mplefit$glm$offset, 
#                     family = object$mplefit$glm$family, control = object$mplefit$glm$control)
#      resdev <- c(resdev, fit$deviance)
#      resdf <- c(resdf, fit$df.residual)
#    }
#
#    df <- c(0, object$mplefit$glm$df.null - object$mplefit$glm$df.residual, df)
#    Rdf <- c(object$mplefit$glm$df.null, resdf,object$mplefit$glm$df.residual, Rdf)
#    df <- n - Rdf
#    logl <- c(-object$mplefit$glm$null.deviance/2, -resdev/2,-object$mplefit$glm$deviance/2, logl)
#    pmm <- c(rep(NA,nvars),NA,pmodemean2)
#  }else{
    df <- c(0, df)
#   Rdf <- c(object$mplefit$glm$df.null, Rdf)
#   logl <- c(-object$mplefit$glm$null.deviance/2, logl)
    Rdf <- c(n, Rdf)
    logl <- c(-n*log(2), logl)
    pmm <- c(NA,pmodemean2)
#  }
  pv <- pchisq(abs(2 * diff(logl)), abs(diff(df)), lower.tail = FALSE)
  table <- data.frame(pmm,c(NA, -diff(Rdf)), c(NA, diff(2 * logl)), 
                      Rdf, -2 * logl, c(NA, pv))
  variables <- paste(deparse(formula(object)), collapse = "\n")
  colnames(table) <- c("RSS", "Df", "Deviance", "Resid. Df", "Resid. Dev", 
                       "Pr(>|Chisq|)")
#  if(k>1 & is.latent(object)){
#    if(!is.null(object$cluster)){
#     rownames(table) <- c("NULL",names(object$coef), "Latent Cluster")
#    }else{
#     rownames(table) <- c("NULL",names(object$coef), "Latent")
#    }
#  }else{
    rownames(table) <- c("NULL", "Model 1:")
#  }
  title <- "Analysis of Variance Table\n"
  topnote <- paste("Model ", format(1), ": ", variables, sep = "", 
                   collapse = "\n")
  structure(table, heading = c(title, topnote), class = c("anova", 
                                                  "data.frame"))
}
