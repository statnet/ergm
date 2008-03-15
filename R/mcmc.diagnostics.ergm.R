mcmc.diagnostics <- function(object, ...) {
  UseMethod("mcmc.diagnostics")
}

mcmc.diagnostics.default <- function(object, ...) {
  stop("An object must be given as an argument ")
}

mcmc.diagnostics.ergm <- function(object, sample="sample",
                                  smooth=TRUE,
                                  r=0.0125, digits=6,
                                  maxplot=1000, verbose=TRUE, center=TRUE,
                                  main="Summary of MCMC samples",  
                                  xlab = "Iterations", ylab = "", ...) {
#
  if(!is.null(object$degeneracy.value) && !is.na(object$degeneracy.value)){
   degeneracy.value <- object$degeneracy.value
   degeneracy.type <- object$degeneracy.type
  }else{
   degeneracy.value <- NULL
   degeneracy.type <- NULL
  }

  if(sample=="missing"){
    component <- "sample"
    statsmatrix.miss <- object[["conditionalsample"]]
    if(missing(mcmc.title)){mcmc.title <- "Summary of the Conditional Samples"} 
  }else{
    component <- sample
  }
  if(component=="sample"&&is.null(object$sample)){
    stop("There is no MCMC sample associated with the object.\n")
  }
  statsmatrix <- object[[component]]
  if(!is.matrix(statsmatrix) || length(dim(statsmatrix))==0){
    stop("There is no ",component," component of the object.\n")
  }
  if(component=="thetasample"){
    if(is.null(object$MCMCtheta)){
     x0 <- rep(0,ncol(statsmatrix))
    }else{
     x0 <- object$MCMCtheta
    }
#   statsmatrix <- sweep(statsmatrix,2,x0,"-")
  }else{
    if(is.null(object$formula) ){
     x0 <- rep(0,ncol(statsmatrix))
    }else{
#     if(!is.latent(object) ){
      if(sample=="missing"){
        x0 <- apply(statsmatrix.miss,2,mean)
      }else{
        x0 <- rep(0,ncol(statsmatrix))
      }
      if(!center){
       x0 <- try(summary(object$formula), silent=TRUE)
       if(!inherits(x0,"try-error")){
        statsmatrix <- sweep(statsmatrix,2,x0,"+")
       }else{
        x0 <- rep(0,ncol(statsmatrix))
       }
      }else{
       if(sample=="missing"){
        statsmatrix <- sweep(statsmatrix,2,x0,"-")
        x0 <- rep(0,ncol(statsmatrix))
       }
      }
      if(ncol(statsmatrix) < length(x0)){
       x0 <- x0[-c(1:(length(x0)-ncol(statsmatrix)))]
      }
#     }
    }
  }
  attributes(statsmatrix)$class <- NULL
  novar <- apply(statsmatrix,2,var)<1e-6
  if(all(novar)){
     warning("All the statistics are the same.\n")
     print(apply(statsmatrix,2,summary.statsmatrix.ergm),scipen=6)
     return(invisible())
  }else{
    statsmatrix <- as.matrix(statsmatrix[,!novar,drop=FALSE])
    x0 <- x0[!novar]
    colnames(statsmatrix) <- colnames(object[[component]])[!novar]

    if(verbose){
     cat("\nCorrelations of sample statistics:\n")
     if(is.null(object$acf)){
      print(ergm.MCMCacf(statsmatrix))
     }else{
      print(object$acf)
     }
    }

    attr(statsmatrix, "mcpar") <- attr(object[[component]], "mcpar")
    if(is.null(attr(statsmatrix, "mcpar"))){
      attr(statsmatrix, "mcpar") <- c(1,nrow(statsmatrix),1)
    }
    attr(statsmatrix, "class") <- "mcmc"
    if(require("coda", quietly = TRUE)) {
     plot.mcmc.ergm(statsmatrix, ask=FALSE, smooth=smooth, 
                    maxplot=maxplot, parallel=object$parallel,
                    x0=x0, 
                    xlab=xlab, ylab=ylab, mcmc.title=main, ...)
    }else{
     warning("For all MCMC diagnostics you need the 'coda' package.")
     return(invisible())
    }
    if(is.null(degeneracy.value) || !is.infinite(degeneracy.value)){
     cat("\nr=0.0125 and 0.9875:\n")
     raft9875 <- ergm.raftery.diag(statsmatrix, r=0.9875, ...)
     raft     <- ergm.raftery.diag(statsmatrix, r=0.0125, ...)
     aaa <- raft9875$resmatrix > raft$resmatrix
     aaa[is.na(aaa)] <- FALSE
     raft$resmatrix[aaa] <- raft9875$resmatrix[aaa]
     simvalues <- attr(statsmatrix, "mcpar")
     if(is.null(simvalues)){
       simvalues <- c(2, nrow(statsmatrix), 1)
     }
     raft$degeneracy.value <- degeneracy.value
     raft$degeneracy.type <- degeneracy.type
     if(verbose){
      print(raft, simvalues=simvalues)
     }
    }else{
     raft <- list(simvalues=NULL,
                  degeneracy.value=degeneracy.value,
                  degeneracy.type=degeneracy.type)
    }
    return(invisible(raft))
  }
}
"ergm.raftery.diag" <-
function (data, q = 0.025, r = 0.005, s = 0.95, converge.eps = 0.001) 
{
    if (is.mcmc.list(data)) 
        return(lapply(data, raftery.diag, q, r, s, converge.eps))
#   data <- as.mcmc(data)
    resmatrix <- matrix(nrow = nvar(data), ncol = 4,
        dimnames = list(varnames(data, 
        allow.null = TRUE), c("M", "N", "Nmin", "I")))
    phi <- qnorm(0.5 * (1 + s))
    nmin <- as.integer(ceiling((q * (1 - q) * phi^2)/r^2))
    if (nmin > nrow(data)) 
      resmatrix <- c("Error", nmin)
    else for (i in 1:nvar(data)) {
      #          First need to find the thinning parameter kthin 
      # 
      if (is.matrix(data)) {
        quant <- quantile(data[, i, drop = TRUE], probs = q)
        dichot <- mcmc(data[, i, drop = TRUE] <= quant,
                       start = 1, end = nrow(data), thin = 1)
      }
      else {
          quant <- quantile(data, probs = q)
          dichot <- mcmc(data <= quant, start = 1,
                         end = nrow(data), thin = 1)
      }
      attr(dichot, "mcpar") <- attr(data, "mcpar")
      kthin <- 0
      bic <- 1
      use <- 1:nrow(data)
      while (bic >= 0) {
        kthin <- kthin + 1
#       testres <- as.vector(ergm.window.mcmc(dichot, thin = kthin))
        testres <- dichot[use <= trunc((length(dichot) - 1)/kthin + 1.5)]
        newdim <- length(testres)
#       testtran <- table(testres[1:(newdim - 2)],
#                         testres[2:(newdim - 1)],
#                         testres[3:newdim])
        aaa <- testres[1:(newdim - 2)] + testres[2:(newdim - 1)]*2 +testres[3:newdim]*4 + 1
        testtran <- array(tabulate(aaa, nbin=9), dim=c(2,2,2),
            dimnames=list(c(TRUE,FALSE), c(TRUE,FALSE),c(TRUE,FALSE)))
        testtran <- array(as.double(testtran), dim = dim(testtran))
        g2 <- 0
        for (i1 in 1:2) {
          for (i2 in 1:2) {
            for (i3 in 1:2) {
              if (testtran[i1, i2, i3] != 0) {
                fitted <- (sum(testtran[i1, i2, 1:2]) * 
                  sum(testtran[1:2, i2, i3]))/(sum(testtran[1:2, i2, 1:2]))
                g2 <- g2 + testtran[i1, i2, i3] * log(testtran[i1, 
                  i2, i3]/fitted) * 2
              }
            }
          }
        }
        bic <- g2 - log(newdim - 2) * 2
      }
      #
      # then need to find length of burn-in and No of iterations 
      # for required precision 
      # 
#     finaltran <- table(testres[1:(newdim - 1)], testres[2:newdim])
      aaa <- testres[1:(newdim - 1)] + testres[2:newdim]*2 + 1
      finaltran <- array(tabulate(aaa, nbin=4), dim=c(2,2),
            dimnames=list(c(TRUE,FALSE), c(TRUE,FALSE)))
      alpha <- finaltran[1, 2]/(finaltran[1, 1] + finaltran[1, 2])
      beta <- finaltran[2, 1]/(finaltran[2, 1] + finaltran[2, 2])
      tempburn <- log((converge.eps * (alpha + beta))/max(alpha, 
          beta))/(log(abs(1 - alpha - beta)))
      tempburn[is.na(tempburn)] <- 0 
      nburn <- as.integer(ceiling(tempburn) * kthin * thin(data))
      tempprec <- ((2 - alpha - beta) * alpha * beta * phi^2)/(((alpha + 
          beta)^3) * r^2)
      tempprec[is.na(tempprec)] <- 0 
      nkeep <- as.integer(ceiling(tempprec) * kthin * thin(data))
      iratio <- (nburn + nkeep)/nmin
      resmatrix[i, 1] <- nburn
      resmatrix[i, 2] <- nkeep + nburn
      resmatrix[i, 3] <- nmin
      resmatrix[i, 4] <- signif(iratio, digits = 3)
    }
    resmatrix[is.na(resmatrix)] <- 0
    y <- list(params = c(r = r, s = s, q = q), resmatrix = resmatrix)
    class(y) <- "raftery.diag.ergm"
    return(y)
}
print.raftery.diag.ergm <- function (x, digits = 3, simvalues=NULL, ...) 
{
    if(is.null(simvalues)){
      simvalues <- c(2, nrow(x$resmatrix), 1)
    }
    cat("\nQuantile (q) =", x$params["q"])
    cat("\nAccuracy (r) = +/-", x$params["r"])
    cat("\nProbability (s) =", x$params["s"], "\n")
    if (!is.na(x$resmatrix[1]) && x$resmatrix[1] == "Error") 
        cat("\nYou need a sample size of at least", x$resmatrix[2], 
            "with these values of q, r and s\n")
    else {
        out <- cbind(x$resmatrix, 
          sweep(matrix(x$resmatrix[,1:2],nrow=nrow(x$resmatrix)),2,simvalues[1:2],"<"))
        for (i in ncol(out)) out[, i] <- format(out[, i], digits = digits)
        out <- rbind(matrix(c("Burn-in ", "Total", "Lower bound ", 
            "Dependence", "enough", "enough", "(M)", "(N)", "(Nmin)",
            "factor (I)", "burn-in?", "samples?"), 
            byrow = TRUE, nrow = 2), out)
        if (!is.null(rownames(x$resmatrix))){ 
            out <- cbind(c("", "", rownames(x$resmatrix)), out)
        }else{
            out <- cbind("", out)
        }
        dimnames(out) <- list(rep("", nrow(out)), rep("", ncol(out)))
        if(ncol(out)>6){
         out[out[,7]=="1",7] <- "yes"
         out[out[,7]=="0",7] <- " no"
         out[out[,2]=="0",7] <- " no"
         out[out[,3]=="0",7] <- " no"
         out[out[,5]=="0",7] <- " no"
        }
        if(ncol(out)>6){
         out[out[,6]=="1",6] <- "yes"
         out[out[,6]=="0",6] <- " no"
         out[out[,2]=="0",6] <- " no"
         out[out[,3]=="0",6] <- " no"
         out[out[,5]=="0",6] <- " no"
        }
        print.default(out, quote = FALSE, ...)
        cat("\n")
    }
    invisible(x)
}
"plot.mcmc.ergm" <- function(x, trace = TRUE, density = TRUE, 
                         smooth = TRUE, bwf, 
                         auto.layout = TRUE, ask = TRUE,
                         xlab = "Iterations", ylab = "",
                         maxplot=1000, parallel=0, x0, mcmc.title="", ...) 
{
  oldpar <- NULL
  on.exit(par(oldpar))
  if (auto.layout) {
    mfrow <- set.mfrow(Nchains = nchain(x), Nparms = nvar(x), 
                       nplots = trace + density)
    oldpar <- par(mfrow = mfrow)
  }
  oldpar <- c(oldpar, par(ask = ask))
  for (i in 1:nvar(x)) {
    y <- x[, i, drop = FALSE]
    attr(y, "class") <- "mcmc"
    if (trace){ 
      traceplot.ergm(y, smooth = smooth, maxplot=maxplot, parallel=parallel,
                     xlab=xlab, ylab=ylab)
      abline(h=x0[i],lty=3, lwd=2)
    }
    if (density){
      if (missing(bwf)) 
        densplot(y)
      else{densplot(y, bwf = bwf)}
      abline(v=x0[i],lty=3, lwd=2)
    }
    if(i==1){mtext(text=mcmc.title, outer=TRUE, line = -1.5)}
  }
}
"traceplot.ergm" <-
function (x, smooth = TRUE, col = 1:6, type = "l",
          xlab = "Iterations", ylab = "", maxplot=1000, parallel=0, ...) 
{
# x <- mcmc.list(x)
# xmcpar <- c(start(x), end(x), thin(x))
  args <- list(...)
  for (j in 1:nvar(x)) {
    xp <- seq(start(x), end(x), thin(x))
    yp <- x
    if(length(yp) > maxplot){
     yp <- as.matrix(yp[seq(1,length(yp),length=maxplot),])
     xp <- xp[seq(1,length(xp),length=maxplot)]
    }
    matplot(xp, yp, xlab = xlab, ylab = ylab, type = type, 
         col = col, ...)
    if(!is.null(parallel) && parallel>0){
     abline(v=xp[round(seq(from=1,to=length(xp),length=parallel+1))],lty=2)
    }
    if (!is.null(varnames(x)) && is.null(list(...)$main)) 
      title(paste("Trace of", varnames(x)[j]))
    if (smooth) {
      scol <- rep(col, length = nchain(x))
      for (k in 1:nchain(x)) lines(lowess(xp, yp[, k]), 
                                   col = scol[k])
    }
  }
}
"set.mfrow" <-
function (Nchains = 1, Nparms = 1, nplots = 1, sepplot = FALSE) 
{
  ## Set up dimensions of graphics window: 
  ## If only density plots OR trace plots are requested, dimensions are: 
  ##	1 x 1	if Nparms = 1 
  ##	1 X 2 	if Nparms = 2 
  ##	2 X 2 	if Nparms = 3 or 4 
  ##	3 X 2 	if Nparms = 5 or 6 or 10 - 12 
  ##	3 X 3 	if Nparms = 7 - 9 or >= 13 
  ## If both density plots AND trace plots are requested, dimensions are: 
  ##	1 x 2	if Nparms = 1 
  ##	2 X 2 	if Nparms = 2 
  ##	3 X 2 	if Nparms = 3, 5, 6, 10, 11, or 12 
  ##	4 x 2	if Nparms otherwise 
  ## If separate plots are requested for each chain, dimensions are: 
  ##	1 x 2	if Nparms = 1 & Nchains = 2 
  ##	2 X 2 	if Nparms = 2 & Nchains = 2 OR Nparms = 1 & Nchains = 3 or 4 
  ##	3 x 2	if Nparms = 3 or >= 5 & Nchains = 2  
  ##		   OR Nchains = 5 or 6 or 10 - 12 (and any Nparms) 
  ##	2 x 3	if Nparms = 2 or 4 & Nchains = 3 
  ##	4 x 2   if Nparms = 4 & Nchains = 2  
  ##		   OR Nchains = 4 & Nparms > 1 
  ##	3 x 3	if Nparms = 3 or >= 5  & Nchains = 3  
  ##		   OR Nchains = 7 - 9 or >= 13 (and any Nparms)
  mfrow <- if (sepplot && Nchains > 1 && nplots == 1) {
    ## Separate plots per chain
    ## Only one plot per variable
    if (Nchains == 2) {
      switch(min(Nparms, 5),
             c(1,2),
             c(2,2),
             c(3,2),
             c(4,2),
             c(3,2))
    }
    else if (Nchains == 3) {
      switch(min(Nparms, 5),
             c(2,2),
             c(2,3),
             c(3,3),
             c(2,3),
             c(3,3))
    }
    else if (Nchains == 4) {
      if (Nparms == 1)
        c(2,2)
      else
        c(4,2)
    }
    else if (any(Nchains == c(5,6,10,11,12)))
      c(3,2)
    else if (any(Nchains == c(7,8,9)) || Nchains >=13)
      c(3,3)
      
  }
  else {
    if (nplots==1) {
      ## One plot per variable
      mfrow <- switch(min(Nparms,13),
                      c(1,1),
                      c(1,2),
                      c(2,2),
                      c(2,2),
                      c(3,2),
                      c(3,2),
                      c(3,3),
                      c(3,3),
                      c(3,3),
                      c(3,2),
                      c(3,2),
                      c(3,2),
                      c(3,3))
    }
    else {
      ## Two plot per variable
      ##
      mfrow <- switch(min(Nparms, 13),
                      c(1,2),
                      c(2,2),
                      c(3,2),
                      c(4,2),
                      c(3,2),
                      c(3,2),
                      c(4,2),
                      c(4,2),
                      c(4,2),
                      c(3,2),
                      c(3,2),
                      c(3,2),
                      c(4,2))
    }
  }
  return(mfrow)
}

#"nvar.mcmc" <- function (x) {
#  if (is.mcmc.object(x)) {
#    if (is.matrix(x)){ncol(x)}else{1}
#  }
#  else if (is.mcmc.list(x)) {
#    if (is.matrix(x[[1]])){ncol(x[[1]])}else{1}
#  }
#  else NULL
#}
##
#"is.mcmc.object" <- function (x){
#  inherits(x, "mcmc")
#}
#
#"is.mcmc.list.object" <- function (x){
#  inherits(x, "mcmc.list")
#}
#
#"varnames.mcmc" <- function (x, allow.null = TRUE)
#{
#  if (!is.mcmc.object(x) && !is.mcmc.list.object(x))
#    return(NULL)
#  y <- if (is.mcmc.object(x))
#    dimnames(x)[[2]]
#  else if (is.mcmc.list.object(x))
#    dimnames(x[[1]])[[2]]
#  if (is.null(y) && !allow.null)
#    y <- paste("var", 1:nvar(x), sep = "")
#  return(y)
#}
