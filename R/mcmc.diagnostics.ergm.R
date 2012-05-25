#  File ergm/R/mcmc.diagnostics.ergm.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
#########################################################################
# The <mcmc.diagnostics.X> functions create diagnostic plots for the
# MCMC sampled statistics of the ergm X and prints the Raftery-Lewis
# diagnostics, indicating whether they are sufficient or not; if X is not
# an ergm, execution will halt
##########################################################################

mcmc.diagnostics <- function(object, ...) {
  UseMethod("mcmc.diagnostics")
}

mcmc.diagnostics.default <- function(object, ...) {
  stop("An ergm object must be given as an argument ")
}

mcmc.diagnostics.ergm <- function(object,
                                  center=TRUE,
                                  curved=TRUE,
                                  vars.per.page=3,...) {
#
  if(!is.null(object$degeneracy.value) && !is.na(object$degeneracy.value)){
   degeneracy.value <- object$degeneracy.value
   degeneracy.type <- object$degeneracy.type
  }else{
   degeneracy.value <- NULL
   degeneracy.type <- NULL
  }

  # Coerce sample objects to mcmc.list. This allows all subsequent
  # operations to assume mcmc.list. The reason [["sample"]] is being
  # used here rather than $sample is because there is an unlikely
  # posibility that $sample doesn't exist but $sample.obs does.
  sm <- if(is.null(object[["sample"]])) NULL else as.mcmc.list(object[["sample"]])
  sm.obs <- if(is.null(object[["sample.obs"]])) NULL else as.mcmc.list(object[["sample.obs"]])

  if(is.null(sm)) stop("MCMC was not run or MCMC sample was not stored.")

  if(!center){
    sm <- sweep.mcmc.list(sm, object$target.stats, "+")
    if(!is.null(sm.obs)) sm.obs <- sweep.mcmc.list(sm.obs, object$target.stats, "+")
  }else{
    # Then sm is already centered, *unless* there is missing data.  In
    # that case, center sm relative to sm.obs and center sm.obs to
    # 0. (The reason sm.obs is used as the reference is that it yields
    # the same result as if the whole network were observed.)
    if(!is.null(sm.obs)){
      sm.obs.mean <- colMeans.mcmc.list(sm.obs)
      sm <- sweep.mcmc.list(sm, sm.obs.mean, "-")
      sm.obs <- sweep.mcmc.list(sm.obs, sm.obs.mean, "-")
    }
  }

  if(curved){
    sm <- do.call(mcmc.list, lapply(sm, ergm.sample.eta2theta, coef=object$coef, etamap=object$etamap))
    if(!is.null(sm.obs)) sm.obs <- do.call(mcmc.list, lapply(sm.obs, ergm.sample.eta2theta, coef=object$coef, etamap=object$etamap))
  }

  cat("Sample statistics summary:\n")
  print(summary(sm))
  if(!is.null(sm.obs)){
    cat("Constrained sample statistics summary:\n")
    print(summary(sm.obs))
  }
  
  # This can probably be improved.
  if(is.null(sm.obs)){
    cat("\nAre sample statistics significantly different from observed?\n")
    ds <- colMeans.mcmc.list(sm) - if(!center) object$target.stats else 0
    sds <- apply(as.matrix(sm),2,sd)
    ns <- effectiveSize(sm)

    cv <-  cov(as.matrix(sm))
    
    z <- ds/sds*sqrt(ns)
    chi2<-t(sqrt(ns)*ds)%*%robust.inverse(cv)%*%(ds*sqrt(ns))
  }else{
    cat("\nAre unconstrained sample statistics significantly different from constrained?\n")
    ds <- colMeans.mcmc.list(sm) - if(!center) colMeans.mcmc.list(sm.obs) else 0
    sds <- apply(as.matrix(sm),2,sd)
    sds.obs <- apply(as.matrix(sm.obs),2,sd)
    ns <- effectiveSize(sm)
    ns.obs <- effectiveSize(sm.obs)

    cv <-  cov(as.matrix(sm))
    cv.obs <-  cov(as.matrix(sm.obs))

    
    z <- ds/sqrt(sds^2/ns+sds.obs^2/ns.obs)
    chi2<-t(ds)%*%robust.inverse(t(cv/sqrt(ns))/sqrt(ns)+t(cv.obs/sqrt(ns.obs))/sqrt(ns.obs))%*%(ds)
  }
  p.z <- pnorm(abs(z),lower.tail=FALSE)*2
  p.chi2 <- pchisq(chi2,nvar(sm),lower.tail=FALSE)
  
  m <- rbind(c(ds,NA),c(z,chi2),c(p.z,p.chi2))
  rownames(m) <- c("diff.","test stat.","P-val.")
  colnames(m) <- c(varnames(sm),"Overall (Chi^2)")
  print(m)

  # End simulated vs. observed test.
  
  cat("\nSample statistics cross-correlations:\n")
  print(crosscorr(sm))
  if(!is.null(sm.obs)){
    cat("Constrained sample statistics cross-correlations:\n")
    print(crosscorr(sm.obs))
  }

  cat("\nSample statistics auto-correlation:\n")
  for(chain in seq_along(sm)){
    ac<-autocorr(sm[[chain]],0:5)
    ac<-sapply(seq_len(ncol(sm[[chain]])),
               function(i) ac[,i,i])
    colnames(ac)<-varnames(sm)
    cat("Chain", chain, "\n")
    print(ac)
  }
  if(!is.null(sm.obs)){
    cat("Constrained sample statistics auto-correlation:\n")
      for(chain in seq_along(sm.obs)){
        ac<-autocorr(sm.obs[[chain]],0:5)
        ac<-sapply(seq_len(ncol(sm.obs[[chain]])),
                   function(i) ac[,i,i])
        colnames(ac)<-varnames(sm.obs)
        cat("Chain", chain, "\n")
        print(ac)
      }
  }

  cat("\nSample statistics burn-in diagnostic (Geweke):\n")
  sm.gw<-geweke.diag(sm)
  for(i in seq_along(sm.gw)){
    cat("Chain", chain, "\n")
    print(sm.gw[[i]])
    cat("P-values (lower = worse):\n")
    print(2*pnorm(abs(sm.gw[[i]]$z),lower.tail=FALSE))
  }
  if(!is.null(sm.obs)){
    cat("Sample statistics burn-in diagnostic (Geweke):\n")
    sm.obs.gw<-geweke.diag(sm.obs)
    for(i in seq_along(sm.obs.gw)){
      cat("Chain", chain, "\n")
      print(sm.obs.gw[[i]])
      cat("P-values (lower = worse):\n")
      print(2*pnorm(abs(sm.obs.gw[[i]]$z),lower.tail=FALSE))
    }
  }
  
  if(require(latticeExtra)){  
    plot.mcmc.list.ergm(sm,main="Sample statistics",vars.per.page=vars.per.page,...)
    if(!is.null(sm.obs)) plot.mcmc.list.ergm(sm.obs,main="Constrained sample statistics",vars.per.page=vars.per.page,...)
  }else{
    message("Package latticeExtra is not installed. Falling back on coda's default MCMC diagnostic plots.")
    plot(sm,...)
    if(!is.null(sm.obs)) plot(sm.obs,...)
  }

  invisible(list(degeneracy.value=degeneracy.value,
                 degeneracy.type=degeneracy.type))
}

plot.mcmc.list.ergm <- function(x, main=NULL, vars.per.page=3,...){
  dp <- update(densityplot(x, panel=function(...){panel.densityplot(...);panel.abline(v=0)}),xlab=NULL,ylab=NULL)
  tp <- update(xyplot.mcmc.list.ergm(x, panel=function(...){panel.xyplot(...);panel.loess(...);panel.abline(0,0)}),xlab=NULL,ylab=NULL)

  library(latticeExtra)

  pages <- ceiling(nvar(x)/vars.per.page)
  
  reordering <- c(rbind(seq_len(nvar(x)),nvar(x)+seq_len(nvar(x))))
  
  print(update(c(tp,dp)[reordering],layout=c(2,vars.per.page),as.table=TRUE,main=main))
}


# Some utility functions:
colMeans.mcmc.list<-function(x,...) colMeans(as.matrix(x),...)

sweep.mcmc.list<-function(x, STATS, FUN="-", check.margin=TRUE, ...){
  for(chain in seq_along(x)){
    x[[chain]] <- sweep(x[[chain]], 2, STATS, FUN, check.margin, ...)
  }
  x
}

## The following function is a modified version of xyplot.mcmc.list
## from the coda R package. The original code is Copyright (C) 2005
## Deepayan Sarkar <Deepayan.Sarkar@R-project.org>, Douglas Bates
## <Douglas.Bates@R-project.org>.
##
## It is incorporated into the ergm package under the terms of the
## GPL v3 license.
##
## I hope that they'll incorproate the changes at some point in the
## future.
xyplot.mcmc.list.ergm <-
    function(x, data = NULL,
             outer = FALSE, groups = !outer,
             aspect = "xy", layout = c(1, ncol(x[[1]])),
             default.scales = list(y = list(relation = "free")),
             type = 'l',
             start = 1, thin = 1,
             main = attr(x, "title"),
             ylab = "",
             ...)
{
    if (!is.R()) {
      stop("This function is not yet available in S-PLUS")
    }
    if (groups && outer) warning("'groups=TRUE' ignored when 'outer=TRUE'")
    datalist <- lapply(x, function(x) as.data.frame(x))
    data <- do.call("rbind", datalist)
    form <-
        if (outer)
            eval(parse(text = paste(paste(lapply(names(data), as.name),
                       collapse = "+"), "~.index | .run")))
        else
            eval(parse(text = paste(paste(lapply(names(data), as.name),
                       collapse = "+"), "~.index")))
##     form <-
##         if (outer)
##             as.formula(paste(paste(names(data),
##                                    collapse = "+"),
##                              "~ index | .run"))
##         else
##             as.formula(paste(paste(names(data),
##                                    collapse = "+"),
##                              "~ index"))
    data[[".index"]] <- seq(from = start(x), by = thin(x), length = nrow(datalist[[1]])) ## repeated
    .run <- gl(length(datalist), nrow(datalist[[1]]))
    if (groups && !outer)
        xyplot(form, data = data,
               outer = TRUE,
               layout = layout,
               groups = .run,
               default.scales = default.scales,
               type = type,
               main = main,
               ylab = ylab,
               ...)
    else
        xyplot(form, data = data,
               outer = TRUE,
               layout = layout,
               default.scales = default.scales,
               type = type,
               main = main,
               ylab = ylab,
               ...)
}
