###############################################################################
# The <summary.ergm> function prints a 'summary of model fit' table and returns
# the components of this table and several others listed below
#
# --PARAMETERS--
#   object     : an ergm object
#
#
# --IGNORED PARAMETERS--
#   digits     : significant digits for the coefficients;default=
#                max(3,getOption("digits")-3), but the hard-coded value is 5
#   correlation: whether the correlation matrix of the estimated parameters
#                should be printed (T or F); default=FALSE
#   covariance : whether the covariance matrix of the estimated parameters
#                should be printed (T or F); default=FALSE
#   eps        : the numerical tolerance inputted to the native R function
#                <format.pval>; the code using 'eps' is now commented out;
#                default=.0001
#
# --RETURNED--
#   ans: a "summary.ergm" object as a list containing the following
#      formula         : object$formula
#      randomeffects   : object$re
#      digits          : the 'digits' inputted to <summary.ergm> or the default
#                        value (despite the fact the digits will be 5)
#      correlation     : the 'correlation' passed to <summary.ergm>
#      degeneracy.value: object$degenarcy.value
#      offset          : object$offset
#      drop            : object$drop
#      covariance      : the 'covariance' passed to <summary.ergm>
#      pseudolikelihood: whether pseudoliklelihood was used (T or F)
#      independence    : whether ?? (T or F)
#      iterations      : object$iterations
#      samplesize      : NA if 'pseudolikelihood'=TRUE, object$samplesize otherwise
#      message         : a message regarding the validity of the standard error
#                        estimates
#      aic             : the AIC goodness of fit measure
#      bic             : the BIC goodness of fit measure
#      coefs           : the dataframe of parameter coefficients and their
#                        standard erros and p-values
#      asycov          : the asymptotic covariance matrix
#      asyse           : the asymptotic standard error matrix
#      senderreceivercorrelation: 'randomeffects' if this is a matrix;
#                        otherwise, the correlation between sender and receiver??
#
################################################################################

summary.ergm <- function (object, ..., 
                          digits = max(3, getOption("digits") - 3),
                          correlation=FALSE, covariance=FALSE,
                          eps=0.0001)
{
# separates out summary and print fns: MSH
  pseudolikelihood <- is.null(object$samplesize) || is.na(object$samplesize)
  independence <- !is.null(object$theta1$independent) && all(object$theta1$independent)
  if(any(is.na(object$coef)) & !is.null(object$mplefit)){
     object$coef[is.na(object$coef)] <-
     object$mplefit$coef[is.na(object$coef)]
  }
  if(is.null(object$hessian) && is.null(object$covar)){
   return()
  }
  nodes<- network.size(object$network)
  dyads<- network.dyadcount(object$network)
  if(is.null(object$covar)){
   asycov <- try(robust.inverse(-object$hessian), silent=TRUE)
   if(inherits(asycov,"try-error")){
    asycov <- diag(1/diag(-object$hessian))
   }
  }else{
   asycov <- object$covar
   if(FALSE & !pseudolikelihood & !independence){
    if(is.directed(object$network)){
     mdyads <- nodes * (nodes-1)
    }else{
     mdyads <- nodes * (nodes-1) / 2
    }
    # Adjust for missing dyads
    asycov <- asycov*mdyads/dyads
    object$mc.se <- object$mc.se*sqrt(mdyads/dyads)
   }
  }
  rownames(asycov) <- names(object$coef)
  colnames(asycov) <- rownames(asycov)
  
  asyse <- diag(asycov)
  asyse[asyse<0|is.infinite(object$coef)|object$offset] <- NA
  asyse <- sqrt(asyse)
  if(any(is.na(asyse)&!object$offset) & !is.null(object$mplefit)){
   if(is.null(object$mplefit$covar)){
    if(!is.null(object$mplefit$covar)){
     mpleasycov <- try(robust.inverse(-object$mplefit$hessian), silent=TRUE)
     if(inherits(mpleasycov,"try-error")){
      mpleasycov <- diag(1/diag(-object$mplefit$hessian))
     }
     asyse[is.na(asyse)] <- sqrt(diag(mpleasycov))[is.na(asyse)]
    }
   }else{
    mpleasycov <- object$mplefit$covar
    asyse[is.na(asyse)] <- sqrt(diag(mpleasycov))[is.na(asyse)]
   }
  }
  asyse <- matrix(asyse, ncol=length(asyse))
  colnames(asyse) <- colnames(asycov)
  
#   MSH changed to make clearer
#   original <- format(object$MCMCtheta, digits = digits)
#   original <- format(object$theta.original, digits = digits)

  ans <- list(formula=object$formula, randomeffects=object$re,
              digits=digits, correlation=correlation,
              degeneracy.value = object$degeneracy.value,
              offset = object$offset,
              drop = object$drop,
              covariance=covariance,
              pseudolikelihood=pseudolikelihood,
              independence=independence,
              iterations=object$iterations[1])

  if(ans$pseudolikelihood){
    ans$samplesize <- NA
  }else{
    ans$samplesize <-  object$samplesize
  }
    
  if(!is.null(ans$randomeffects)){ 
   if(!is.matrix(ans$randomeffects)){
    ans$senderreceivercorrelation<-ans$randomeffects
   }else{
    corr <- ans$randomeffects[1,2]/sqrt(ans$randomeffects[1,1]*ans$randomeffects[2,2])
    corr <- max(min(1,corr),-1)
    ans$senderreceivercorrelation<-corr
   }
  }

  nodes<- network.size(object$network)
  dyads<- network.dyadcount(object$network)
  if(!is.null(object$Z.mkl)){
    p <- ncol(object$Z.mkl)
  }else{
    p <- 0
  }
  if(!is.null(object$cluster)){
    df <- length(object$coef) + object$ngroups*(p+2) - 1 # ng-1 + ng *p + ng
  }else{
    df <- length(object$coef) + (nodes - (p + 1)/2) * p
  }
  rdf <- dyads - df
  tval <- object$coef / asyse
  pval <- 2 * pt(q=abs(tval), df=rdf, lower.tail=FALSE)

# values <- format(object$coef,digits=digits)
# names <- names(values)
# names(values) <- NULL
# casyse<-format(asyse, digits=digits)
# cpval<-format(pval, digits=digits)
# cmc.se <- format(object$mc.se,digits=digits)

# cmc.se[object$offset] <- NA
# cpval[object$offset]  <- NA
# casyse[object$offset] <- NA

#  count <- 1
#  templist <- NULL
#  while (count <= length(names))
#   {
#    templist <- append(templist,c(values[count],
#         casyse[count],cpval[count],cmc.se[count]))
#    count <- count+1
#   }
#
#   tempmatrix <- matrix(templist, ncol=4,byrow=TRUE)
#   tempmatrix[,c(1:2,4)] <- format(tempmatrix[,c(1:2,4)], digits=digits, 
#                                   print.gap=2)
#   tempmatrix[,3] <- format.pval(as.numeric(tempmatrix[,3]), digits = 3, eps= 1e-4)
#   colnames(tempmatrix) <- c("estimate","s.e.","p-value","MCMC s.e.")
#   rownames(tempmatrix) <- names

  count <- 1
  templist <- NULL
  while (count <= length(names(object$coef)))
    {
     templist <- append(templist,c(object$coef[count],
          asyse[count],object$mc.se[count],pval[count]))
     count <- count+1
    }

  tempmatrix <- matrix(templist, ncol=4,byrow=TRUE)
  colnames(tempmatrix) <- c("Estimate", "Std. Error", "MCMC s.e.", "p-value")
  rownames(tempmatrix) <- names(object$coef)

  devtext <- "Deviance:"
  if (!independence) {
    if (pseudolikelihood) {
      devtext <- "Pseudo-deviance:"
      ans$message <- "\nWarning:  The standard errors are based on naive pseudolikelihood and are suspect.\n"
    } 
    else if(any(is.na(object$mc.se))) {
      ans$message <- "\nWarning:  The standard errors are suspect due to possible poor convergence.\n"
    }
  } else {
    ans$message <- "\nFor this model, the pseudolikelihood is the same as the likelihood.\n"
  }
  ans$devtable <- c("",apply(cbind(paste(format(c("   Null", 
            "Residual", ""), width = 8), devtext), 
            format(c(object$null.deviance,
                     -2*object$mle.lik, 
                     object$null.deviance+2*object$mle.lik),
                digits = 5), " on",
            format(c(dyads, rdf, df),
                digits = 5)," degrees of freedom\n"), 
            1, paste, collapse = " "),"\n")

  
  ans$aic <- -2*object$mle.lik + 2*df
  ans$bic <- -2*object$mle.lik + log(dyads)*df
  
  ans$coefs <- as.data.frame(tempmatrix)
  ans$asycov <- asycov
  ans$asyse <- asyse
  class(ans) <- "summary.ergm"
  ans
}

