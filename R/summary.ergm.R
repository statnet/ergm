###############################################################################
# The <summary.ergm> function prints a 'summary of model fit' table and returns
# the components of this table and several others listed below
#
# --PARAMETERS--
#   object     : an ergm object
#
#
# --IGNORED PARAMETERS--
#   ...        : used for flexibility
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
                          total.variation=TRUE,
                          eps=0.0001)
{
# separates out summary and print fns: MSH
  pseudolikelihood <- is.null(object$samplesize) || is.na(object$samplesize)
  independence <- is.dyad.independent(object)
  if(any(is.na(object$coef)) & !is.null(object$mplefit)){
     object$coef[is.na(object$coef)] <-
     object$mplefit$coef[is.na(object$coef)]
  }
  if(is.null(object$hessian) && is.null(object$covar)){
   return()
  }
  nodes<- network.size(object$network)
  dyads<- network.dyadcount(object$network)
  mc.se<- object$mc.se
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
    mc.se <- mc.se*sqrt(mdyads/dyads)
   }
  }
  rownames(asycov) <- names(object$coef)
  colnames(asycov) <- rownames(asycov)
  
  asyse <- diag(asycov)
  if(total.variation & any(!is.na(mc.se))){
   asyse[!is.na(mc.se)] <- asyse[!is.na(mc.se)]+mc.se[!is.na(mc.se)]^2
  }
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

  ans <- list(formula=object$formula,
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
    
  nodes<- network.size(object$network)
  dyads<- network.dyadcount(object$network)
  df <- length(object$coef)

  rdf <- dyads - df
  tval <- object$coef / asyse
  pval <- 2 * pt(q=abs(tval), df=rdf, lower.tail=FALSE)

# Convert to % error
  if(any(!is.na(mc.se))){
#  mc.se[!is.na(mc.se)] <- formatC(round(100*mc.se[!is.na(mc.se)]*mc.se[!is.na(mc.se)]/(asyse[!is.na(mc.se)]*asyse[!is.na(mc.se)])),digits=0,width=5,format="d")
   mc.se[!is.na(mc.se)] <- round(100*mc.se[!is.na(mc.se)]*mc.se[!is.na(mc.se)]/(asyse[!is.na(mc.se)]*asyse[!is.na(mc.se)]))
  }

# values <- format(object$coef,digits=digits)
# names <- names(values)
# names(values) <- NULL
# casyse<-format(asyse, digits=digits)
# cpval<-format(pval, digits=digits)
# cmc.se <- format(mc.se,digits=digits)

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
          asyse[count],mc.se[count],pval[count]))
     count <- count+1
    }

  tempmatrix <- matrix(templist, ncol=4,byrow=TRUE)
  colnames(tempmatrix) <- c("Estimate", "Std. Error", "MCMC %", "p-value")
  rownames(tempmatrix) <- names(object$coef)

  devtext <- "Deviance:"
  if (!independence || object$reference!="Bernoulli") {
    if (pseudolikelihood) {
      devtext <- "Pseudo-deviance:"
      ans$message <- "\nWarning:  The standard errors are based on naive pseudolikelihood and are suspect.\n"
    } 
    else if(any(is.na(mc.se))) {
      ans$message <- "\nWarning:  The standard errors are suspect due to possible poor convergence.\n"
    }
  } else {
    ans$message <- "\nFor this model, the pseudolikelihood is the same as the likelihood.\n"
  }
  llk<-try(logLik(object,...), silent=TRUE)

  if(!inherits(llk,"try-error")){
  
    ans$devtable <- c("",apply(cbind(paste(format(c("   Null", 
                                                    "Residual", ""), width = 8), devtext), 
                                     format(c(object$null.deviance,
                                              -2*llk, 
                                              object$null.deviance+2*llk),
                                            digits = 5), " on",
                                     format(c(dyads, rdf, df),
                                            digits = 5)," degrees of freedom\n"), 
                               1, paste, collapse = " "),"\n")
    
    
    ans$aic <- AIC(llk)
    ans$bic <- BIC(llk)
  }else ans$objname<-deparse(substitute(object))
  
  ans$coefs <- as.data.frame(tempmatrix)
  ans$asycov <- asycov
  ans$asyse <- asyse
  class(ans) <- "summary.ergm"
  ans
}

