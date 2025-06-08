#  File R/control.ergm3.R in package ergm, part of the Statnet suite of
#  packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
#' @rdname control.ergm
#'
#' @description
#' `control.ergm3()` is a wrapper that sets the defaults to use
#' algorithms and settings circa \pkg{ergm} 3.11 (the last release
#' before 4.0).
#'
#' @export control.ergm3
control.ergm3<-function(
                        MCMLE.termination=c("Hummel", "confidence", "Hotelling", "precision", "none"),
                        MCMLE.effectiveSize=NULL
                        ){}

.ce.args <- formals(control.ergm)
.ce3.args <- formals(control.ergm3)
.ce.args[names(.ce3.args)] <- .ce3.args
formals(control.ergm3) <- .ce.args
body(control.ergm3) <- body(control.ergm)

rm(.ce.args, .ce3.args)
