#  File R/ergm-defunct.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
# The last home for functions to removed from ergm.

#' @name ergm-defunct
#' @usage sociality(object, ...)
#' @title Functions that have been removed from this package
#' @description Functions that have been removed after a period of deprecation.
#' @param ... Arguments to defunct functions.
#' @keywords internal

NULL

# The following were defunct-ed on 2018-04-07.

#' @rdname ergm-defunct
robust.inverse <- function (...) .Defunct("MASS::ginv")

# The following were defunct-ed on 2019-03-07.
#' @rdname ergm-defunct
plot.network.ergm <- function(...) .Defunct("latentnet::plot.ergmm()")

#' @rdname ergm-defunct
ergm.getterms<-function(...) .Defunct("statnet.common::list_rhs.formula() and statnet.common::eval_lhs.formula()")

#' @rdname ergm-defunct
plot.mcmc.list.ergm <- function(...) .Defunct("ergm_plot.mcmc.list()")

#' @rdname ergm-defunct
plot.ergm <- function(...) .Defunct("mcmc.diagnostics(x,...)")

#' @rdname ergm-defunct
summary.statistics <- function(...) .Defunct("summary_formula()")

#' @rdname ergm-defunct
ergm.checkargs <- function(...) .Defunct("check.ErgmTerm")

#' @rdname ergm-defunct
ergm.checkbipartite <- function(...) .Defunct("check.ErgmTerm")

#' @rdname ergm-defunct
ergm.checkdirected <- function(...) .Defunct("check.ErgmTerm")

#' @rdname ergm-defunct
summary.gof <- function(...) .Defunct("print.gof")

#' @rdname ergm-defunct
ergm.getMCMCsample <- function(...) .Defunct("ergm_MCMC_sample")

#' @rdname ergm-defunct
ergm.MHP.table <- function(...) .Defunct("ergm_proposal_table()")

#' @rdname ergm-defunct
MHproposal <- function(...) .Defunct("ergm_proposal()")

#' @rdname ergm-defunct
MHproposal.character <- function(...) .Defunct("ergm_proposal()")

#' @rdname ergm-defunct
MHproposal.ergm <- function(...) .Defunct("ergm_proposal()")

#' @rdname ergm-defunct
MHproposal.formula <- function(...) .Defunct("ergm_proposal()")

#' @rdname ergm-defunct
ergm.init.methods <- function(...) .Defunct(msg="Function ergm.init.methods() has been deprecated in favor of specifying init_methods in InitErgmReference.*() functions, and has no effect.")

#' @rdname ergm-defunct
ergm.ConstraintImplications <- function(...) .Defunct(msg="Function ergm.ConstraintImplications() has been deprecated in favor of specifying the implications in the InitErgmConstraint.*() functions, and has no effect.")

#' @rdname ergm-defunct
ergm.mcmcslave <- function(...) .Defunct("ergm_MCMC_slave")

# The following were defunct-ed on 2019-08-21.

#' @rdname ergm-defunct
ergm.update.formula <- function(...) .Defunct("statnet.common::nonsimp_update.formula")

#' @rdname ergm-defunct
remove.offset.formula <- function(...) .Defunct("statnet.common::filter_rhs.formula")

#' @rdname ergm-defunct
network.update<-function(...) .Defunct("update.network")

#' @rdname ergm-defunct
ergm.getmodel <- function(...) .Defunct("ergm_model")

#' @rdname ergm-defunct
ergm.getglobalstats <- function(...) .Defunct("summary.ergm_model")

#' @rdname ergm-defunct
as.edgelist.compressed<-function(...) .Defunct(msg="No longer used.")

#' @rdname ergm-defunct
as.network.uncompressed<-function(...) .Defunct(msg="No longer used.")

# The following were defunct-ed 2020-09-25.

#' @rdname ergm-defunct
standardize.network <- function(...) .Defunct(msg=paste0("Obviated by improvements to ", sQuote("network"), " package."))

#' @rdname ergm-defunct
newnw.extract<-function(...) .Defunct('ergm_state "API"')

# The following were defunct-ed 2022-06-21.

#' @rdname ergm-defunct
san.ergm <- function(...) .Defunct(msg="Removed due to no meaningful use case.")

# The following were defunct-ed 2022-09-10.
#' @rdname ergm-defunct
is.inCH <- function(...) .Defunct("shrink_into_CH")
