#  File R/ergm-defunct.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################
# The last home for functions to removed from ergm.

#' @name ergm-defunct
#' @title Functions that have been removed from this package
#' @usage
#'
#' robust.inverse(...)
#'
#' plot.network.ergm(...)
#'
#' ergm.getterms(...)
#'
#' plot.mcmc.list.ergm(...)
#'
#' plot.ergm(...)
#'
#' summary.statistics(...)
#'
#' ergm.checkargs(...)
#'
#' ergm.checkbipartite(...)
#'
#' ergm.checkdirected(...)
#'
#' summary.gof(...)
#'
#' ergm.getMCMCsample(...)
#'
#' ergm.MHP.table(...)
#'
#' MHproposal(...)
#'
#' MHproposal.character(...)
#'
#' MHproposal.ergm(...)
#'
#' MHproposal.formula(...)
#'
#' ergm.init.methods(...)
#'
#' ergm.ConstraintImplications(...)
#'
#' ergm.mcmcslave(...)
#'
#' ergm.update.formula(...)
#'
#' remove.offset.formula(...)
#'
#' network.update(...)
#'
#' ergm.getmodel(...)
#'
#' ergm.getglobalstats(...)
#'
#' as.edgelist.compressed(...)
#'
#' as.network.uncompressed(...)
#'
#' standardize.network(...)
#'
#' newnw.extract(...)
#'
#' san.ergm(...)
#'
#' is.inCH(...)
#'
#' as.rlebdm.ergm(...)
#'
#' @description Functions that have been removed after a period of deprecation.
#' @param ... Arguments to defunct functions.
#' @details
#' `robust.inverse()`: use `MASS::ginv()`.
#'
#' `plot.network.ergm()`: use `latentnet::plot.ergmm()`.
#'
#' `ergm.getterms()`: use `statnet.common::list_rhs.formula()` and `statnet.common::eval_lhs.formula()`.
#'
#' `plot.mcmc.list.ergm()`: use `ergm_plot.mcmc.list()`.
#'
#' `plot.ergm()`: use `mcmc.diagnostics()`.
#'
#' `summary.statistics()`: use `summary_formula()`.
#'
#' `ergm.checkargs()`: use `check.ErgmTerm()`.
#'
#' `ergm.checkbipartite()`: use `check.ErgmTerm()`.
#'
#' `ergm.checkdirected()`: use `check.ErgmTerm()`.
#'
#' `summary.gof()`: use `print.gof()`.
#'
#' `ergm.getMCMCsample()`: use `ergm_MCMC_sample()`.
#'
#' `ergm.MHP.table()`: use `ergm_proposal_table()`.
#'
#' `MHproposal()`: use `ergm_proposal()`.
#'
#' `MHproposal.character()`: use `ergm_proposal()`.
#'
#' `MHproposal.ergm()`: use `ergm_proposal()`.
#'
#' `MHproposal.formula()`: use `ergm_proposal()`.
#'
#' `ergm.init.methods()`: Initial methods are now specified in `InitErgmReference.*()` functions.
#'
#' `ergm.ConstraintImplications()`: Implications are now specified in the `InitErgmConstraint.*()` functions.
#'
#' `ergm.mcmcslave()`: use `ergm_MCMC_slave()`.
#'
#' `ergm.update.formula()`: use `statnet.common::nonsimp_update.formula()`.
#'
#' `remove.offset.formula()`: use `statnet.common::filter_rhs.formula()`.
#'
#' `network.update()`: use `update.network()`.
#'
#' `ergm.getmodel()`: use `ergm_model()`.
#'
#' `ergm.getglobalstats()`: use `summary.ergm_model()`.
#'
#' `as.edgelist.compressed()`: no longer used
#'
#' `as.network.uncompressed()`: no longer used
#'
#' `standardize.network()`: obviated by improvements to `network` package.
#'
#' `newnw.extract()`: use `ergm_state` "API"
#'
#' `san.ergm()`: removed due to no meaningful use case
#'
#' `is.inCH()`: use `shrink_into_CH()`.
#'
#' `as.rlebdm.ergm()`: no longer used
#'
#' @aliases robust.inverse plot.network.ergm ergm.getterms plot.mcmc.list.ergm plot.ergm summary.statistics ergm.checkargs ergm.checkbipartite ergm.checkdirected summary.gof ergm.getMCMCsample ergm.MHP.table MHproposal MHproposal.character MHproposal.ergm MHproposal.formula ergm.init.methods ergm.ConstraintImplications ergm.mcmcslave ergm.update.formula remove.offset.formula network.update ergm.getmodel ergm.getglobalstats as.edgelist.compressed as.network.uncompressed standardize.network newnw.extract san.ergm is.inCH as.rlebdm.ergm
#'
#' @keywords internal
NULL
