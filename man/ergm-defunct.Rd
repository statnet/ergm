% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ergm-defunct.R
\name{ergm-defunct}
\alias{ergm-defunct}
\alias{robust.inverse}
\alias{plot.network.ergm}
\alias{ergm.getterms}
\alias{plot.mcmc.list.ergm}
\alias{plot.ergm}
\alias{summary.statistics}
\alias{ergm.checkargs}
\alias{ergm.checkbipartite}
\alias{ergm.checkdirected}
\alias{summary.gof}
\alias{ergm.getMCMCsample}
\alias{ergm.MHP.table}
\alias{MHproposal}
\alias{MHproposal.character}
\alias{MHproposal.ergm}
\alias{MHproposal.formula}
\alias{ergm.init.methods}
\alias{ergm.ConstraintImplications}
\alias{ergm.mcmcslave}
\alias{ergm.update.formula}
\alias{remove.offset.formula}
\alias{network.update}
\alias{ergm.getmodel}
\alias{ergm.getglobalstats}
\alias{as.edgelist.compressed}
\alias{as.network.uncompressed}
\alias{standardize.network}
\alias{newnw.extract}
\alias{san.ergm}
\alias{is.inCH}
\alias{as.rlebdm.ergm}
\alias{degreepopularity-ergmTerm}
\alias{odegreepopularity-ergmTerm}
\alias{idegreepopularity-ergmTerm}
\alias{hammingmix-ergmTerm}
\alias{offset.info.formula}
\alias{InitErgmTerm.degreepopularity}
\alias{InitErgmTerm.idegreepopularity}
\alias{InitErgmTerm.odegreepopularity}
\title{Functions that have been removed from this package}
\usage{
robust.inverse(...)

plot.network.ergm(...)

ergm.getterms(...)

plot.mcmc.list.ergm(...)

plot.ergm(...)

summary.statistics(...)

ergm.checkargs(...)

ergm.checkbipartite(...)

ergm.checkdirected(...)

summary.gof(...)

ergm.getMCMCsample(...)

ergm.MHP.table(...)

MHproposal(...)

MHproposal.character(...)

MHproposal.ergm(...)

MHproposal.formula(...)

ergm.init.methods(...)

ergm.ConstraintImplications(...)

ergm.mcmcslave(...)

ergm.update.formula(...)

remove.offset.formula(...)

network.update(...)

ergm.getmodel(...)

ergm.getglobalstats(...)

as.edgelist.compressed(...)

as.network.uncompressed(...)

standardize.network(...)

newnw.extract(...)

san.ergm(...)

is.inCH(...)

as.rlebdm.ergm(...)

offset.info.formula(...)

InitErgmTerm.degreepopularity(...)

InitErgmTerm.idegreepopularity(...)

InitErgmTerm.odegreepopularity(...)
}
\arguments{
\item{...}{Arguments to defunct functions.}
}
\description{
Functions that have been removed after a period of deprecation.
}
\details{
\code{robust.inverse()}: use \code{MASS::ginv()}.

\code{plot.network.ergm()}: use \code{latentnet::plot.ergmm()}.

\code{ergm.getterms()}: use \code{statnet.common::list_rhs.formula()} and \code{statnet.common::eval_lhs.formula()}.

\code{plot.mcmc.list.ergm()}: use \code{ergm_plot.mcmc.list()}.

\code{plot.ergm()}: use \code{mcmc.diagnostics()}.

\code{summary.statistics()}: use \code{summary_formula()}.

\code{ergm.checkargs()}: use \code{check.ErgmTerm()}.

\code{ergm.checkbipartite()}: use \code{check.ErgmTerm()}.

\code{ergm.checkdirected()}: use \code{check.ErgmTerm()}.

\code{summary.gof()}: use \code{print.gof()}.

\code{ergm.getMCMCsample()}: use \code{ergm_MCMC_sample()}.

\code{ergm.MHP.table()}: use \code{ergm_proposal_table()}.

\code{MHproposal()}: use \code{ergm_proposal()}.

\code{MHproposal.character()}: use \code{ergm_proposal()}.

\code{MHproposal.ergm()}: use \code{ergm_proposal()}.

\code{MHproposal.formula()}: use \code{ergm_proposal()}.

\code{ergm.init.methods()}: Initial methods are now specified in \verb{InitErgmReference.*()} functions.

\code{ergm.ConstraintImplications()}: Implications are now specified in the \verb{InitErgmConstraint.*()} functions.

\code{ergm.mcmcslave()}: use \code{ergm_MCMC_slave()}.

\code{ergm.update.formula()}: use \code{statnet.common::nonsimp_update.formula()}.

\code{remove.offset.formula()}: use \code{statnet.common::filter_rhs.formula()}.

\code{network.update()}: use \code{update.network()}.

\code{ergm.getmodel()}: use \code{ergm_model()}.

\code{ergm.getglobalstats()}: use \code{summary.ergm_model()}.

\code{as.edgelist.compressed()}: no longer used

\code{as.network.uncompressed()}: no longer used

\code{standardize.network()}: obviated by improvements to \code{network} package.

\code{newnw.extract()}: use \code{ergm_state} "API"

\code{san.ergm()}: removed due to no meaningful use case

\code{is.inCH()}: use \code{shrink_into_CH()}.

\code{as.rlebdm.ergm()}: no longer used

\code{offset.info.formula()}: no longer used

\code{degreepopularity}, \code{odegreepopularity}, \code{idegreepopularity}: use the corresponding \code{degree1.5} term

\code{hammingmix}: use \code{hamming(...):nodemix(...)} for example
}
\keyword{internal}
