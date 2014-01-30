%  File man/ergm-internal.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2013 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergm-internal}
\title{Internal ergm Objects}


\alias{InitErgmTerm.conddegmple}
\alias{InitErgmTerm.indices}
\alias{as.edgelist.network}
\alias{ergm.ConstraintImplications}
\alias{ergm.MCMC.packagenames}
\alias{ergm.MHP.table}
\alias{ergm.init.methods}
\alias{which.package.InitFunction}
% MH proposals ought to get a documentation entry eventually.

\alias{InitMHP.CondDegree}
\alias{InitMHP.CondDegreeDist}
\alias{InitMHP.CondInDegreeDist}
\alias{InitMHP.CondOutDegreeDist}
\alias{InitMHP.ConstantEdges}
\alias{InitMHP.HammingConstantEdges}
\alias{InitMHP.HammingTNT}
\alias{InitMHP.TNT}
\alias{InitMHP.TNT10}
\alias{InitMHP.nobetweengroupties}
\alias{InitMHP.randomtoggle}
\alias{InitMHP.randomtoggleNonObserved}
\alias{InitMHP.CondB1Degree}
\alias{InitMHP.CondB2Degree}
\alias{InitMHP.CondInDegree}
\alias{InitMHP.CondOutDegree}
\alias{InitMHP.blockdiag}
\alias{InitMHP.blockdiagNonObserved}
\alias{InitMHP.blockdiagTNT}
\alias{InitWtMHP.DescRank}
\alias{InitWtMHP.DescRankEquivalent}
\alias{InitWtMHP.StdNormal}
\alias{InitWtMHP.StdNormalRank}

\alias{MHproposal}
\alias{MHproposal.MHproposal}
\alias{MHproposal.NULL}
\alias{MHproposal.character}
\alias{MHproposal.ergm}
\alias{MHproposal.formula}

\alias{approx.hotelling.diff.test}
\alias{append.rhs.formula}
\alias{as.edgelist}
\alias{as.edgelist.compressed}
\alias{as.network.uncompressed}
\alias{central.network}
\alias{check.ErgmTerm}
\alias{coef.length.model}
\alias{coef.sublength.model}
\alias{colMeans.mcmc.list}
\alias{compress.data.frame}
\alias{control.ergm.toplevel}
\alias{control.simulate.ergm.toplevel}
\alias{copy.named}
\alias{degreedist}
\alias{degreedistfactor}
\alias{delete.isolates}
\alias{drawcircle}
\alias{drawpie}
\alias{dspartnerdist}
\alias{ergm.Cprepare}
\alias{ergm.Cprepare.el}
\alias{ergm.Cprepare.miss}
\alias{ergm.MCMCse}
\alias{ergm.MCMCse.lognormal}
\alias{ergm.MCMLE}
\alias{ergm.MPIcluster.started}
\alias{ergm.bounddeg}
\alias{ergm.bridge.0.llk}
\alias{ergm.bridge.hammingstart.llk}
\alias{ergm.bridge.preproc}
\alias{ergm.checkargs}
\alias{ergm.checkbipartite}
\alias{ergm.checkconstraints.model}
\alias{ergm.checkdegeneracy}
\alias{ergm.checkdirected}
\alias{ergm.checkextreme.model}
\alias{ergm.coefficient.degeneracy}
\alias{ergm.compute.degeneracy}
\alias{ergm.curved.statsmatrix}
\alias{ergm.degeneracy}
\alias{ergm.design}
\alias{ergm.estimate}
\alias{ergm.eta}
\alias{ergm.etagrad}
\alias{ergm.etagradmult}
\alias{ergm.etamap}
\alias{ergm.geodesicmatrix}
\alias{ergm.geodesicmatrix.edgelist}
\alias{ergm.geodistdist}
\alias{ergm.geodistn}
\alias{ergm.getCluster}
\alias{ergm.getMCMCsample}
\alias{ergm.getglobalstats}
\alias{ergm.getmodel}
\alias{ergm.getnetwork}
\alias{ergm.getterms}
\alias{ergm.independencemodel}
\alias{ergm.initialfit}
\alias{ergm.logisticdeviance}
\alias{ergm.logitreg}
\alias{ergm.mahalanobis}
\alias{ergm.mapl}
\alias{ergm.maple}
\alias{ergm.mcmcslave}
\alias{ergm.mple}
\alias{ergm.nodegeodesics}
\alias{ergm.pairgeodesic}
\alias{ergm.pen.glm}
\alias{ergm.phase12}
\alias{ergm.pl}
\alias{ergm.plinfo}
\alias{ergm.reviseinit}
\alias{ergm.rhs.formula}
\alias{ergm.robmon}
\alias{ergm.sample.eta2theta}
\alias{ergm.sample.tomcmc}
\alias{ergm.stepping}
\alias{ergm.stocapprox}
\alias{ergm.stopCluster}
\alias{ergm.stopCluster.MPIcluster}
\alias{ergm.stopCluster.default}
\alias{ergm.sufftoprob}
\alias{ergm.t.summary}
\alias{ergm.update.formula}
\alias{espartnerdist}
\alias{get.InitErgm.fname}
\alias{get.node.attr}
\alias{init.methods}
\alias{is.ergm}
\alias{is.invertible}
\alias{largest.components}
\alias{llik.fun}
\alias{llik.fun.EF}
\alias{llik.fun.median}
\alias{llik.fun.obs}
\alias{llik.fun.obs.robust}
\alias{llik.fun2}
\alias{llik.fun3}
\alias{llik.grad}
\alias{llik.grad.obs}
\alias{llik.grad2}
\alias{llik.grad3}
\alias{llik.hessian}
\alias{llik.hessian.naive}
\alias{llik.hessian.obs}
\alias{llik.hessian2}
\alias{llik.info3}
\alias{llik.mcmcvar3}
\alias{logistftest}
\alias{logp.rank.ergm}
\alias{mcmc.diagnostics.default}
\alias{midarrow}
\alias{model.matrix.pen.glm}
\alias{model.transform.formula}
\alias{mvmodel}
\alias{mvmodel.default}
\alias{mvmodel.ergm}
\alias{mvmodel.formula}
\alias{mk.conlist}
\alias{myLibLoc}
\alias{newnw.extract}
\alias{nologLik.message}
\alias{nvattr.copy.network}
\alias{ostar2deg}
\alias{paste.and}
\alias{plot.mcmc.list.ergm}
\alias{robust.inverse}
\alias{rspartnerdist}
\alias{set.offset.formula}
\alias{sociality.default}
\alias{sociality.ergm}
\alias{sociality.formula}
\alias{sociality.network}
\alias{sort.data.frame}
\alias{summary.statsmatrix.ergm}
\alias{sweep.mcmc.list}
\alias{term.list.formula}
\alias{twopathdist}
\alias{updatemodel.ErgmTerm}
\alias{vector.namesmatch}
\alias{xyplot.mcmc.list.ergm}

\alias{print.control.list}
\alias{set.control.class}
\alias{check.control.class}

\alias{prune.conlist}
\alias{standardize.network}
\alias{unset.offset.formula}

\alias{is.dyad.independent.NULL}
\alias{is.dyad.independent.conlist}
\alias{is.dyad.independent.ergm.model}


\alias{get.free.dyads}
\alias{get.miss.dyads}
\alias{geweke.diag.ar}
\alias{invert.network}

\alias{mk.edge.to.pos.lasttoggle.f}
\alias{as.edgelist.matrix}

\alias{offset.info.formula}
\alias{remove.offset.formula}


\description{
Internal \code{\link[=ergm-package]{ergm}} functions.}

\details{
Most of these are not to be called by the user (or in some cases are just
waiting for proper documentation to be written :).}

\seealso{ergm-package}
\keyword{internal}
