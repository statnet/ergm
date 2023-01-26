/*  File src/wtMCMC.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2023 Statnet Commons
 */
#ifndef _ERGM_WTMCMC_H_
#define _ERGM_WTMCMC_H_

#include "ergm_wtedgetree.h"
#include "ergm_wtchangestat.h"
#include "ergm_wtMHproposal.h"
#include "ergm_wtmodel.h"
#include "ergm_wtstate.h"

#include "ergm_wttype_defs_common.h"

#define DISPATCH_MCMC_wrapper WtMCMC_wrapper
#define DISPATCH_MCMCSample WtMCMCSample
#define DISPATCH_MetropolisHastings WtMetropolisHastings
#define DISPATCH_MCMCPhase12 WtMCMCPhase12
#define DISPATCH_MCMCSamplePhase12 WtMCMCSamplePhase12

#include "MCMC.h.template.do_not_include_directly.h"
#endif
