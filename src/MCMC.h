/*  File src/MCMC.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#ifndef _MCMC_H_
#define _MCMC_H_

#include "ergm_constants.h"
#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_MHproposal.h"
#include "ergm_model.h"
#include "ergm_state.h"

#include "ergm_type_defs_common.h"

#define DISPATCH_MCMC_wrapper MCMC_wrapper
#define DISPATCH_MCMCSample MCMCSample
#define DISPATCH_MetropolisHastings MetropolisHastings
#define DISPATCH_MCMCPhase12 MCMCPhase12
#define DISPATCH_MCMCSamplePhase12 MCMCSamplePhase12

#include "MCMC.h.template.do_not_include_directly.h"

#endif
