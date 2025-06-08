/*  File inst/include/ergm_wtMHproposal.h in package ergm, part of the Statnet
 *  suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#ifndef _ERGM_WTMHPROPOSAL_H_
#define _ERGM_WTMHPROPOSAL_H_

#include "ergm_wtedgetree.h"

#include "ergm_edgetype_set_double.h"

#include "inc/ergm_MHproposal.h.template.do_not_include_directly.h"

#include "ergm_edgetype_unset.h"

#define WtMH_I_FN(a) void (a) (WtMHProposal *MHp, WtNetwork *nwp)
#define WtMH_U_FN(a) void (a) (Vertex tail, Vertex head, double weight, WtMHProposal *MHp, WtNetwork *nwp, double edgestate)
#define WtMH_P_FN(a) void (a) (WtMHProposal *MHp, WtNetwork *nwp)
#define WtMH_F_FN(a) void (a) (WtMHProposal *MHp, WtNetwork *nwp)
#define WtMH_X_FN(a) void (a) (unsigned int type, void *data, WtMHProposal* MHp, WtNetwork* nwp)

#endif
