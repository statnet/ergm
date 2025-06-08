/*  File inst/include/ergm_MHproposal.h in package ergm, part of the Statnet
 *  suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#ifndef _ERGM_MHPROPOSAL_H_
#define _ERGM_MHPROPOSAL_H_

#include "ergm_edgetree.h"

#include "ergm_edgetype_set_binary.h"

#include "inc/ergm_MHproposal.h.template.do_not_include_directly.h"

#include "ergm_edgetype_unset.h"

#define MH_I_FN(a) void a (MHProposal *MHp, Network *nwp)
#define MH_U_FN(a) void a (Vertex tail, Vertex head, MHProposal *MHp, Network *nwp, Rboolean edgestate)
#define MH_P_FN(a) void a (MHProposal *MHp, Network *nwp)
#define MH_F_FN(a) void a (MHProposal *MHp, Network *nwp)
#define MH_X_FN(a) void a (unsigned int type, void *data, MHProposal* MHp, Network* nwp)

/* Implementation of TNT log ratio for the three cases.
   The parameters are as follows. Let 
   D = number of dyads
   P = probability of drawing from the set of edges
   
   Then, the arguments are:

   E = number of edges
   Q = 1-P
   DP = D*P
   DO = D*P/(1-P)
 */

/* Thanks to Robert Goudie for pointing out an error in an earlier
   version of this sampler when proposing to go from E==0 to E==1 or
   vice versa.  Note that this happens extremely rarely unless the
   network is small or the parameter values lead to extremely sparse
   networks. */

// Select edge.
#define TNT_LR_E(E, Q, DP, DO) (log(((E)==1 ? 1.0/((DP) + (Q)) : (E) / ((DO) + (E)))))
// Select dyad, get edge.
#define TNT_LR_DE(E, Q, DP, DO) (log(((E)==1 ? 1.0/((DP) + (Q)) : (E) / ((DO) + (E)))))
// Select dyad, get nonedge.
#define TNT_LR_DN(E, Q, DP, DO) (log(((E)==0 ? (DP) + (Q) : 1.0 + (DO)/((E) + 1))))

#endif 
