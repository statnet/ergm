/*  File src/changestats_test.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2017 Statnet Commons
 */
#ifndef _ERGM_CHANGESTATS_AUXNET_H_
#define _ERGM_CHANGESTATS_AUXNET_H_

#include "ergm_changestat_auxnet.h"
#include "ergm_dyad_hashmap.h"

typedef struct{StoreDyadSet *nwp; double *ref_el;} StoreDyadSetAndRefEL;

I_CHANGESTAT_FN(i__isociomatrix);
U_CHANGESTAT_FN(u__isociomatrix);
F_CHANGESTAT_FN(f__isociomatrix);

I_CHANGESTAT_FN(i__discord_net_DyadSet);
U_CHANGESTAT_FN(u__discord_net_DyadSet);
F_CHANGESTAT_FN(f__discord_net_DyadSet);

I_CHANGESTAT_FN(i__intersect_net_DyadSet);
U_CHANGESTAT_FN(u__intersect_net_DyadSet);
F_CHANGESTAT_FN(f__intersect_net_DyadSet);

I_CHANGESTAT_FN(i__intersect_net_toggles_in_list_DyadSet);
U_CHANGESTAT_FN(u__intersect_net_toggles_in_list_DyadSet);
F_CHANGESTAT_FN(f__intersect_net_toggles_in_list_DyadSet);

I_CHANGESTAT_FN(i__union_net_DyadSet);
U_CHANGESTAT_FN(u__union_net_DyadSet);
F_CHANGESTAT_FN(f__union_net_DyadSet);

#endif
