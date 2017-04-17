#ifndef CHANGESTATS_NETOP_H
#define CHANGESTATS_NETOP_H

#include "edgetree.h"
#include "changestat.h"
#include "storage.h"

typedef struct{Network nw; double *ref_el;} StoreNetAndRefEL;

I_CHANGESTAT_FN(i__isociomatrix);
U_CHANGESTAT_FN(u__isociomatrix);
F_CHANGESTAT_FN(f__isociomatrix);

I_CHANGESTAT_FN(i__discord_net);
U_CHANGESTAT_FN(u__discord_net);
F_CHANGESTAT_FN(f__discord_net);

I_CHANGESTAT_FN(i__intersect_net);
U_CHANGESTAT_FN(u__intersect_net);
F_CHANGESTAT_FN(f__intersect_net);

I_CHANGESTAT_FN(i__intersect_net_toggles_in_list);
U_CHANGESTAT_FN(u__intersect_net_toggles_in_list);
F_CHANGESTAT_FN(f__intersect_net_toggles_in_list);

I_CHANGESTAT_FN(i__union_net);
U_CHANGESTAT_FN(u__union_net);
F_CHANGESTAT_FN(f__union_net);


#endif
