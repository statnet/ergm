/*  File inst/include/ergm_dyad_hashmap_utils.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2026 Statnet Commons
 */
#ifndef _ERGM_DYAD_HASHMAP_UTILS_H_
#define _ERGM_DYAD_HASHMAP_UTILS_H_

#include "ergm_dyad_hashmap.h"
#include "ergm_edgetree.h"

/* Utility function declarations. */
void PrintDyadMapUInt(StoreDyadMapUInt *h);
void PrintStrictDyadMapUInt(StoreStrictDyadMapUInt *h);
void PrintDyadSet(StoreDyadSet *h);
void PrintStrictDyadSet(StoreStrictDyadSet *h);
StoreDyadSet *NetworkToDyadSet(Network *nwp);
StoreStrictDyadSet *NetworkToStrictDyadSet(Network *nwp);

#endif // _ERGM_DYAD_HASHMAP_H_
