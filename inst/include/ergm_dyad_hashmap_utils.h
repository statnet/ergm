#ifndef _ERGM_DYAD_HASHMAP_UTILS_H_
#define _ERGM_DYAD_HASHMAP_UTILS_H_

#include "ergm_dyad_hashmap.h"
#include "ergm_edgetree.h"

/* Utility function declarations. */
void PrintDyadMapUInt(StoreDyadMapUInt *h);
void PrintDyadSet(StoreDyadSet *h);
StoreDyadSet *NetworkToDyadSet(Network *nwp);

#endif // _ERGM_DYAD_HASHMAP_H_
