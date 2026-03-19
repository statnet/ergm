#ifndef _WTMHPROPOSALS_H_
#define _WTMHPROPOSALS_H_

#include "ergm_wtMHproposal.h"
#include "ergm_dyadgen.h"
#include "ergm_wtMHproposal_changestat.h"

typedef struct{DyadGen *gen; WtModel* m;} StoreDyadGenAndWtModel;

#define INIT_DYADGEN_AND_WTMODEL(el)                    \
  ALLOC_STORAGE(1, StoreDyadGenAndWtModel, storage);    \
  storage->gen = DyadGenInitializeR(MHp->R, nwp, el);   \
  GET_CHANGESTATS_MODEL(storage->m);

#define DESTROY_DYADGEN_AND_WTMODEL             \
  GET_STORAGE(StoreDyadGenAndWtModel, storage); \
  DyadGenDestroy(storage->gen);

#endif // _WTMHPROPOSALS_H_
