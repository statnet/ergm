/*  File inst/include/ergm_MHproposal_bd.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#ifndef _ERGM_MHPROPOSAL_BD_H_
#define _ERGM_MHPROPOSAL_BD_H_

#include "ergm_MHproposal.h"
#include "ergm_Rutil.h"

typedef struct DegreeBoundstruct {
  int attrcount;
  int fBoundDegByAttr;
  int *attribs;
  int *maxout;
  int *minout;
  int *maxin;
  int *minin;
} DegreeBound;

DegreeBound* DegreeBoundInitializeR(SEXP MHpR, Network *nwp);

DegreeBound* DegreeBoundInitialize(int *attribs, int *maxout, int *maxin,
			   int *minout, int *minin, int condAllDegExact,
			   int attriblength, Network *nwp);

void DegreeBoundDestroy(DegreeBound *bd);

int CheckTogglesValid(DegreeBound *bd, MHProposal *MHp, Network *nwp);
int CheckConstrainedTogglesValid(DegreeBound *bd, MHProposal *MHp, Network *nwp);

#define BD_LOOP(bd, proc) BD_COND_LOOP(bd, {proc}, TRUE, 1)

#define BD_COND_LOOP(bd, proc, cond, tryfactor)				\
  unsigned int trytoggle;						\
  for(trytoggle = 0; trytoggle < MAX_TRIES*tryfactor; trytoggle++){	\
    {proc}								\
    if(cond) break;							\
  }									\
  if(trytoggle>=MAX_TRIES*tryfactor){					\
    MHp->toggletail[0]=MH_FAILED;					\
    MHp->togglehead[0]=MH_UNSUCCESSFUL;					\
  }else	if(!CheckTogglesValid(bd, MHp, nwp)){				\
    MHp->toggletail[0]=MH_FAILED;					\
    MHp->togglehead[0]=MH_CONSTRAINT;                                   \
  }

#endif
