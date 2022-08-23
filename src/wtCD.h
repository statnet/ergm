/*  File src/wtCD.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2022 Statnet Commons
 */
#ifndef WTCD_H
#define WTCD_H

#include "wtMCMC.h"

#define CD_UNDOS_ALLOC \
  Vertex *undotail = R_calloc(MHp->ntoggles * INTEGER(CDparams)[0] * INTEGER(CDparams)[1], Vertex); \
  Vertex *undohead = R_calloc(MHp->ntoggles * INTEGER(CDparams)[0] * INTEGER(CDparams)[1], Vertex); \
  double *undoweight = R_calloc(MHp->ntoggles * INTEGER(CDparams)[0] * INTEGER(CDparams)[1], double);

#define CD_UNDOS_PASS undotail, undohead, undoweight

#define CD_UNDOS_RECEIVE Vertex *undotail, Vertex *undohead, double *undoweight

#define CD_PROP_TOGGLE_PROVISIONAL                                      \
  Vertex t=MHp->toggletail[i], h=MHp->togglehead[i];                    \
  double w=MHp->toggleweight[i];                                        \
  undotail[ntoggled]=t;                                                 \
  undohead[ntoggled]=h;                                                 \
  undoweight[ntoggled]=WtGetEdge(MHp->toggletail[i], MHp->togglehead[i], nwp); \
  ntoggled++;                                                           \
  WtSetEdge(t, h, w, nwp);

#define CD_PROP_UNDO_TOGGLE(idvar)                              \
  Vertex t = undotail[idvar], h = undohead[idvar];              \
  double w = undoweight[idvar];                                 \
  /* FIXME: This should be done in one call, but it's very easy \
     to make a fencepost error here. */                         \
  WtSetEdge(t, h, w, nwp);

#define DISPATCH_CD_wrapper WtCD_wrapper
#define DISPATCH_CDSample WtCDSample
#define DISPATCH_CDStep WtCDStep

#include "CD.h.template.do_not_include_directly.h"

#endif
