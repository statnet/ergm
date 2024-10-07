/*  File src/CD.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2024 Statnet Commons
 */
#ifndef CD_H
#define CD_H

#include "MCMC.h"

#define CD_UNDOS_ALLOC \
  Vertex *undotail = R_calloc(MHp->ntoggles * INTEGER(CDparams)[0] * INTEGER(CDparams)[1], Vertex); \
  Vertex *undohead = R_calloc(MHp->ntoggles * INTEGER(CDparams)[0] * INTEGER(CDparams)[1], Vertex);

#define CD_UNDOS_PASS undotail, undohead

#define CD_UNDOS_RECEIVE Vertex *undotail, Vertex *undohead

#define CD_PROP_TOGGLE_PROVISIONAL                                      \
  Vertex t=MHp->toggletail[i], h=MHp->togglehead[i];                    \
  undotail[ntoggled]=t;                                                 \
  undohead[ntoggled]=h;                                                 \
  ntoggled++;                                                           \
  ToggleEdge(t, h, nwp);

#define CD_PROP_UNDO_TOGGLE(idvar)                              \
  Vertex t = undotail[idvar], h = undohead[idvar];              \
  /* FIXME: This should be done in one call, but it's very easy \
     to make a fencepost error here. */                         \
  ToggleEdge(t, h, nwp);

#define DISPATCH_CD_wrapper CD_wrapper
#define DISPATCH_CDSample CDSample
#define DISPATCH_CDStep CDStep

#include "CD.h.template.do_not_include_directly.h"

#endif
