/*  File src/MHproposals.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2023 Statnet Commons
 */
#ifndef MHPROPOSALS_H
#define MHPROPOSALS_H

#include "ergm_MHproposal.h"
#include "ergm_MHproposal_bd.h"
#include "ergm_dyadgen.h"

typedef struct{DyadGen *gen; DegreeBound *bd;} StoreDyadGenAndDegreeBound;

#endif 
