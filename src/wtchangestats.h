/*  File src/wtchangestats.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2024 Statnet Commons
 */
#ifndef WTCHANGESTATS_H
#define WTCHANGESTATS_H

#include "ergm_wtedgetree.h"
#include "ergm_wtchangestat.h"
#include "ergm_storage.h"
#include "ergm_Rutil.h"


/********************  Utility macros    ***********/

// 0 = untransformed
// 1 = square root
#define TRANSFORM_DYADVAL(y, transcode) (transcode==0? y : transcode==1? sqrt(y): 0)

#endif
