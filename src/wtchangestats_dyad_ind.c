/*  File src/wtchangestats_dyad_ind.c in package ergm, part of the Statnet
 *  suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#include "ergm_wtchangestat.h"
#include "ergm_edgetype_set_double.h"
#define ECHANGE(a) ((a) * (weight - edgestate))
#define ECHANGE1 (weight - edgestate)
#define SVARIANT(a) a ## _sum
#include "changestats_dyad_ind.c.template.do_not_include_directly.h"
#undef ECHANGE
#undef ECHANGE1
#undef SVARIANT

#define ECHANGE(a) ((a) * ((weight!=0) - (edgestate!=0)))
#define ECHANGE1 ((weight!=0) - (edgestate!=0))
#define SVARIANT(a) a ## _nonzero
#include "changestats_dyad_ind.c.template.do_not_include_directly.h"
