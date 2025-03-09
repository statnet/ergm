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
