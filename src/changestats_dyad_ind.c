#include "ergm_changestat.h"
#include "ergm_edgetype_set_binary.h"
#define ECHANGE(a) (edgestate ? -(a) : +(a))
#define ECHANGE1 (edgestate ? -1 : +1)
#define SVARIANT(a) a
#include "changestats_dyad_ind.c.template.do_not_include_directly.h"

S_CHANGESTAT_FN(s_edges) {
  CHANGE_STAT[0] = N_EDGES;
}
