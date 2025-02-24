#define EDGETYPE(FUN) Wt ## FUN

#include "ergm_edgetree_inline_template.do_not_include_directly.h"

#undef EDGETYPE

/*****************
 int WtGetEdge

Get weighted edge value. Return 0 if edge does not exist.
*****************/
static inline double WtGetEdge (Vertex tail, Vertex head, WtNetwork *nwp)
{
  ENSURE_TH_ORDER;

  Edge oe=WtEdgetreeSearch(tail,head,nwp->outedges);
  if(oe) return nwp->outedges[oe].weight;
  else return 0;
}
