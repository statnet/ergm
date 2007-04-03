#include "changestats_duration.h"

/*****************
 void d_D_off
 This gives the change in mean off-duration for all non-existant edges
 (It does not yet work.)
*****************/
void d_D_off (int ntoggles, Vertex *heads, Vertex *tails, 
	           ModelTerm *mtp, Network *nwp) 
{
  int i;
  
  *(mtp->dstats) = 0.0;
  for (i=0; i < ntoggles; i++)
  {
    
    
    if (i+1 < ntoggles)
      ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
  }
  i--;                                                          
  while (--i >= 0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_D_dyad
 This gives the change in mean duration for all dyads (where dyads without
 an edge are considered to have duration zero)
*****************/
void d_D_dyad (int ntoggles, Vertex *heads, Vertex *tails, 
	           ModelTerm *mtp, Network *nwp) 
{
  int edgeflag, i;
  Vertex ii, jj;
  double ndyads;
  long int elapsed;

  *(mtp->dstats) = 0.0;
  ndyads = nwp->nnodes * (nwp->nnodes - 1);
  if (!nwp->directed_flag) ndyads /= 2;
  for (i=0; i < ntoggles; i++)
  {
    edgeflag = (EdgetreeSearch(ii=heads[i], jj=tails[i], nwp->outedges) != 0);
    if (edgeflag) {
      elapsed = ElapsedTime(ii,jj,nwp);
      *(mtp->dstats) += (double) (nwp->nedges - 1.0 - elapsed) / ndyads;
    } else {
      /* If we gain an edge, the change is always negative */
      *(mtp->dstats) += (double)(nwp->nedges + 1.0) / ndyads;
    }
    if (i+1 < ntoggles)
      ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i >= 0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_D_edge
 This gives the change in mean duration for existant edges
*****************/
void d_D_edge (int ntoggles, Vertex *heads, Vertex *tails, 
	           ModelTerm *mtp, Network *nwp) 
{
  int edgeflag, i;
  Vertex ii, jj;
  Edge e=nwp->nedges;
  long int elapsed;
  double md = mean_duration(nwp);

  *(mtp->dstats) = 0.0;
  for (i=0; i < ntoggles; i++)
  {
    edgeflag = (EdgetreeSearch(ii=heads[i], jj=tails[i], nwp->outedges) != 0);
    if (edgeflag) {
      elapsed = ElapsedTime(ii,jj,nwp);
      /* Make sure to prevent division by zero */
      *(mtp->dstats) += e==1 ? -(double)elapsed : (md - elapsed)/(e-1.0);      
    } else {
      *(mtp->dstats) += (1.0-md)/(e+1.0);
    }
    if (i+1 < ntoggles)
      ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i >= 0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

double mean_duration(Network *nwp) 
{
  Vertex i, j;  
  Edge k, e=nwp->nedges;
  double sum=0.0;
    
  /* Repeated use of FindithEdge makes for inefficient code at the moment. */
  for (k=1; k <= e; k++) {
    FindithEdge(&i, &j, k, nwp);
    sum += ElapsedTime(i, j, nwp);
  }
  /* Protect against division by zero */
  return e==0 ? 0.0 : sum/(double)e;
}

