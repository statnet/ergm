/*
 *  File ergm/src/plinfo.c
 *  Part of the statnet package, http://statnet.org
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) in
 *    http://statnet.org/attribution
 *
 *  Copyright 2012 the statnet development team
 */
#include "plinfo.h"

void plinfo_wrapper (int *tails, int *heads, int *dnedges,
		     int *dn, int *dflag, int *nterms, char **funnames,
		     char **sonames, double *inputs,  
		     double *responsevec, double *covmat,
		     int *start, int *end)
{
  Network nw;
  Vertex n_nodes = (Vertex) *dn; 
  Edge n_edges = (Edge) *dnedges;
  int directed_flag = *dflag;
  Model *m;
  Vertex bip=0;  /* Assumes bipartite is irrelevant; is this true? */

  GetRNGstate(); /* Necessary for use of R random number generator */
  nw=NetworkInitialize(tails, heads, n_edges, n_nodes, directed_flag, bip, 1, 0, NULL);
  m=ModelInitialize(*funnames, *sonames, &inputs, *nterms);
  
  plinfoInitialize(responsevec, covmat,(Vertex*)start,(Vertex*)end, &nw, m);
  ModelDestroy(m);
  NetworkDestroy(&nw);
  PutRNGstate(); /* Must be called after GetRNGstate before returning to R */
}

void plinfoInitialize (double *responsevec, double *covmat,
                       Vertex *start, Vertex *end,
                       Network *nwp, Model *m)
{
  int l, outflag, inflag = 0;
  Edge offset1 = (Edge) (nwp->nnodes * (nwp->nnodes - 1))/2;
  Edge offset2 = (Edge) offset1 * m->n_stats;
  Vertex i, j, dc;
  ModelTerm *mtp;
  
  dc = 0;
  for(i=1; i < nwp->nnodes; i++){
    for(j = i+1; j <= nwp->nnodes; j++){
      dc++;
      if((dc >= (*start)) & (dc <= (*end))){
        if (nwp->directed_flag){
          *(responsevec + offset1) = inflag = 
          (EdgetreeSearch(i, j, nwp->inedges) != 0);
        }
        *responsevec++ = outflag = (EdgetreeSearch(i, j, nwp->outedges) != 0);
        
        for (mtp=m->termarray; mtp < m->termarray + m->n_terms; mtp++){
          mtp->dstats = covmat;
          (*(mtp->d_func))(1, &i, &j, mtp, nwp);
          /* dstats are returned for TOGGLED edges; for MPLE, we need values
             reflecting a change from 0->1 instead.  Thus, we have to change 
             the sign of dstats if the edge exists. */
          if (outflag) {
            for (l = 0; l < mtp->nstats; l++)
              mtp->dstats[l] = -mtp->dstats[l];
          }
          if (nwp->directed_flag) {
            /* then get the network stats for these stats */
            mtp->dstats = covmat + offset2;
            (*(mtp->d_func))(1, &j, &i, mtp, nwp);
            if (inflag) {
              for (l=0; l<mtp->nstats; l++)
                mtp->dstats[l] = -mtp->dstats[l];
            }
          }
          covmat += mtp->nstats;    
        }
      }
    }
  }
}


