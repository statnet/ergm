/*  This is a collection of functions used to calculate diagnostic 
    statistics for dynamic networks. */

#include "infection.h"

/* These #defines are not really necessary but may make the code a bit
   easier to read.  They come at a price (the possibility of hard-to-track
   errors).   */
#define DMATRIX(a,b) (dmatrix[(a)+(offset)*(b)])
#define EDGE(a,b) (edge[(a)+(*nedge)*(b)])
#define CHANGE(a,b) (change[(a)+(*nchange)*(b)])
#define DISSOLVE(a,b) (dissolve[(a)+(*ndissolve)*(b)])
#define OMATRIX(a,b) (omatrix[(a)+(maxo)*(b)])
#define DEGMIXMAT(a,b) (degmixmat[(a)+(*nnodes)*(b)])

void Prevalence (int *nnodes,
      int *nedge, int *edge, int *ntimestep, int *nfem,
      int *ntotal, int *nchange, int *change, int *ndissolve, int *dissolve,
      double *betarate, int *infected, int *prev) {
  Vertex f, m, alter;
  Edge e;
  Vertex *id, *od;
  Edge i, j, ne = *nedge;
  int time;
  int bipartite = *nfem;
  double *heads, *tails;
  double beta=*betarate;
  TreeNode *ie, *oe;
  Network nw, nw0;
  
  /* Set up a statnet Network using existing edgeTree code 
     Must coerce heads and tails to double, */
  heads = (double *) malloc(sizeof(double) * ne);
  tails = (double *) malloc(sizeof(double) * ne);
  for (i=0; i < ne; i++) {
    heads[i] = (double) EDGE(i, 0);
    tails[i] = (double) EDGE(i, 1);
  }
  nw0 = NetworkInitialize(heads, tails, ne, *nnodes, 0, bipartite);
  free (heads);
  free (tails);

  nw = nw0;
  ie=nw.inedges;
  oe=nw.outedges;
  id=nw.indegree;
  od=nw.outdegree;
  
  /* Step through time one click at a time */
  for (time=j=0; time <= *ntimestep; time++) {
    /* Update the infection vector */
    for (i=0; i < *nfem; i++) {
     /* step through outedges of i  */
     if(infected[i]){
      for(e = EdgetreeMinimum(oe, i+1);
	(alter = oe[e].value) != 0;
	e = EdgetreeSuccessor(oe, e)){
	     if(!infected[alter]){
		     if(unif_rand() < beta/od[i]){infected[alter]=1;}
	     }
         }
     }
    }
    for (; i < *nnodes; i++) {
     /* step through outedges of i  */
     if(infected[i]){
      for(e = EdgetreeMinimum(ie, i+1);
	(alter = ie[e].value) != 0;
	e = EdgetreeSuccessor(ie, e)){
	     if(!infected[alter]){
		     if(unif_rand() < beta/id[i]){infected[alter]=1;}
	     }
         }
     }
    }
    /* Toggle the edges at this timestep */
    if (time < *ntimestep) {
      for(; CHANGE(j,0) == time; j++) {
        ToggleEdge(CHANGE(j, 1), CHANGE(j, 2), &nw); 
      }
    }
  }
  for (i=0; i < *nnodes; i++) {
	  prev[0]=prev[0]+infected[i];
  }
  /* Free memory used by network object before returning */  
  NetworkDestroy (&nw);
  NetworkDestroy (&nw0);
}
