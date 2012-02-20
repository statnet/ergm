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
      int *nedge, int *edge, int *ntimestep, int *nfem, int *nseeds,
      int *ntotal, int *nchange, int *change, int *ndissolve, int *dissolve,
      int *randomseeds, double *betarate, int *infected, int *totinfected,
      int *nsim, int *prev) {
  Vertex alter=0;
  Edge e;
  Vertex *id, *od;
  Edge i, j, ne = *nedge;
  int k, time, ndyads, rane;
  int bipartite = *nfem;
  int *sinfected;
  double beta=*betarate;
  Network nw;
  
  sinfected = (int *) malloc(sizeof(int) * (*nnodes));
//    Rprintf("initial bipartite %d edges %d tails[i] %f heads[i] %f\n", bipartite,ne,
//		           tails[i-1],heads[i-1]);

  for (k=0; k < *nsim; k++) {
    /* R's serialization of matrixes is column-major, so this works: */
   nw = NetworkInitialize(edge, edge+*nedge, ne,
                          *nnodes, 0, bipartite, 0, 0, NULL);
   id=nw.indegree;
   od=nw.outdegree;
   ndyads = bipartite*(*nnodes-bipartite);
   if(*randomseeds){
      // Sample numdissolved edges without replacement
      ndyads = *nnodes;
      for (i = 0; i < ndyads; i++){sinfected[i] = i;}
      for (i = 0; i < (*nseeds); i++) {
        rane = ndyads * unif_rand();
        infected[i] = sinfected[rane] + 1;
        sinfected[rane] = sinfected[--ndyads];
      }
      for (i = 0; i < *nnodes; i++){sinfected[i] = 0;}
//  Rprintf("\n", infected[i]);
      for (i=0; i < (*nseeds); i++) {
//  Rprintf(" %d", infected[i]);
	sinfected[infected[i]] = 1;
      }
   }else{
    for (i=0; i < *nnodes; i++) {
      sinfected[i] = infected[i];
    }
   }
  
   /* Step through time one click at a time */
   for (time=j=0; time <= *ntimestep; time++) {
//Rprintf("time %d \n",time);
    /* Update the infection vector */
    for (i=0; i < *nfem; i++) {
     /* step through outedges of i  */
     if(sinfected[i]){
      for(e = EdgetreeMinimum(nw.outedges, i+1);
	(alter = nw.outedges[e].value) != 0;
	e = EdgetreeSuccessor(nw.outedges, e)){
	     if(!sinfected[alter-1]){
 	       if(unif_rand() < beta){
//  Rprintf("f time %d i %d sinfected %d alter %d sinfected %d beta %f\n",time,i,sinfected[i],alter,sinfected[alter],beta);
		       sinfected[alter-1]=1;}
//       if(unif_rand() < beta/od[i+1]){sinfected[alter-1]=1;}
 	     }
          }
      }
//    Rprintf("time %d i %d sinfected %d beta %f\n",time,i,sinfected[i]);
    }
    for (; i < *nnodes; i++) {
//    Rprintf("males time %d i %d sinfected %d\n",time,i,sinfected[i]);
    /* step through outedges of i  */
    if(sinfected[i]){
     for(e = EdgetreeMinimum(nw.inedges, i+1);
       (alter = nw.inedges[e].value) != 0;
       e = EdgetreeSuccessor(nw.inedges, e)){
          if(!sinfected[alter-1]){
//          if(unif_rand() < beta/id[i]){sinfected[alter-1]=1;}
            if(unif_rand() < beta){
//    Rprintf("m time %d i %d sinfected %d alter %d sinfected %d beta %f\n",time,i,sinfected[i],alter,sinfected[alter],beta);
		    sinfected[alter-1]=1;}
          }
        }
    }
   }
//Rprintf("time %d i %d sinfected %d alter %d sinfected %d beta %f\n",i,sinfected[i],alter,sinfected[alter],beta);
    /* Toggle the edges at this timestep */
    if (time < *ntimestep) {
     for(; CHANGE(j,0) == time; j++) {
        ToggleEdge(CHANGE(j, 1), CHANGE(j, 2), &nw); 
     }
     // End time step toggle
    }
   // End of time step 
   }
   NetworkDestroy (&nw);
   // Next k
   for (i=0; i < *nnodes; i++) {
     prev[k]=prev[k]+sinfected[i];
   }
   for (i=0; i < *nnodes; i++) {
     totinfected[i] = totinfected[i] + sinfected[i];
   }
//   Rprintf("k %d edges %d prev %d \n",k,nw.nedges,prev[k]);
// if (k < *nsim) {
//  NetworkDestroy (&nw);
//  nw = NetworkInitialize(tails, heads, ne, maxedges, *nnodes, 0, bipartite);
//  id=nw.indegree;
//  od=nw.outdegree;
// }
  }
  for (i=0; i < *nnodes; i++) {
      infected[i] = sinfected[i];
  }
  /* Free memory used by network object before returning */  
  free (sinfected);
//NetworkDestroy (&nw);
}
void PrevalenceWithBernoulliOption(int *nnodes,
      int *nedge, int *edge, int *ntimestep, int *nfem,
      int *ntotal, int *nchange, int *change, int *ndissolve, int *dissolve,
      int *bernoulli, double *betarate, int *infected, int *nsim, int *prev) {
  Vertex alter=0;
  Edge e;
  Vertex *id, *od;
  Edge i, j, ne = *nedge, nwedge;
  int k, time, ndyads, rane;
  int bipartite = *nfem;
  Vertex *btails, *bheads;
  int *sinfected, *bsort;
  double beta=*betarate;
  Network nw, nws;
  
  sinfected = (int *) malloc(sizeof(int) * (*nnodes));
//    Rprintf("initial bipartite %d edges %d tails[i] %f heads[i] %f\n", bipartite,ne,
//		           tails[i-1],heads[i-1]);
  nws = NetworkInitialize(edge, edge+*nedge, ne,
                          *nnodes, 0, bipartite, 0, 0, NULL);
  if(*bernoulli){
    nw = NetworkInitialize(edge, edge+*nedge, ne,
                           *nnodes, 0, bipartite, 0, 0, NULL);
  }else{
    nw = nws;
  }
  id=nw.indegree;
  od=nw.outdegree;
  ndyads = bipartite*(nws.nnodes-bipartite);
  bsort = (int *) malloc(sizeof(int) * ndyads);

  for (k=0; k < *nsim; k++) {
   for (i=0; i < *nnodes; i++) {
     sinfected[i] = infected[i];
   }
  
   /* Step through time one click at a time */
   for (time=j=0; time <= *ntimestep; time++) {
//Rprintf("time %d \n",time);
    /* Update the infection vector */
    for (i=0; i < *nfem; i++) {
     /* step through outedges of i  */
     if(sinfected[i]){
      for(e = EdgetreeMinimum(nw.outedges, i+1);
	(alter = nw.outedges[e].value) != 0;
	e = EdgetreeSuccessor(nw.outedges, e)){
	     if(!sinfected[alter-1]){
//  Rprintf("time %d i %d sinfected %d alter %d sinfected %d beta %f\n",time,i,sinfected[i],alter,sinfected[alter],beta);
 	       if(unif_rand() < beta/od[i+1]){sinfected[alter-1]=1;}
 	     }
          }
      }
//    Rprintf("time %d i %d sinfected %d beta %f\n",time,i,sinfected[i]);
    }
    for (; i < *nnodes; i++) {
//    Rprintf("males time %d i %d sinfected %d\n",time,i,sinfected[i]);
    /* step through outedges of i  */
    if(sinfected[i]){
     for(e = EdgetreeMinimum(nw.inedges, i+1);
       (alter = nw.inedges[e].value) != 0;
       e = EdgetreeSuccessor(nw.inedges, e)){
//    Rprintf("time %d i %d sinfected %d alter %d sinfected %d beta %f\n",time,i,sinfected[i],alter,sinfected[alter],beta);
          if(!sinfected[alter-1]){
//    Rprintf("time %d i %d sinfected %d alter %d sinfected %d beta %f\n",time,i,sinfected[i],alter,sinfected[alter],beta);
            if(unif_rand() < beta/id[i]){sinfected[alter-1]=1;}
          }
        }
    }
   }
//Rprintf("time %d i %d sinfected %d alter %d sinfected %d beta %f\n",i,sinfected[i],alter,sinfected[alter],beta);
    /* Toggle the edges at this timestep */
    if (time < *ntimestep) {
     for(; CHANGE(j,0) == time; j++) {
        ToggleEdge(CHANGE(j, 1), CHANGE(j, 2), &nws); 
     }
     if(*bernoulli){
      // Sample numdissolved edges without replacement
      ndyads = bipartite*(nws.nnodes-bipartite);
      btails = (Vertex *) malloc(sizeof(Vertex) * nws.nedges);
      bheads = (Vertex *) malloc(sizeof(Vertex) * nws.nedges);
      for (i = 0; i < ndyads; i++){bsort[i] = i;}
      for (i = 0; i < nws.nedges; i++) {
	rane = ndyads * unif_rand();
	bheads[i] = bsort[rane] + 1;
	bsort[rane] = bsort[--ndyads];
      }
      for (i=0; i < nws.nedges; i++) {
//    Rprintf("i %d sort[i] %f ",i, bheads[i]);
	rane = (double)(((Vertex)(bheads[i]/bipartite)));
	btails[i] = bheads[i] - bipartite*rane;
	bheads[i] = rane + bipartite;
//    Rprintf("i %d btails[i] %f bheads[i] %f\n",i, btails[i],bheads[i]);
      }
//    Rprintf("final k %d time %d bipartite %d edges %d btails[i] %f bheads[i] %f\n",k, time, bipartite,nws.nedges,
//		           btails[i-1],bheads[i-1]);
      NetworkDestroy (&nw);
      nwedge=nws.nedges;
      Network nw;
      nw = NetworkInitialize(btails, bheads, nwedge,
                             *nnodes, 0, bipartite, 0, 0, NULL);
//    Rprintf("network reinitalized for Bernoulli bipartite %d edges %d\n", bipartite,nw.nedges);
      id=nw.indegree;
      od=nw.outdegree;
      free (btails);
      free (bheads);
     }
     // End time step toggle
    }
   // End of time step 
   }
   // Next k
   for (i=0; i < *nnodes; i++) {
     prev[k]=prev[k]+sinfected[i];
   }
//   Rprintf("k %d edges %d prev %d \n",k,nw.nedges,prev[k]);
   if (k < *nsim) {
    NetworkDestroy (&nw);
    nw = NetworkInitialize(edge, edge+*nedge, ne,
                           *nnodes, 0, bipartite, 0, 0, NULL);
    id=nw.indegree;
    od=nw.outdegree;
   }
  }
  /* Free memory used by network object before returning */  
  free (sinfected);
  NetworkDestroy (&nw);
}
