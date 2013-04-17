/*  File src/MHproposals_degree.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2013 Statnet Commons
 */
#include "MHproposals_degree.h"
#include "changestat.h"

/* Shorthand. */
#define Mhead (MHp->togglehead)
#define Mtail (MHp->toggletail)

/* 
void MH_CondDegreeTetrad

   Select two edges with no nodes in common, A1-A2 and B1-B2, s.t. A1-B2 and B1-A2 are not edges, and propose to replace the former two by the latter two.
 */
void MH_CondDegreeTetrad(MHproposal *MHp, Network *nwp)  {  
  Vertex A1, A2, B1, B2;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=4;    
    return;
  }

//nstats = N_CHANGE_STATS;
//matchvaltail = INPUT_PARAM[tail-1+2*nstats];
//matchvalhead = INPUT_PARAM[head-1+2*nstats];
  do{
    /*
      The reason the following randomization is needed is that, given
      the a tetrad with (only) A1-A2 and B1-B2 having edges, there are
      two degree-preserving proposals that can be made: to replace the
      edges with A1-B2 and B1-A2 or to replace with A1-B1 and A2-B2.
      
      GetRandEdge always returns tail<head for undirected networks, so
      a sampler that only uses GetRandEdge(&A1, &A2, nwp) misses out
      on potential proposals. This causes it to become trapped in
      bipartite or near-bipartite configurations. Bipartite networks
      are already bipartite, so they are not affected.

      Swapping A1 and A2 half the time allows either of the above
      proposals to be considered.
    */
    if(!nwp->directed_flag && !nwp->bipartite && unif_rand()<0.5) GetRandEdge(&A2, &A1, nwp);
    else GetRandEdge(&A1, &A2, nwp);

    GetRandEdge(&B1, &B2, nwp);
    //Rprintf("A1 %d A2 %d B1 %d B2 %d\n",A1,A2,B1,B2); 
  }while(A1==B1 || A1==B2 || A2==B1 || A2==B2 || 
	 (nwp->directed_flag ? 
	  (IS_OUTEDGE(A1, B2) || IS_OUTEDGE(B1, A2)) : // Directed
	  (IS_UNDIRECTED_EDGE(A1,B2) || IS_UNDIRECTED_EDGE(B1,A2)) // Undirected
	  ));
  //Rprintf("in A1 %d A2 %d B1 %d B2 %d\n",A1,A2,B1,B2); 
  if(nwp->directed_flag){
    Mtail[0]=A1; Mhead[0]=A2;
    Mtail[1]=A1; Mhead[1]=B2;
    Mtail[2]=B1; Mhead[2]=B2;
    Mtail[3]=B1; Mhead[3]=A2;
  }else{
    Mtail[0]=MIN(A1,A2); Mhead[0]=MAX(A1,A2);
    Mtail[1]=MIN(A1,B2); Mhead[1]=MAX(A1,B2);
    Mtail[2]=MIN(B1,B2); Mhead[2]=MAX(B1,B2);
    Mtail[3]=MIN(B1,A2); Mhead[3]=MAX(B1,A2);
  }
}

/* 
void MH_CondDegreeHexad

   Select three edges A1->A2, B1->B2, C1->C2 at random and rotate them to 
   A1->B2, B1->C2, and C1->A2.

   Note that while all the *1s need to be distinct and all the *2s
   need to be distinct and all dyads to be created must be empty and
   meaningful (i.e. not loops), it's OK if A1=C2, B1=A2, and/or
   C1=B2. Indeed, "reversing" the cyclical triangle is the one
   operation that tetradic toggles can't accomplish.

   Note that this must *never* be called for undirected networks.

 */
void MH_CondDegreeHexad(MHproposal *MHp, Network *nwp)  {  
  Vertex A1, A2, B1, B2, C1, C2;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=6;
    return;
  }

  GetRandEdge(&A1, &A2, nwp);

  do{
    GetRandEdge(&B1, &B2, nwp);
  }while(B1==A1 || B2==A1 || B2==A2 || EdgetreeSearch(A1, B2, nwp->outedges));

  do{
    GetRandEdge(&C1, &C2, nwp);
  }while(C1==A1 || C1==B1 || C1==A2 || C2==A2 || C2==B2 || C2==B1 || EdgetreeSearch(B1, C2, nwp->outedges) || EdgetreeSearch(C1, A2, nwp->outedges));

  Mtail[0]=A1; Mhead[0]=A2;
  Mtail[1]=A1; Mhead[1]=B2;
  Mtail[2]=B1; Mhead[2]=B2;
  Mtail[3]=B1; Mhead[3]=C2;
  Mtail[4]=C1; Mhead[4]=C2;
  Mtail[5]=C1; Mhead[5]=A2;
}

void MH_CondDegree(MHproposal *MHp, Network *nwp)  {  
  
  if(MHp->ntoggles == 0) { /* Initialize CondDeg by */
      MHp->ntoggles= nwp->directed_flag ? 6 : 4;
    return;
  }

  if(nwp->directed_flag && unif_rand() > 0.9){ /* Do the tetrad or hexad proposal. Undirected networks don't need the hexad.*/
    MHp->ntoggles=6;
    MH_CondDegreeHexad(MHp, nwp);
  }else{
    MHp->ntoggles=4;
    MH_CondDegreeTetrad(MHp, nwp);
  }
}

void MH_CondOutDegree(MHproposal *MHp, Network *nwp)  {  
  Vertex A1, A2, B2;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=2;
    return;
  }

  do{
    GetRandEdge(&A1, &A2, nwp);
    B2 = 1 + unif_rand() * nwp->nnodes;
  }while(A1==B2 || A2==B2 || EdgetreeSearch(A1, B2, nwp->outedges));

//Rprintf("A1 %d A2 %d B1 %d B2 %d\n",A1,A2,B1,B2); 
//Rprintf("in A1 %d A2 %d B1 %d B2 %d\n",A1,A2,B1,B2); 
  Mtail[0]=A1; Mhead[0]=A2;
  Mtail[1]=A1; Mhead[1]=B2;
}

void MH_CondInDegree(MHproposal *MHp, Network *nwp)  {  
  Vertex A1, A2, B1;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=2;
    return;
  }

  do{
    GetRandEdge(&A1, &A2, nwp);
    B1 = 1 + unif_rand() * nwp->nnodes;
  }while(A1==B1 || A2==B1 || EdgetreeSearch(B1, A2, nwp->outedges));

//Rprintf("A1 %d A2 %d B1 %d B2 %d\n",A1,A2,B1,B2); 
//Rprintf("in A1 %d A2 %d B1 %d B2 %d\n",A1,A2,B1,B2); 
  Mtail[0]=A1; Mhead[0]=A2;
  Mtail[1]=B1; Mhead[1]=A2;
}


void MH_CondB1Degree(MHproposal *MHp, Network *nwp)  {  
  Vertex A1, A2, B2;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=2;
    return;
  }

  do{
    GetRandEdge(&A1, &A2, nwp);
    B2 = 1 + nwp->bipartite + unif_rand() * (nwp->nnodes-nwp->bipartite);
  }while(A2==B2 || EdgetreeSearch(A1, B2, nwp->outedges));

//Rprintf("A1 %d A2 %d B1 %d B2 %d\n",A1,A2,B1,B2); 
//Rprintf("in A1 %d A2 %d B1 %d B2 %d\n",A1,A2,B1,B2); 
  Mtail[0]=A1; Mhead[0]=A2;
  Mtail[1]=A1; Mhead[1]=B2;
}

void MH_CondB2Degree(MHproposal *MHp, Network *nwp)  {  
  Vertex A1, A2, B1;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=2;
    return;
  }

  do{
    GetRandEdge(&A1, &A2, nwp);
    B1 = 1 + unif_rand() * nwp->bipartite;
  }while(A1==B1 || A2==B1 || EdgetreeSearch(B1, A2, nwp->outedges));

//Rprintf("A1 %d A2 %d B1 %d B2 %d\n",A1,A2,B1,B2); 
//Rprintf("in A1 %d A2 %d B1 %d B2 %d\n",A1,A2,B1,B2); 
  Mtail[0]=A1; Mhead[0]=A2;
  Mtail[1]=B1; Mhead[1]=A2;
}
