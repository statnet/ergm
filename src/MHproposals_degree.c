/*
 *  File ergm/src/MHproposals_degree.c
 *  Part of the statnet package, http://statnet.org
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) in
 *    http://statnet.org/attribution
 *
 *  Copyright 2012 the statnet development team
 */
#include "MHproposals_degree.h"

/* Shorthand. */
#define Mhead (MHp->togglehead)
#define Mtail (MHp->toggletail)

void MH_CondDegreeSimpleTetrad(MHproposal *MHp, Network *nwp)  {  
  Vertex A1, A2, B1, B2;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=4;    
    return;
  }
  GetRandEdge(&A1, &A2, nwp);
  GetRandEdge(&B1, &B2, nwp);
  while(A1==B1 || A1==B2 || A2==B1 || A2==B2 || EdgetreeSearch(A1, B2, nwp->outedges) || EdgetreeSearch(B1, A2, nwp->outedges)){
    GetRandEdge(&A1, &A2, nwp);
    GetRandEdge(&B1, &B2, nwp);
//Rprintf("A1 %d A2 %d B1 %d B2 %d\n",A1,A2,B1,B2); 
  }
//Rprintf("in A1 %d A2 %d B1 %d B2 %d\n",A1,A2,B1,B2); 
  Mtail[0]=A1; Mhead[0]=A2;
  Mtail[1]=A1; Mhead[1]=B2;
  Mtail[2]=B1; Mhead[2]=B2;
  Mtail[3]=B1; Mhead[3]=A2;
}

/* 
void MH_CondDegreeSimple

   Select three edges A1->A2, B1->B2, C1->C2 at random and rotate them to 
   A1->B2, B1->C2, and C1->A2.

   Note that while all the *1s need to be distinct and all the *2s
   need to be distinct, all dyads to be created must be empty and
   meaningful (i.e. not loops), it's OK if A1=C2, B1=A2, and/or C1=B2.

 */
void MH_CondDegreeSimpleHexad(MHproposal *MHp, Network *nwp)  {  
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

void MH_CondDegreeSimple (MHproposal *MHp, Network *nwp)  {  
  
  if(MHp->ntoggles == 0) { /* Initialize CondDeg by */
	                   /* Choosing Hexad or Tetrad */
    MHp->ntoggles=6;
    return;
  }

  if( unif_rand() > 0.9 &&  nwp->directed_flag ){
    MHp->ntoggles=6;
    MH_CondDegreeSimpleHexad(MHp, nwp);
  }else{
    MHp->ntoggles=4;
    MH_CondDegreeSimpleTetrad(MHp, nwp);
  }
}
