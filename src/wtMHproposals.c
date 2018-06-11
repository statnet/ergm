/*  File src/wtMHproposals.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2018 Statnet Commons
 */
#include "wtMHproposals.h"
#include "ergm_rlebdm.h"

/* Shorthand. */
#define Mtail (MHp->toggletail)
#define Mhead (MHp->togglehead)
#define Mweight (MHp->toggleweight)

/*********************
 void MH_StdNormal

 Default MH algorithm for a standard-normal-reference ERGM
*********************/
void MH_StdNormal(WtMHproposal *MHp, WtNetwork *nwp)  {  
  double oldwt;
  
  if(MHp->ntoggles == 0) { // Initialize StdNormal 
    MHp->ntoggles=1;
    return;
  }
  
  GetRandDyad(Mtail, Mhead, nwp);
  
  oldwt = WtGetEdge(Mtail[0],Mhead[0],nwp);

  const double propsd = 0.2; // This ought to be tunable.

  Mweight[0] = rnorm(oldwt, propsd);    
  
  // Symmetric proposal, but depends on the reference measure
  MHp->logratio += -(Mweight[0]*Mweight[0]-oldwt*oldwt)/2;
}

/*********************
 void MH_Unif

 Default MH algorithm for continuous-uniform-reference ERGM
*********************/
void MH_Unif(WtMHproposal *MHp, WtNetwork *nwp)  {  
  double oldwt;
  static int a, b;
  
  if(MHp->ntoggles == 0) { // Initialize Unif 
    MHp->ntoggles=1;
    a = MHp->inputs[0];
    b = MHp->inputs[1];
    return;
  }
  
  GetRandDyad(Mtail, Mhead, nwp);
  
  oldwt = WtGetEdge(Mtail[0],Mhead[0],nwp);

  do{
    Mweight[0] = runif(a,b);
  }while(Mweight[0]==oldwt);

  MHp->logratio += 0; // h(y) is uniform and the proposal is symmetric
}


/*********************
 void MH_UnifNonObserved

 Missing data MH algorithm for continuous-uniform-reference ERGM.
*********************/
void MH_UnifNonObserved(WtMHproposal *MHp, WtNetwork *nwp)  {  
  static Edge nmissing;
  static int a, b;
  
  if(MHp->ntoggles == 0) { // Initialize Unif 
    MHp->ntoggles=1;
    a = MHp->inputs[0];
    b = MHp->inputs[1];
    nmissing = MHp->inputs[2];
    return;
  }

  if(nmissing==0){
    *Mtail = MH_FAILED;
    *Mhead = MH_IMPOSSIBLE;
    return;
  }


  // Note that missing edgelist is indexed from 0 but the first
  // element of MHp->inputs is the number of missing edges.
  Edge rane = 1 + unif_rand() * nmissing;

  Mtail[0]=MHp->inputs[rane+2];
  Mhead[1]=MHp->inputs[nmissing+rane+2];

  double oldwt = WtGetEdge(Mtail[0],Mhead[0],nwp);

  do{
    Mweight[0] = runif(a,b);
  }while(Mweight[0]==oldwt);

  MHp->logratio += 0; // h(y) is uniform and the proposal is symmetric
}


/*********************
 void MH_DiscUnif

 Default MH algorithm for discrete-uniform-reference ERGM
*********************/
void MH_DiscUnif(WtMHproposal *MHp, WtNetwork *nwp)  {  
  double oldwt;
  static int a, b;
  
  if(MHp->ntoggles == 0) { // Initialize DiscUnif 
    MHp->ntoggles=1;
    a = MHp->inputs[0];
    b = MHp->inputs[1];
    return;
  }
  
  GetRandDyad(Mtail, Mhead, nwp);
  
  oldwt = WtGetEdge(Mtail[0],Mhead[0],nwp);

  do{
    Mweight[0] = floor(runif(a,b+1));
  }while(Mweight[0]==oldwt);

  MHp->logratio += 0; // h(y) is uniform and the proposal is symmetric
}


/*********************
 void MH_DiscUnifNonObserved

 Missing data MH algorithm for discrete-uniform-reference ERGM.
*********************/
void MH_DiscUnifNonObserved(WtMHproposal *MHp, WtNetwork *nwp)  {  
  static Edge nmissing;
  static int a, b;
  
  if(MHp->ntoggles == 0) { // Initialize DiscUnif 
    MHp->ntoggles=1;
    a = MHp->inputs[0];
    b = MHp->inputs[1];
    nmissing = MHp->inputs[2];
    return;
  }

  if(nmissing==0){
    *Mtail = MH_FAILED;
    *Mhead = MH_IMPOSSIBLE;
    return;
  }


  // Note that missing edgelist is indexed from 0 but the first
  // element of MHp->inputs is the number of missing edges.
  Edge rane = 1 + unif_rand() * nmissing;

  Mtail[0]=MHp->inputs[rane+2];
  Mhead[1]=MHp->inputs[nmissing+rane+2];

  double oldwt = WtGetEdge(Mtail[0],Mhead[0],nwp);

  do{
    Mweight[0] = floor(runif(a,b+1));
  }while(Mweight[0]==oldwt);

  MHp->logratio += 0; // h(y) is uniform and the proposal is symmetric
}

/********************
   void MH_DistRLE
   Propose ONLY edges on an RLE-compressed list
***********************/
void MH_DistRLE(WtMHproposal *MHp, WtNetwork *nwp) 
{  
  static RLEBDM1D r;
  static double *inputs;

  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    inputs = MHp->inputs;
    r = unpack_RLEBDM1D(&inputs, nwp->nnodes);
    return;
  }
  
  if(r.ndyads==0){ /* No dyads to toggle. */
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_IMPOSSIBLE;
    return;
  }

  GetRandRLEBDM1D_RS(Mtail, Mhead, &r);
  double oldwt = WtGetEdge(Mtail[0],Mhead[0],nwp);

  do{
    switch((unsigned int) *inputs){
    case 0:
      Mweight[0] = runif(inputs[1],inputs[2]);
      break;
    case 1:
      Mweight[0] = floor(runif(inputs[1],inputs[2]+1));
      break;
    case 2:
      Mweight[0] = rnorm(inputs[1],inputs[2]);
      break;
    case 3:
      Mweight[0] = rpois(inputs[1]);
      break;
    case 4:
      Mweight[0] = rbinom(inputs[1],inputs[2]);
      break;
    }
  }while(Mweight[0]==oldwt);
  
  // MHp->logratio += 0; // h(y) is uniform and the proposal is symmetric
}

