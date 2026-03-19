/*  File src/wtMHproposals.c in package ergm, part of the Statnet suite of
 *  packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2026 Statnet Commons
 */
#include "wtMHproposals.h"
#include "ergm_MHstorage.h"

/*********************
 void MH_DiscUnifTwice

 MH algorithm for discrete-uniform-reference ERGM, twice; mainly used for testing.
*********************/
void MH_DiscUnif2(WtMHProposal *MHp, WtNetwork *nwp)  {  
  double oldwt;
  static int a, b;
  
  if(MHp->ntoggles == 0) { // Initialize DiscUnif 
    MHp->ntoggles=2;
    a = MHp->inputs[0];
    b = MHp->inputs[1];
    return;
  }
  
  GetRandDyad(Mtail, Mhead, nwp);
  
  oldwt = WtGETWT(Mtail[0],Mhead[0]);

  do{
    Mweight[0] = floor(runif(a,b+1));
  }while(Mweight[0]==oldwt);

  do{
    GetRandDyad(Mtail+1, Mhead+1, nwp);
    
    oldwt = WtGETWT(Mtail[1],Mhead[1]);
    
    do{
      Mweight[1] = floor(runif(a,b+1));
    }while(Mweight[1]==oldwt);
  }while(Mtail[0]==Mtail[1] && Mhead[0]==Mhead[1]);
  
  MHp->logratio += 0; // h(y) is uniform and the proposal is symmetric
}

/********************
   void MH_DistRLE
   Propose ONLY edges on an RLE-compressed list
***********************/
WtMH_I_FN(Mi_Dist) {
  INIT_DYADGEN_AND_WTMODEL(FALSE);
  MHp->ntoggles = storage->gen->ndyads!=0 ? 1 : MH_FAILED;
}

WtMH_P_FN(Mp_Dist) {  
  GET_STORAGE(StoreDyadGenAndWtModel, storage);

  DyadGenRandDyad(Mtail, Mhead, storage->gen);
  double oldwt = WtGETWT(Mtail[0], Mhead[0]);

  do{
    switch((unsigned int) *MH_IINPUTS){
    case 0:
      Mweight[0] = runif(MH_DINPUTS[0],MH_DINPUTS[1]);
      break;
    case 1:
      Mweight[0] = floor(runif(MH_DINPUTS[0],MH_DINPUTS[1]+1));
      break;
    case 2:
      Mweight[0] = rnorm(oldwt, MH_DINPUTS[0]);
      // Symmetric proposal, but depends on the reference measure
      MHp->logratio = -(Mweight[0]*Mweight[0]-oldwt*oldwt)/2;
      break;
    case 3:
      Mweight[0] = rpois(MH_DINPUTS[0]);
      break;
    case 4:
      Mweight[0] = rbinom(MH_DINPUTS[0],MH_DINPUTS[1]);
      break;
    }
  }while(Mweight[0] == oldwt);
  
  WtCHECK_CHANGESTATS(storage->m);
}

WtMH_F_FN(Mf_Dist) {
  DESTROY_DYADGEN_AND_WTMODEL;
}
