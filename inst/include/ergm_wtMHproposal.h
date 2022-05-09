/*  File inst/include/ergm_wtMHproposal.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2022 Statnet Commons
 */
#ifndef _ERGM_WTMHPROPOSAL_H_
#define _ERGM_WTMHPROPOSAL_H_

#include "ergm_wtedgetree.h"
#include "R_ext/Rdynload.h"

#define NO_EDGE       0x00 /*these four used in realocateWithReplacement */
#define OLD_EDGE      0x01 
#define NEW_EDGE      0x02 
#define CAN_IGNORE    (OLD_EDGE | NEW_EDGE)  

/* Maximum tries (up to an MH-specific constant). */
#define MAX_TRIES 5000

/* MH_* proposal failed codes. */
/* Tails: */
#define MH_FAILED 0
/* Heads: */
#define MH_UNRECOVERABLE 0
#define MH_IMPOSSIBLE 1
#define MH_UNSUCCESSFUL 2
#define MH_CONSTRAINT 3

/* "Quit" threshold for unsuccessful proposals as a fraction of steps. */
#define MH_QUIT_UNSUCCESSFUL 0.05


/* Macros to test for logical inequality (XOR) and logical equality (XNOR). */
#define XOR(a,b) (((a)==0) != ((b)==0))
#define XNOR(a,b) (((a)==0) == ((b)==0))

/*  Notes on WtMHProposal type:
   An Weighted MH proposal function must take two arguments:  a pointer to an 
   Weighted MHProposal structure, which holds all the information regarding the
   MH proposal; and a pointer to an array of WtNetwork structures, which 
   contain the network(s).  
   
   Each Weighted MH proposal function should check to see whether ntoggles==0
   upon being called.  If so, the Weighted MH proposal function should change
   the value of ntoggles to be the largest possible number of toggles
   required, so that this amount of memory can be allocated.
*/


/* *** don't forget tail-> head */

typedef struct WtMHProposalstruct {
  SEXP R;
  void (*i_func)(struct WtMHProposalstruct*, WtNetwork*);
  void (*p_func)(struct WtMHProposalstruct*, WtNetwork*);
  void (*u_func)(Vertex tail, Vertex head, double weight, struct WtMHProposalstruct*, WtNetwork*);
  void (*f_func)(struct WtMHProposalstruct*, WtNetwork*);
  void (*x_func)(unsigned int type, void *data, struct WtMHProposalstruct*, WtNetwork*);
  Edge ntoggles;
  Vertex *toggletail;
  Vertex *togglehead;
  double *toggleweight;
  double logratio;
  int status;
  double *inputs; /* may be used if needed, ignored if not. */
  int *iinputs; /* may be used if needed, ignored if not. */
  void *storage;
  void **aux_storage;
  unsigned int n_aux;
  unsigned int *aux_slots;
} WtMHProposal;

WtMHProposal *WtMHProposalInitialize(SEXP pR, WtNetwork *nwp, void **aux_storage);

void WtMHProposalDestroy(WtMHProposal *MH, WtNetwork *nwp);

/* Helper macros */
#define MH_DINPUTS MHp->inputs
#define MH_INPUTS MH_DINPUTS
#define MH_IINPUTS MHp->iinputs

#define Mtail (MHp->toggletail)
#define Mhead (MHp->togglehead)
#define Mweight (MHp->toggleweight)

#define WtMH_I_FN(a) void (a) (WtMHProposal *MHp, WtNetwork *nwp)
#define WtMH_U_FN(a) void (a) (Vertex tail, Vertex head, double weight, WtMHProposal *MHp, WtNetwork *nwp)
#define WtMH_P_FN(a) void (a) (WtMHProposal *MHp, WtNetwork *nwp)
#define WtMH_F_FN(a) void (a) (WtMHProposal *MHp, WtNetwork *nwp)
#define WtMH_X_FN(a) void (a) (unsigned int type, void *data, WtMHProposal* MHp, WtNetwork* nwp)

#endif 



