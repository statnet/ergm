/*  File inst/include/ergm_MHproposal.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2020 Statnet Commons
 */
#ifndef _ERGM_MHPROPOSAL_H_
#define _ERGM_MHPROPOSAL_H_

#include "ergm_edgetree.h"
#include "R_ext/Rdynload.h"

typedef struct DegreeBoundstruct {
  int attrcount;
  int fBoundDegByAttr;
  int *attribs;
  int *maxout;
  int *minout;
  int *maxin;
  int *minin;
} DegreeBound;

DegreeBound* DegreeBoundInitialize(int *attribs, int *maxout, int *maxin, 
			   int *minout, int *minin, int condAllDegExact, 
			   int attriblength, Network *nwp);

void DegreeBoundDestroy(DegreeBound *bd);


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

/*  Notes on MHProposal type:
   An MH proposal function must take two arguments:  a pointer to an 
   MHProposal structure, which holds all the information regarding the
   MH proposal; and a pointer to an array of Network structures, which 
   contain the network(s).  
   
   Each MH proposal function should check to see whether ntoggles==0
   upon being called.  If so, the MH proposal function should change
   the value of ntoggles to be the largest possible number of toggles
   required, so that this amount of memory can be allocated.
*/


/* *** don't forget tail-> head */

typedef struct MHProposalstruct {
  void (*func)(struct MHProposalstruct*, Network*);
  Edge ntoggles;
  Vertex *toggletail;
  Vertex *togglehead;
  double logratio;
  int status;
  DegreeBound *bd;
  Network **discord;
  double *inputs; /* may be used if needed, ignored if not. */
} MHProposal;


MHProposal *MHProposalInitialize(
	     char *MHProposaltype, char *MHProposalpackage, 
	     double *inputs,
	     int fVerbose,
	     Network *nwp, 
	     int *attribs, int *maxout, int *maxin, 
	     int *minout, int *minin, int condAllDegExact, 
	     int attriblength);

void MHProposalDestroy(MHProposal *MHp);

int CheckTogglesValid(MHProposal *MHp, Network *nwp);
int CheckConstrainedTogglesValid(MHProposal *MHp, Network *nwp);

#define BD_LOOP(proc) BD_COND_LOOP({proc}, TRUE, 1)

#define BD_COND_LOOP(proc, cond, tryfactor)				\
  unsigned int trytoggle;						\
  for(trytoggle = 0; trytoggle < MAX_TRIES*tryfactor; trytoggle++){	\
    {proc}								\
    if(cond) break;							\
  }									\
  if(trytoggle>=MAX_TRIES*tryfactor){					\
    MHp->toggletail[0]=MH_FAILED;					\
    MHp->togglehead[0]=MH_UNSUCCESSFUL;					\
  }else	if(!CheckTogglesValid(MHp,nwp)){				\
    MHp->toggletail[0]=MH_FAILED;					\
    MHp->togglehead[0]=MH_CONSTRAINT;                                   \
  }									

/* Helper macros */
#define MH_INPUTS MHp->inputs

#define Mtail (MHp->toggletail)
#define Mhead (MHp->togglehead)

#define MH_P_FN(a) void (a) (MHProposal *MHp, Network *nwp)

/* Implementation of TNT log ratio for the three cases.
   The parameters are as follows. Let 
   D = number of dyads
   P = probability of drawing from the set of edges
   
   Then, the arguments are:

   E = number of edges
   Q = 1-P
   DP = D*P
   DO = D*P/(1-P)
 */

/* Thanks to Robert Goudie for pointing out an error in an earlier
   version of this sampler when proposing to go from E==0 to E==1 or
   vice versa.  Note that this happens extremely rarely unless the
   network is small or the parameter values lead to extremely sparse
   networks. */

// Select edge.
#define TNT_LR_E(E, Q, DP, DO) (log(((E)==1 ? 1.0/((DP) + (Q)) : (E) / ((DO) + (E)))))
// Select dyad, get edge.
#define TNT_LR_DE(E, Q, DP, DO) (log(((E)==1 ? 1.0/((DP) + (Q)) : (E) / ((DO) + (E)))))
// Select dyad, get nonedge.
#define TNT_LR_DN(E, Q, DP, DO) (log(((E)==0 ? (DP) + (Q) : 1.0 + (DO)/((E) + 1))))

#endif 



