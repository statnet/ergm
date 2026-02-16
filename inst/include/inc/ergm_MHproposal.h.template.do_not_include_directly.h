/*  File inst/include/inc/ergm_MHproposal.h.template.do_not_include_directly.h
 *  in package ergm, part of the Statnet suite of packages for network
 *  analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2026 Statnet Commons
 */

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

#ifdef __cplusplus
extern "C" {
#endif

/*  Notes on ETYPE(MHProposal) type:
   An Weighted MH proposal function must take two arguments:  a pointer to an
   Weighted MHProposal structure, which holds all the information regarding the
   MH proposal; and a pointer to an array of ETYPE(Network) structures, which
   contain the network(s).

   Each Weighted MH proposal function should check to see whether ntoggles==0
   upon being called.  If so, the Weighted MH proposal function should change
   the value of ntoggles to be the largest possible number of toggles
   required, so that this amount of memory can be allocated.
*/


/* *** don't forget tail-> head */

typedef struct ETYPE(MHProposalstruct) {
  SEXP R;
  void (*i_func)(struct ETYPE(MHProposalstruct)*, ETYPE(Network)*);
  void (*p_func)(struct ETYPE(MHProposalstruct)*, ETYPE(Network)*);
  void (*u_func)(Vertex, Vertex, IFEWT(EWTTYPE,) struct ETYPE(MHProposalstruct)*, ETYPE(Network)*, EWTTYPE);
  void (*f_func)(struct ETYPE(MHProposalstruct)*, ETYPE(Network)*);
  void (*x_func)(unsigned int, void *, struct ETYPE(MHProposalstruct)*, ETYPE(Network)*);
  Edge ntoggles;
  Vertex *toggletail;
  Vertex *togglehead;
  IFEWT(EWTTYPE *toggleweight;)
  double logratio;
  unsigned int ninputs;
  double *inputs;
  unsigned int niinputs;
  int *iinputs;
  void *storage;
  void **aux_storage;
  unsigned int n_aux;
  unsigned int *aux_slots;
} ETYPE(MHProposal);

ETYPE(MHProposal) *ETYPE(MHProposalInitialize)(SEXP pR, ETYPE(Network) *nwp, void **aux_storage);

void ETYPE(MHProposalDestroy)(ETYPE(MHProposal) *MH, ETYPE(Network) *nwp);

#ifdef __cplusplus
}
#endif

/* Helper macros */
#define MH_N_DINPUTS MHp->ninputs
#define MH_DINPUTS MHp->inputs
#define MH_N_INPUTS MH_N_DINPUTS
#define MH_INPUTS MH_DINPUTS
#define MH_N_IINPUTS MHp->niinputs
#define MH_IINPUTS MHp->iinputs

#define Mtail (MHp->toggletail)
#define Mhead (MHp->togglehead)
#define Mweight (MHp->toggleweight)
