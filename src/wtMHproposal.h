#ifndef WTMHPROPOSAL_H
#define WTMHPROPOSAL_H

#include "wtedgetree.h"
#include "R_ext/Rdynload.h"

#define NO_EDGE       0x00 /*these four used in realocateWithReplacement */
#define OLD_EDGE      0x01 
#define NEW_EDGE      0x02 
#define CAN_IGNORE    (OLD_EDGE | NEW_EDGE)  

/* Maximum tries (up to an MH-specific constant). */
#define MAX_TRIES 5000

/* MH_* proposal failed codes. */
/* Heads: */
#define MH_FAILED 0
/* Tails: */
#define MH_UNRECOVERABLE 0
#define MH_IMPOSSIBLE 1
#define MH_UNSUCCESSFUL 2


/* Macros to test for logical inequality (XOR) and logical equality (XNOR). */
#define XOR(a,b) (((a)==0) != ((b)==0))
#define XNOR(a,b) (((a)==0) == ((b)==0))

/*  Notes on MHproposal type:
   An MH proposal function must take two arguments:  a pointer to an 
   MHproposal structure, which holds all the information regarding the
   MH proposal; and a pointer to an array of WtNetwork structures, which 
   contain the network(s).  
   
   Each MH proposal function should check to see whether ntoggles==0
   upon being called.  If so, the MH proposal function should change
   the value of ntoggles to be the largest possible number of toggles
   required, so that this amount of memory can be allocated.
*/
typedef struct WtMHproposalstruct {
  void (*func)(struct WtMHproposalstruct*, WtNetwork*);
  Edge ntoggles;
  Vertex *toggletail;
  Vertex *togglehead;
  double *toggleweight;
  double logratio;
  int status;
  double *inputs; /* may be used if needed, ignored if not. */
  /* int multiplicity; Is this needed? I removed all references to
       'multiplicity' everywhere */
} WtMHproposal;


void WtMH_init(WtMHproposal *MH, 
	     char *MHproposaltype, char *MHproposalpackage, 
	       double *inputs,
	     int fVerbose,
	     WtNetwork *nwp);

void WtMH_free(WtMHproposal *MH);

#endif 



