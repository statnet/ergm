#ifndef MHproposal_H
#define MHproposal_H

#include "edgetree.h"
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
   MH proposal; and a pointer to an array of Network structures, which 
   contain the network(s).  
   
   Each MH proposal function should check to see whether ntoggles==0
   upon being called.  If so, the MH proposal function should change
   the value of ntoggles to be the largest possible number of toggles
   required, so that this amount of memory can be allocated.
*/


/* *** don't forget tail-> head */

typedef struct MHproposalstruct {
  void (*func)(struct MHproposalstruct*, DegreeBound*, Network*);
  Edge ntoggles;
  Vertex *toggletail;
  Vertex *togglehead;
  double ratio;
  int status;
  double *inputs; /* may be used if needed, ignored if not. */
  /* int multiplicity; Is this needed? I removed all references to
       'multiplicity' everywhere */
} MHproposal;


void MH_init(MHproposal *MH, 
	     char *MHproposaltype, char *MHproposalpackage, 
	     int fVerbose,
	     Network *nwp, DegreeBound *bd);

void MH_free(MHproposal *MH);

int CheckTogglesValid(MHproposal *MHp, DegreeBound *bd, Network *nwp);
int CheckConstrainedTogglesValid(MHproposal *MHp, DegreeBound *bd, Network *nwp);
#endif 



