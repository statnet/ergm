/*  File src/wtchangestat.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2014 Statnet Commons
 */
#ifndef WTCHANGESTAT_H
#define WTCHANGESTAT_H

#include "wtedgetree.h"

typedef struct WtModelTermstruct {
  void (*d_func)(Edge, Vertex*, Vertex*, double*, struct WtModelTermstruct*, WtNetwork*);
  void (*s_func)(struct WtModelTermstruct*, WtNetwork*);
  double *attrib; /* Ptr to vector of covariates (if necessary; generally unused) */
  int nstats;   /* Number of change statistics to be returned */
  double *dstats; /* ptr to change statistics returned */
  int ninputparams; /* Number of input parameters passed to function */
  double *inputparams; /* ptr to input parameters passed */
  double *statcache; /* vector of the same length as dstats */
} WtModelTerm;


/* binomial coefficient macro: */
#define CHOOSE(n,r) ((n)<(r) ? (0) : (my_choose((double)(n),(int)(r)))) 

/* Macros to test for logical inequality (XOR) and logical equality (XNOR). */
#define XOR(a,b) (((a)==0) != ((b)==0))
#define XNOR(a,b) (((a)==0) == ((b)==0))

/* macros that tell whether a particular edge exists */
#define IS_OUTEDGE(a,b) (WtEdgetreeSearch((a),(b),nwp->outedges)!=0?1:0)
#define IS_INEDGE(a,b) (WtEdgetreeSearch((a),(b),nwp->inedges)!=0?1:0)
#define IS_UNDIRECTED_EDGE(a,b) IS_OUTEDGE(MIN(a,b), MAX(a,b))

/* The OUTVAL and INVAL macros give the "other endnode" of edge e, depending
   on whether it is an in-edge or an out-edge.  Presumably the first endnode
   of the edge is already known in this context. */
#define OUTVAL(e) (nwp->outedges[(e)].value)
#define INVAL(e) (nwp->inedges[(e)].value)

/* The OUTWT and INWT macros give the weight of edge e, depending
   on whether it is an in-edge or an out-edge.  Presumably the first endnode
   of the edge is already known in this context. */
#define OUTWT(e) (nwp->outedges[(e)].weight)
#define INWT(e) (nwp->inedges[(e)].weight)


/* Return the Edge number of the smallest-labelled neighbor of the node 
   labelled "a".  Or, return the Edge number of the next-largest neighbor 
   starting from the pointer "e", which points to a node in an edgetree. 
   Mostly, these are utility macros used by the STEP_THROUGH_OUTEDGES 
   and STEP_THROUGH_INEDGES macros. */
#define MIN_OUTEDGE(a) (WtEdgetreeMinimum(nwp->outedges, (a)))
#define MIN_INEDGE(a) (WtEdgetreeMinimum(nwp->inedges, (a)))
#define NEXT_OUTEDGE(e) (WtEdgetreeSuccessor(nwp->outedges,(e)))
#define NEXT_INEDGE(e) (WtEdgetreeSuccessor(nwp->inedges,(e)))

/* Return each of the out-neighbors or in-neighbors, one at a time,
   of node a.  At each iteration of the loop, the variable v gives the node 
   number of the corresponding neighbor.  The e variable, which should be
   initialized as type Edge, is merely the looping variable. */
#define STEP_THROUGH_OUTEDGES(a,e,v) for((e)=MIN_OUTEDGE(a);((v)=OUTVAL(e))!=0;(e)=NEXT_OUTEDGE(e))
#define STEP_THROUGH_INEDGES(a,e,v) for((e)=MIN_INEDGE(a);((v)=INVAL(e))!=0;(e)=NEXT_INEDGE(e))

// These are "declaring" versions of the above, optimized for use in EXEC_TROUGH_*EDGES macros.
#define STEP_THROUGH_OUTEDGES_DECL(a,e,v) for(Edge e=MIN_OUTEDGE(a);OUTVAL(e)!=0;e=NEXT_OUTEDGE(e))
#define STEP_THROUGH_INEDGES_DECL(a,e,v) for(Edge e=MIN_INEDGE(a);INVAL(e)!=0;e=NEXT_INEDGE(e))

/* Instead of stepping through execute "subroutine" for each neighbor
   automatically adapting to undirected networks. w gets the weight of
   the edge in question. */
/* NOTE: For some reason, GCC complains when it encounters multiple
   declaration in a single statement inside the subroutine, so
   double v1, v2;
   might fail, while
   double v1;
   double v2;
   works.*/
#define EXEC_THROUGH_OUTEDGES(a,e,v,w,subroutine) if(DIRECTED){ STEP_THROUGH_OUTEDGES_DECL(a,e,v) {Vertex v=OUTVAL(e); double w=OUTWT(e); subroutine} } else { EXEC_THROUGH_EDGES(a,e,v,w,subroutine) }
#define EXEC_THROUGH_INEDGES(a,e,v,w,subroutine) if(DIRECTED){ STEP_THROUGH_INEDGES_DECL(a,e,v) {Vertex v=INVAL(e); double w=INWT(e); subroutine} } else { EXEC_THROUGH_EDGES(a,e,v,w,subroutine) }
#define EXEC_THROUGH_EDGES(a,e,v,w,subroutine) { STEP_THROUGH_OUTEDGES_DECL(a,e,v) {Vertex v=OUTVAL(e); double w=OUTWT(e); subroutine}  STEP_THROUGH_INEDGES_DECL(a,e,v) {Vertex v=INVAL(e); double w=INWT(e); subroutine} }

/* Non-adaptive versions of the above. (I.e. ForceOUT/INEDGES.) */
#define EXEC_THROUGH_FOUTEDGES(a,e,v,w,subroutine) STEP_THROUGH_OUTEDGES_DECL(a,e,v) {Vertex v=OUTVAL(e); double w=OUTWT(e); subroutine}
#define EXEC_THROUGH_FINEDGES(a,e,v,w,subroutine) STEP_THROUGH_INEDGES_DECL(a,e,v) {Vertex v=INVAL(e); double w=INWT(e); subroutine}


/* Get and set the weight of the (a,b) edge. */
#define GETWT(a,b) (WtGetEdge(a,b,nwp))
#define SETWT(a,b,w) (WtSetEdge(a,b,w,nwp))

#define N_NODES (nwp->nnodes) /* Total number of nodes in the network */
#define N_DYADS (DYADCOUNT(nwp->nnodes,nwp->bipartite,nwp->directed_flag))
#define OUT_DEG (nwp->outdegree) /* Vector of length N_NODES giving current outdegrees */
#define IN_DEG (nwp->indegree) /* Vector of length N_NODES giving current indegrees */
#define DIRECTED (nwp->directed_flag) /* 0 if network is undirected, 1 if directed */
#define N_EDGES (nwp->nedges) /* Total number of edges in the network currently */

/* 0 if network is not bipartite, otherwise number of first node of second type */
#define BIPARTITE (nwp->bipartite)

/* Used for internal purposes:  assigning the next in- and out-edge when
   needed */
#define NEXT_INEDGE_NUM (nwp->next_inedge)
#define NEXT_OUTEDGE_NUM (nwp->next_outedge)

/* Vector of change statistics to be modified by the function*/
#define CHANGE_STAT (mtp->dstats)
/* Number of change statistics required by the current term */
#define N_CHANGE_STATS (mtp->nstats)

/* Vector of values passed via "inputs" from R */
#define INPUT_PARAM (mtp->inputparams)
#define N_INPUT_PARAMS (mtp->ninputparams) /* Number of inputs passed */

#define TOGGLEIND toggleind_var

/* macro to set all changestats to zero at start of function */
#define ZERO_ALL_CHANGESTATS() for(unsigned int TOGGLEIND=0; TOGGLEIND<N_CHANGE_STATS; TOGGLEIND++) CHANGE_STAT[TOGGLEIND]=0.0

/* Cycle through all toggles proposed for the current step, then
   make the current toggle in case of more than one proposed toggle, then
   undo all of the toggles to reset the original network state.  */
#define FOR_EACH_TOGGLE() for(unsigned int TOGGLEIND=0; TOGGLEIND<ntoggles; TOGGLEIND++)
/* The idea here is to essentially swap the contents of the proposed
   weights with the current weights, and then swap them back when
   done. */
#define TAIL (tail_var)
#define HEAD (head_var)
#define NEWWT (newwt_var)
#define OLDWT (oldwt_var)

#define GETOLDTOGGLEINFO() Vertex TAIL=tails[TOGGLEIND], HEAD=heads[TOGGLEIND]; double OLDWT=GETWT(TAIL,HEAD);
#define GETTOGGLEINFO() GETOLDTOGGLEINFO(); double NEWWT=weights[TOGGLEIND];

/* SETWT_WITH_BACKUP(a) must be called _after_ GETTOGGLEINFO! */
#define SETWT_WITH_BACKUP() {SETWT(TAIL,HEAD,NEWWT); weights[TOGGLEIND]=OLDWT;}
#define UNDO_SETWT() {GETOLDTOGGLEINFO(); SETWT(TAIL,HEAD,weights[TOGGLEIND]); weights[TOGGLEIND]=OLDWT;}
#define SETWT_IF_MORE_TO_COME() {if(TOGGLEIND+1<ntoggles) {SETWT_WITH_BACKUP()}}
#define UNDO_PREVIOUS_SETWTS() for(unsigned int TOGGLEIND=0; TOGGLEIND+1<ntoggles; TOGGLEIND++){UNDO_SETWT()}
/* Brings together the above operations:
   For each toggle:
      Get the current edge weight.
      Calculate the change.
      Back up the current edge weight by swapping weight[i] with current edge weight.
   For each toggle:
      Undo the changes by swapping them back. */
#define EXEC_THROUGH_TOGGLES(subroutine){ZERO_ALL_CHANGESTATS();FOR_EACH_TOGGLE(){ GETTOGGLEINFO(); {subroutine}; SETWT_IF_MORE_TO_COME();}; UNDO_PREVIOUS_SETWTS();}

#define SAMEDYAD(a1,b1,a2,b2) (DIRECTED? a1==a2 && b1==b2 : MIN(a1,b1)==MIN(a2,b2) && MAX(a1,b1)==MAX(a2,b2))

#define GETOLDWT(a,b) (SAMEDYAD(TAIL,HEAD,a,b)?OLDWT:GETWT(a,b))
#define GETNEWWT(a,b) (SAMEDYAD(TAIL,HEAD,a,b)?NEWWT:GETWT(a,b))
#define GETNEWWTOLD(a,b,old) (SAMEDYAD(TAIL,HEAD,a,b)?NEWWT:(old))


/****************************************************/
/* changestat function prototypes, 
   plus a few supporting function prototypes */
#define WtD_CHANGESTAT_FN(a) void (a) (Edge ntoggles, Vertex *tails, Vertex *heads, double *weights, WtModelTerm *mtp, WtNetwork *nwp)
#define WtS_CHANGESTAT_FN(a) void (a) (WtModelTerm *mtp, WtNetwork *nwp)

/* This macro wraps two calls to an s_??? function with toggles
   between them. */
#define D_FROM_S							\
  {									\
    (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */		\
    memcpy(mtp->statcache,mtp->dstats,N_CHANGE_STATS*sizeof(double));	\
    /* Note: This cannot be abstracted into EXEC_THROUGH_TOGGLES. */	\
    FOR_EACH_TOGGLE() { GETTOGGLEINFO(); SETWT_WITH_BACKUP(); }		\
    (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */		\
    for(unsigned int i=0; i<N_CHANGE_STATS; i++)			\
      mtp->dstats[i] -= mtp->statcache[i];				\
    FOR_EACH_TOGGLE() { UNDO_SETWT(); }					\
  }

/* This macro constructs a function that wraps D_FROM_S. */
#define WtD_FROM_S_FN(a) WtD_CHANGESTAT_FN(a) D_FROM_S

/* Not often used */
#define INPUT_ATTRIB (mtp->attrib)

#endif
