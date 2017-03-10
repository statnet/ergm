/*  File src/wtchangestat.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2013 Statnet Commons
 */
#ifndef WTCHANGESTAT_H
#define WTCHANGESTAT_H

#include "wtedgetree.h"

typedef struct WtModelTermstruct {
  void (*c_func)(Vertex, Vertex, double, struct WtModelTermstruct*, WtNetwork*);
  void (*d_func)(Edge, Vertex*, Vertex*, double*, struct WtModelTermstruct*, WtNetwork*);
  void (*i_func)(struct WtModelTermstruct*, WtNetwork*);
  void (*u_func)(Vertex, Vertex, double, struct WtModelTermstruct*, WtNetwork*);
  void (*f_func)(struct WtModelTermstruct*, WtNetwork*);
  void (*s_func)(struct WtModelTermstruct*, WtNetwork*);
  double *attrib; /* Ptr to vector of covariates (if necessary; generally unused) */
  int nstats;   /* Number of change statistics to be returned */
  double *dstats; /* ptr to change statistics returned */
  int ninputparams; /* Number of input parameters passed to function */
  double *inputparams; /* ptr to input parameters passed */
  double *statcache; /* vector of the same length as dstats */
  void *storage; /* optional space for persistent storage */
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

/* Exectute through all edges (nonzero values) in the network. */
#define EXEC_THROUGH_NET_EDGES(a,b,e,w,subroutine) for(Vertex a=1; a <= N_NODES; a++)  EXEC_THROUGH_FOUTEDGES(a, e, b, w, {subroutine});

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

/* Get the number of tails and the number of heads consistently for both bipartite and unipartite networks. */
#define N_TAILS (BIPARTITE ? BIPARTITE : N_NODES)
#define N_HEADS (BIPARTITE ? N_NODES-BIPARTITE : N_NODES)

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
#define ZERO_ALL_CHANGESTATS() for(unsigned int TOGGLEIND=0; TOGGLEIND<N_CHANGE_STATS; TOGGLEIND++) memset(CHANGE_STAT, 0, sizeof(double)*N_CHANGE_STATS);

/* Cycle through all toggles proposed for the current step, then
   make the current toggle in case of more than one proposed toggle, then
   undo all of the toggles to reset the original network state.  */
#define FOR_EACH_TOGGLE for(unsigned int TOGGLEIND=0; TOGGLEIND<ntoggles; TOGGLEIND++)
/* The idea here is to essentially swap the contents of the proposed
   weights with the current weights, and then swap them back when
   done. */
#define TAIL (tail_var)
#define HEAD (head_var)
#define NEWWT (newwt_var)
#define OLDWT (oldwt_var)

#define GETOLDTOGGLEINFO() Vertex TAIL=tails[TOGGLEIND], HEAD=heads[TOGGLEIND]; double OLDWT=GETWT(TAIL,HEAD);
#define GETTOGGLEINFO() GETOLDTOGGLEINFO(); double NEWWT=weights[TOGGLEIND];
#define GETNEWTOGGLEINFO() Vertex TAIL=tails[TOGGLEIND], HEAD=heads[TOGGLEIND]; double NEWWT=weights[TOGGLEIND];

/* SETWT_WITH_BACKUP(a) must be called _after_ GETTOGGLEINFO! */
#define SETWT_WITH_BACKUP() {SETWT(TAIL,HEAD,NEWWT); weights[TOGGLEIND]=OLDWT;}
#define UNDO_SETWT() {GETOLDTOGGLEINFO(); SETWT(TAIL,HEAD,weights[TOGGLEIND]); weights[TOGGLEIND]=OLDWT;}
#define IF_MORE_TO_COME if(TOGGLEIND+1<ntoggles)
#define SETWT_IF_MORE_TO_COME() {IF_MORE_TO_COME{SETWT_WITH_BACKUP();}}
#define UNDO_PREVIOUS for(int TOGGLEIND=ntoggles-2; TOGGLEIND>=0; TOGGLEIND--)
#define UNDO_PREVIOUS_SETWTS() {UNDO_PREVIOUS{UNDO_SETWT();}}
/* Brings together the above operations:
   For each toggle:
      Get the current edge weight.
      Calculate the change.
      Back up the current edge weight by swapping weight[i] with current edge weight.
   For each toggle:
      Undo the changes by swapping them back. */
#define EXEC_THROUGH_TOGGLES(subroutine){FOR_EACH_TOGGLE{ GETTOGGLEINFO(); {subroutine}; SETWT_IF_MORE_TO_COME();}; UNDO_PREVIOUS_SETWTS();}

#define SAMEDYAD(a1,b1,a2,b2) (DIRECTED? a1==a2 && b1==b2 : MIN(a1,b1)==MIN(a2,b2) && MAX(a1,b1)==MAX(a2,b2))

#define GETOLDWT(a,b) (SAMEDYAD(TAIL,HEAD,a,b)?OLDWT:GETWT(a,b))
#define GETNEWWT(a,b) (SAMEDYAD(TAIL,HEAD,a,b)?NEWWT:GETWT(a,b))
#define GETNEWWTOLD(a,b,old) (SAMEDYAD(TAIL,HEAD,a,b)?NEWWT:(old))

#define ALLOC_STORAGE(nmemb, stored_type, store_into) stored_type *store_into = (stored_type *) (mtp->storage = calloc(nmemb, sizeof(stored_type)));
#define GET_STORAGE(stored_type, store_into) stored_type *store_into = (stored_type *) mtp->storage;

#define ALLOC_AUX_STORAGE(nmemb, stored_type, store_into) stored_type *store_into = (stored_type *) (nwp->aux_storage[(unsigned int) INPUT_PARAM[0]] = calloc(nmemb, sizeof(stored_type)));
#define GET_AUX_STORAGE(stored_type, store_into) stored_type *store_into = (stored_type *) nwp->aux_storage[(unsigned int) INPUT_PARAM[0]];
#define GET_AUX_STORAGE_NUM(stored_type, store_into, ind) stored_type *store_into = (stored_type *) nwp->aux_storage[(unsigned int) INPUT_PARAM[ind]];

/* Allocate a sociomatrix as auxiliary storage. */
#define ALLOC_AUX_SOCIOMATRIX(stored_type, store_into)			\
  /* Note: the following code first sets up a 2D array indexed from 0, then shifts all pointers by -1 so that sm[t][h] would work for vertex IDs. */ \
  ALLOC_AUX_STORAGE(N_TAILS, stored_type*, store_into);			\
  Dyad sm_size = BIPARTITE? N_TAILS*N_HEADS : DIRECTED ? N_NODES*N_NODES : N_NODES*(N_NODES+1)/2; /* For consistency, and possible future capabilities, include diagonal: */ \
  ALLOC_STORAGE(sm_size, stored_type, data); /* A stored_type* to data. */ \
  Dyad pos = 0;	  /* Start of the next row's data in the data vector. */ \
  for(Vertex t=0; t<N_TAILS; t++){                                      \
  /* First set up the pointer to the right location in the data vector, */ \
  if(BIPARTITE){							\
  store_into[t] = data+pos - N_TAILS; /* This is so that store_into[t][h=BIPARTITE] would address the 0th element of that row. */ \
  pos += N_HEADS;							\
  }else if(DIRECTED){							\
    store_into[t] = data+pos;						\
    pos += N_HEADS;							\
  }else{ /* Undirected. */						\
    store_into[t] = data+pos - t; /* tail <= head, so this is so that store_into[t][h=t] would address the 0th element of that row. */ \
    pos += N_HEADS-t+1; /* Each row has N_NODES - t + 1 elements (including diagonal). */ \
  }									\
  store_into[t]--; /* Now, shift the pointer by -1. */			\
  }									\
									\
  store_into--; /* Shift the pointer array by -1. */			\
  nwp->aux_storage[(unsigned int) INPUT_PARAM[0]] = store_into; /* This is needed to make sure the pointer array itself is updated. */

/* Free a sociomatrix in auxiliary storage. */
#define FREE_AUX_SOCIOMATRIX						\
  unsigned int myslot = (unsigned int) INPUT_PARAM[0];			\
  /* If we hadn't shifted the pointers by -1, this would not have been necessary. */ \
  GET_AUX_STORAGE(void*, sm);						\
  free(sm + 1);								\
  nwp->aux_storage[myslot] = NULL;					\
  /* nwp->storage was not shifted, so it'll be freed automatically. */	


#define INIT_STORAGE(stored_type, store_into, initialization_code)	\
  stored_type *store_into;						\
  if(!mtp->storage){							\
    store_into = (stored_type *) (mtp->storage = malloc(sizeof(stored_type))); \
    {initialization_code};						\
  }else store_into = (stored_type *) mtp->storage;

/****************************************************/
/* changestat function prototypes, 
   plus a few supporting function prototypes */
#define WtC_CHANGESTAT_FN(a) void (a) (Vertex tail, Vertex head, double weight, WtModelTerm *mtp, WtNetwork *nwp)
#define WtD_CHANGESTAT_FN(a) void (a) (Edge ntoggles, Vertex *tails, Vertex *heads, double *weights, WtModelTerm *mtp, WtNetwork *nwp)
#define WtI_CHANGESTAT_FN(a) void (a) (WtModelTerm *mtp, WtNetwork *nwp)
#define WtU_CHANGESTAT_FN(a) void (a) (Vertex tail, Vertex head, double weight, WtModelTerm *mtp, WtNetwork *nwp)
#define WtF_CHANGESTAT_FN(a) void (a) (WtModelTerm *mtp, WtNetwork *nwp)
#define WtS_CHANGESTAT_FN(a) void (a) (WtModelTerm *mtp, WtNetwork *nwp)

/* This macro wraps two calls to an s_??? function with toggles
   between them. */
#define D_FROM_S							\
  {									\
    (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */		\
    memcpy(mtp->statcache,mtp->dstats,N_CHANGE_STATS*sizeof(double));	\
    /* Note: This cannot be abstracted into EXEC_THROUGH_TOGGLES. */	\
    FOR_EACH_TOGGLE{ GETTOGGLEINFO(); SETWT_WITH_BACKUP(); }		\
    (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */		\
    for(unsigned int i=0; i<N_CHANGE_STATS; i++)			\
      mtp->dstats[i] -= mtp->statcache[i];				\
    FOR_EACH_TOGGLE{ UNDO_SETWT(); }					\
  }

/* This macro constructs a function that wraps D_FROM_S. */
#define WtD_FROM_S_FN(a) WtD_CHANGESTAT_FN(a) D_FROM_S

/* Not often used */
#define INPUT_ATTRIB (mtp->attrib)

#endif
