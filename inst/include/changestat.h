/*  File inst/include/changestat.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2013 Statnet Commons
 */
#ifndef CHANGESTAT_H
#define CHANGESTAT_H

#include "edgetree.h"
#include "edgelist.h"

typedef struct ModelTermstruct {
  void (*c_func)(Vertex, Vertex, struct ModelTermstruct*, Network*);
  void (*d_func)(Edge, Vertex*, Vertex*, struct ModelTermstruct*, Network*);
  void (*i_func)(struct ModelTermstruct*, Network*);
  void (*u_func)(Vertex, Vertex, struct ModelTermstruct*, Network*);
  void (*f_func)(struct ModelTermstruct*, Network*);
  void (*s_func)(struct ModelTermstruct*, Network*);
  double *attrib; /* Ptr to vector of covariates (if necessary; generally unused) */
  int nstats;   /* Number of change statistics to be returned */
  double *dstats; /* ptr to change statistics returned */
  int ninputparams; /* Number of input parameters passed to function */
  double *inputparams; /* ptr to input parameters passed */
  double *statcache; /* vector of the same length as dstats */
  void *storage; /* optional space for persistent storage */
} ModelTerm;


/* binomial coefficient function and macro: */
double my_choose(double n, int r);
#define CHOOSE(n,r) ((n)<(r) ? (0) : (my_choose((double)(n),(int)(r)))) 

/* Comparison macro for doubles: */
#define EQUAL(a,b) (fabs((a)-(b))<0.0000001)

/* Macros to test for logical inequality (XOR) and logical equality (XNOR). */
#define XOR(a,b) (((a)==0) != ((b)==0))
#define XNOR(a,b) (((a)==0) == ((b)==0))

/****************************************************
 Macros to make life easier when writing C code for change statistics:  */

/* return number of tail and head node in the directed node pair
   tail -> head of the selected toggle */
#define TAIL(a) (tails[(a)])
#define HEAD(a) (heads[(a)])

/* tell whether a particular edge exists */
#define IS_OUTEDGE(a,b) (EdgetreeSearch((a),(b),nwp->outedges)!=0?1:0)
#define IS_INEDGE(a,b) (EdgetreeSearch((a),(b),nwp->inedges)!=0?1:0)
#define IS_UNDIRECTED_EDGE(a,b) IS_OUTEDGE(MIN(a,b), MAX(a,b))

/* Return the Edge number of the smallest-labelled neighbor of the node 
   labelled "a".  Or, return the Edge number of the next-largest neighbor 
   starting from the pointer "e", which points to a node in an edgetree. 
   Mostly, these are utility macros used by the STEP_THROUGH_OUTEDGES 
   and STEP_THROUGH_INEDGES macros. */
#define MIN_OUTEDGE(a) (EdgetreeMinimum(nwp->outedges, (a)))
#define MIN_INEDGE(a) (EdgetreeMinimum(nwp->inedges, (a)))
#define NEXT_OUTEDGE(e) (EdgetreeSuccessor(nwp->outedges,(e)))
#define NEXT_INEDGE(e) (EdgetreeSuccessor(nwp->inedges,(e)))

/* Return each of the out-neighbors or in-neighbors, one at a time,
   of node a.  At each iteration of the loop, the variable v gives the node 
   number of the corresponding neighbor.  The e variable, which should be
   initialized as type Edge, is merely the looping variable. */
#define STEP_THROUGH_OUTEDGES(a,e,v) for((e)=MIN_OUTEDGE(a);((v)=OUTVAL(e))!=0;(e)=NEXT_OUTEDGE(e))
#define STEP_THROUGH_INEDGES(a,e,v) for((e)=MIN_INEDGE(a);((v)=INVAL(e))!=0;(e)=NEXT_INEDGE(e))

/* The OUTVAL and INVAL macros give the "other endnode" of edge e, depending
   on whether it is an in-edge or an out-edge.  Presumably the first endnode
   of the edge is already known in this context. */
#define OUTVAL(e) (nwp->outedges[(e)].value)
#define INVAL(e) (nwp->inedges[(e)].value)

/* Change the status of the (a,b) edge:  Add it if it's absent, or 
   delete it if it's present. */
#define TOGGLE(a,b) (ToggleEdge((a),(b),nwp));
#define TOGGLE_DISCORD(a,b) (ToggleEdge((a),(b),nwp+1));

#define N_NODES (nwp->nnodes) /* Total number of nodes in the network */
#define N_DYADS (DYADCOUNT(nwp->nnodes,nwp->bipartite,nwp->directed_flag))
#define OUT_DEG (nwp->outdegree) /* Vector of length N_NODES giving current outdegrees */
#define IN_DEG (nwp->indegree) /* Vector of length N_NODES giving current indegrees */
#define DIRECTED (nwp->directed_flag) /* 0 if network is undirected, 1 if directed */
#define N_EDGES (nwp->nedges) /* Total number of edges in the network currently */

/* 0 if network is not bipartite, otherwise number of nodes of the first type (the first node of the second type has Vertex index BIPARTITE+1 */
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

/* Set all changestats to zero at start of function */
#define ZERO_ALL_CHANGESTATS(a) for((a)=0; (a)<N_CHANGE_STATS; (a)++) CHANGE_STAT[(a)]=0.0

/* Cycle through all toggles proposed for the current step, then
   make the current toggle in case of more than one proposed toggle, then
   undo all of the toggles to reset the original network state.  */


/* *** don't forget tail-> head, so these functions now toggle (tails, heads), instead of (heads, tails) */

#define FOR_EACH_TOGGLE(a) for((a)=0; (a)<ntoggles; (a)++)
#define IF_MORE_TO_COME(a) if((a)+1<ntoggles)
#define TOGGLE_IF_MORE_TO_COME(a) IF_MORE_TO_COME(a){TOGGLE(tails[(a)],heads[(a)])}
#define TOGGLE_DISCORD_IF_MORE_TO_COME(a) IF_MORE_TO_COME(a){TOGGLE_DISCORD(tails[(a)],heads[(a)])}
#define UNDO_PREVIOUS(a) (a)--; while(--(a)>=0)
#define UNDO_PREVIOUS_TOGGLES(a)  UNDO_PREVIOUS(a){TOGGLE(tails[(a)],heads[(a)])}
#define UNDO_PREVIOUS_DISCORD_TOGGLES(a) UNDO_PREVIOUS(a){TOGGLE(tails[(a)],heads[(a)]); TOGGLE_DISCORD(tails[(a)],heads[(a)])}

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

/****************************************************/
/* changestat function prototypes */

/* *** don't forget tail -> head, so this prototype now accepts tails first, not heads first */

#define CHANGESTAT_FN(a) void (a) (Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp)

/* NB:  CHANGESTAT_FN is now deprecated (replaced by D_CHANGESTAT_FN) */
#define C_CHANGESTAT_FN(a) void (a) (Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp)
#define D_CHANGESTAT_FN(a) void (a) (Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp)
#define I_CHANGESTAT_FN(a) void (a) (ModelTerm *mtp, Network *nwp)
#define U_CHANGESTAT_FN(a) void (a) (Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp)
#define F_CHANGESTAT_FN(a) void (a) (ModelTerm *mtp, Network *nwp)
#define S_CHANGESTAT_FN(a) void (a) (ModelTerm *mtp, Network *nwp)

/* This macro wraps two calls to an s_??? function with toggles
   between them. */
#define D_FROM_S							\
  {									\
    (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */		\
    memcpy(mtp->statcache,mtp->dstats,N_CHANGE_STATS*sizeof(double));	\
    /* Note: This cannot be abstracted into EXEC_THROUGH_TOGGLES. */	\
    int j;								\
    FOR_EACH_TOGGLE(j) TOGGLE(TAIL(j),HEAD(j));				\
    (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */		\
    for(unsigned int i=0; i<N_CHANGE_STATS; i++)			\
      mtp->dstats[i] -= mtp->statcache[i];				\
    FOR_EACH_TOGGLE(j) TOGGLE(TAIL(j),HEAD(j));				\
  }

/* This macro constructs a function that wraps D_FROM_S. */
#define D_FROM_S_FN(a) D_CHANGESTAT_FN(a) D_FROM_S

/* Not often used */
#define INPUT_ATTRIB (mtp->attrib)

#endif
