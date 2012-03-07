/*
 *  File ergm/src/changestat.h
 *  Part of the statnet package, http://statnet.org
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) in
 *    http://statnet.org/attribution
 *
 *  Copyright 2012 the statnet development team
 */
#ifndef CHANGESTAT_H
#define CHANGESTAT_H

#include "edgetree.h"
#include "edgelist.h"

typedef struct ModelTermstruct {
	void (*d_func)(Edge, Vertex*, Vertex*, struct ModelTermstruct*, Network*);
  	void (*s_func)(struct ModelTermstruct*, Network*);
        void (*t_func)(struct ModelTermstruct*, Network*);
	double *attrib; /* Ptr to vector of covariates (if necessary; generally unused) */
	int nstats;   /* Number of change statistics to be returned */
	double *dstats; /* ptr to change statistics returned */
	int ninputparams; /* Number of input parameters passed to function */
	double *inputparams; /* ptr to input parameters passed */
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
#define N_DYADS (nwp->directed_flag?(nwp->nnodes*(nwp->nnodes-1)):nwp->bipartite?nwp->bipartite*(nwp->nnodes-nwp->bipartite):(nwp->nnodes*(nwp->nnodes-1)/2))
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

/* Set all changestats to zero at start of function */
#define ZERO_ALL_CHANGESTATS(a) for((a)=0; (a)<N_CHANGE_STATS; (a)++) CHANGE_STAT[(a)]=0.0

/* Cycle through all toggles proposed for the current step, then
   make the current toggle in case of more than one proposed toggle, then
   undo all of the toggles to reset the original network state.  */


/* *** don't forget tail-> head, so these functions now toggle (tails, heads), instead of (heads, tails) */

#define FOR_EACH_TOGGLE(a) for((a)=0; (a)<ntoggles; (a)++)
#define TOGGLE_IF_MORE_TO_COME(a) if((a)+1<ntoggles) TOGGLE(tails[(a)],heads[(a)])
#define TOGGLE_DISCORD_IF_MORE_TO_COME(a) if((a)+1<ntoggles) TOGGLE_DISCORD(tails[(a)],heads[(a)])
#define UNDO_PREVIOUS_TOGGLES(a) (a)--; while(--(a)>=0) TOGGLE(tails[(a)],heads[(a)])
#define UNDO_PREVIOUS_DISCORD_TOGGLES(a) (a)--; while(--(a)>=0) {TOGGLE(tails[(a)],heads[(a)]); TOGGLE_DISCORD(tails[(a)],heads[(a)])}

/****************************************************/
/* changestat function prototypes */

/* *** don't forget tail -> head, so this prototype now accepts tails first, not heads first */

#define CHANGESTAT_FN(a) void (a) (Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp)

/* NB:  CHANGESTAT_FN is now deprecated (replaced by D_CHANGESTAT_FN) */
#define D_CHANGESTAT_FN(a) void (a) (Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp)
#define T_CHANGESTAT_FN(a) void (a) (ModelTerm *mtp, Network *nwp)
#define S_CHANGESTAT_FN(a) void (a) (ModelTerm *mtp, Network *nwp)

/* Not often used */
#define INPUT_ATTRIB (mtp->attrib)

#endif
