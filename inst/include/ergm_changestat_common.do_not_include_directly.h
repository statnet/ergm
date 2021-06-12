/*  File inst/include/ergm_changestat_common.do_not_include_directly.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2021 Statnet Commons
 */

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

/* The OUTVAL and INVAL macros give the "other endnode" of edge e, depending
   on whether it is an in-edge or an out-edge.  Presumably the first endnode
   of the edge is already known in this context. */
#define OUTVAL(e) (nwp->outedges[(e)].value)
#define INVAL(e) (nwp->inedges[(e)].value)

#define N_NODES (nwp->nnodes) /* Total number of nodes in the network */
#define N_DYADS (DYADCOUNT(nwp))
#define OUT_DEG (nwp->outdegree) /* Vector of length N_NODES giving current outdegrees */
#define IN_DEG (nwp->indegree) /* Vector of length N_NODES giving current indegrees */
#define DIRECTED (nwp->directed_flag) /* 0 if network is undirected, 1 if directed */
#define N_EDGES (EDGECOUNT(nwp)) /* Total number of edges in the network currently */

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
#define INPUT_PARAM DINPUT_PARAM
#define N_INPUT_PARAMS N_DINPUT_PARAMS /* Number of inputs passed */
#define DINPUT_PARAM (mtp->inputparams)
#define N_DINPUT_PARAMS (mtp->ninputparams) /* Number of inputs passed */
#define IINPUT_PARAM (mtp->iinputparams)
#define N_IINPUT_PARAMS (mtp->niinputparams) /* Number of inputs passed */

/* Set all changestats to zero at start of function: takes arbitrary arguments, for backwards compatibility. */
#define ZERO_ALL_CHANGESTATS(...) memset(CHANGE_STAT, 0, N_CHANGE_STATS*sizeof(double))

/* Not often used */
#define DINPUT_ATTRIB (mtp->attrib)
#define INPUT_ATTRIB DINPUT_ATTRIB
#define IINPUT_ATTRIB (mtp->iattrib)
