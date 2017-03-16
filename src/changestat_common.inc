/*  File src/changestat.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2013 Statnet Commons
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

/* Set all changestats to zero at start of function: takes arbitrary arguments, for backwards compatibility. */
#define ZERO_ALL_CHANGESTATS(...) memset(CHANGE_STAT, 0, N_CHANGE_STATS*sizeof(double))

/* Storage utilities */

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

/* Not often used */
#define INPUT_ATTRIB (mtp->attrib)

