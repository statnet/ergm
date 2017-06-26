/****************************************************
 Macros to make life easier when writing C code for change statistics:  */

/* The OUTVAL and INVAL macros give the "other endnode" of edge e, depending
   on whether it is an in-edge or an out-edge.  Presumably the first endnode
   of the edge is already known in this context. */
#define MN_OUTVAL(sn, e) (sn->onwp->outedges[(e)].value)
#define MN_INVAL(sn, e) (sn->onwp->inedges[(e)].value)

#define MN_N_NODES(sn) (sn->onwp->nnodes) /* Total number of nodes in the network */
#define MN_N_DYADS(sn) (DYADCOUNT(sn->onwp->nnodes,sn->onwp->bipartite,sn->onwp->directed_flag))
#define MN_OUT_DEG(sn) (sn->onwp->outdegree) /* Vector of length N_NODES giving current outdegrees */
#define MN_IN_DEG(sn) (sn->onwp->indegree) /* Vector of length N_NODES giving current indegrees */
#define MN_DIRECTED(sn) (sn->onwp->directed_flag) /* 0 if network is undirected, 1 if directed */
#define MN_N_EDGES(sn) (sn->onwp->nedges) /* Total number of edges in the network currently */

/* 0 if network is not bipartite, otherwise number of nodes of the first type (the first node of the second type has Vertex index BIPARTITE+1 */
#define MN_BIPARTITE(sn) (sn->onwp->bipartite)

/* Get the number of tails and the number of heads consistently for both bipartite and unipartite networks. */
#define MN_N_TAILS(sn) (BIPARTITE(sn) ? BIPARTITE(sn) : N_NODES(sn))
#define MN_N_HEADS(sn) (BIPARTITE(sn) ? N_NODES(sn)-BIPARTITE(sn) : N_NODES(sn))

/* Used for internal purposes:  assigning the next in- and out-edge when
   needed */
#define MN_NEXT_INEDGE_NUM(sn) (sn->onwp->next_inedge)
#define MN_NEXT_OUTEDGE_NUM(sn) (sn->onwp->next_outedge)

/* Vector of change statistics to be modified by the function*/
#define CHANGE_STAT (mtp->dstats)
/* Number of change statistics required by the current term */
#define N_CHANGE_STATS (mtp->nstats)

/* Vector of values passed via "inputs" from R */
#define INPUT_PARAM (mtp->inputparams)
#define N_INPUT_PARAMS (mtp->ninputparams) /* Number of inputs passed */

/* Set all changestats to zero at start of function: takes arbitrary arguments, for backwards compatibility. */
#define ZERO_ALL_CHANGESTATS(...) memset(CHANGE_STAT, 0, N_CHANGE_STATS*sizeof(double))

/* Not often used */
#define INPUT_ATTRIB (mtp->attrib)

