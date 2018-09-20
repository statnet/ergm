/****************************************************
 Macros to make life easier when writing C code for change statistics:  */

/* The OUTVAL and INVAL macros give the "other endnode" of edge e, depending
   on whether it is an in-edge or an out-edge.  Presumably the first endnode
   of the edge is already known in this context. */
#define ML_OUTVAL(ll, e) (ll->onwp->outedges[(e)].value)
#define ML_INVAL(ll, e) (ll->onwp->inedges[(e)].value)

#define ML_N_NODES(ll) (ll->onwp->nnodes) /* Total number of nodes in the network */
#define ML_N_DYADS(ll) (DYADCOUNT(ll->onwp))
#define ML_OUT_DEG(ll) (ll->onwp->outdegree) /* Vector of length N_NODES giving current outdegrees */
#define ML_IN_DEG(ll) (ll->onwp->indegree) /* Vector of length N_NODES giving current indegrees */
#define ML_DIRECTED(ll) (ll->onwp->directed_flag) /* 0 if network is undirected, 1 if directed */
#define ML_N_EDGES(ll) (EDGECOUNT(ll->onwp)) /* Total number of edges in the network currently */

/* 0 if network is not bipartite, otherwise number of nodes of the first type (the first node of the second type has Vertex index BIPARTITE+1 */
#define ML_BIPARTITE(ll) (ll->onwp->bipartite)

/* Get the number of tails and the number of heads consistently for both bipartite and unipartite networks. */
#define ML_N_TAILS(ll) (BIPARTITE(ll) ? ML_BIPARTITE(ll) : ML_N_NODES(ll))
#define ML_N_HEADS(ll) (BIPARTITE(ll) ? ML_N_NODES(ll)-ML_BIPARTITE(ll) : ML_N_NODES(ll))

/* Used for internal purposes:  assigning the next in- and out-edge when
   needed */
#define ML_NEXT_INEDGE_NUM(ll) (ll->onwp->next_inedge)
#define ML_NEXT_OUTEDGE_NUM(ll) (ll->onwp->next_outedge)

