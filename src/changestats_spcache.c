#include "ergm_changestat.h"
#include "ergm_wtedgetree.h"
#include "ergm_storage.h"

// Increment the value of the edge by inc.
static inline void WtIncEdge (Vertex tail, Vertex head, double inc, WtNetwork *nwp) 
{
  ENSURE_TH_ORDER;

  if(inc==0) return;
  
  // Find the out-edge
  Edge oe=WtEdgetreeSearch(tail,head,nwp->outedges);
  if(oe){
    // Find the corresponding in-edge.
    Edge ie=WtEdgetreeSearch(head,tail,nwp->inedges);
    double weight = nwp->outedges[oe].weight + inc;
    if(weight==0){
      // If the function is to set the edge value to 0, just delete it.
      WtDeleteEdgeFromTrees(tail,head,nwp);
    }else{
      nwp->inedges[ie].weight=nwp->outedges[oe].weight = weight;
    }
  }else{
    // Otherwise, create a new edge with that weight 0+inc==inc.
    WtAddEdgeToTrees(tail,head,inc,nwp);
  }
}

/* Construct and maintain a directed weighted network whose (i,j)
   value is the number of directed two-paths from i to j. */

I_CHANGESTAT_FN(i__otp_wtnet){
  ALLOC_AUX_STORAGE(1, WtNetwork, wtnwp);
  *wtnwp = WtNetworkInitialize(NULL, NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, 0, 0, NULL);
  EXEC_THROUGH_NET_EDGES_PRE(i, j, e1, { // Since i->j
      EXEC_THROUGH_FOUTEDGES_PRE(j, e2, k, { // and j->k
	  if(i!=k)
	    WtIncEdge(i,k,1,wtnwp); // increment i->k.
	});
    });
}

U_CHANGESTAT_FN(u__otp_wtnet){
  GET_AUX_STORAGE(WtNetwork, wtnwp);
  int echange = (IS_OUTEDGE(tail, head) == 0) ? 1 : -1;

  {
    // Update all t->h->k two-paths.
    EXEC_THROUGH_FOUTEDGES_PRE(head, e, k, {
	if(tail!=k)
	  WtIncEdge(tail,k,echange,wtnwp);
      });
  }
  {
    // Update all k->t->h two-paths.
    EXEC_THROUGH_FINEDGES_PRE(tail, e, k, {
	if(k!=head)
	  WtIncEdge(k,head,echange,wtnwp);
      });
  }
}

F_CHANGESTAT_FN(f__otp_wtnet){
  GET_AUX_STORAGE(WtNetwork, wtnwp);

  WtNetworkDestroy(wtnwp);
}


/* Construct and maintain an undirected weighted network whose (i,j)
   value is the number of outgoing shared partners of i and j. */

I_CHANGESTAT_FN(i__osp_wtnet){
  ALLOC_AUX_STORAGE(1, WtNetwork, wtnwp);
  *wtnwp = WtNetworkInitialize(NULL, NULL, NULL, 0, N_NODES, FALSE, BIPARTITE, 0, 0, NULL); // Always undirected.
  EXEC_THROUGH_NET_EDGES_PRE(i, j, e1, { // Since i->j
      EXEC_THROUGH_FINEDGES_PRE(j, e2, k, { // and k->j
	  if(i<k) // Don't double-count.
	    WtIncEdge(i,k,1,wtnwp); // increment i-k.
	});
    });
}

U_CHANGESTAT_FN(u__osp_wtnet){
  GET_AUX_STORAGE(WtNetwork, wtnwp);
  int echange = (IS_OUTEDGE(tail, head) == 0) ? 1 : -1;

  // Update all t->h<-k shared partners.
  EXEC_THROUGH_FINEDGES_PRE(head, e, k, {
      if(tail!=k)
	WtIncEdge(tail,k,echange,wtnwp);
    });
}

F_CHANGESTAT_FN(f__osp_wtnet){
  GET_AUX_STORAGE(WtNetwork, wtnwp);

  WtNetworkDestroy(wtnwp);
}

/* Construct and maintain an undirected weighted network whose (i,j)
   value is the number of incoming shared partners of i and j. */

I_CHANGESTAT_FN(i__isp_wtnet){
  ALLOC_AUX_STORAGE(1, WtNetwork, wtnwp);
  *wtnwp = WtNetworkInitialize(NULL, NULL, NULL, 0, N_NODES, FALSE, BIPARTITE, 0, 0, NULL); // Always undirected.
  EXEC_THROUGH_NET_EDGES_PRE(i, j, e1, { // Since i->j
      EXEC_THROUGH_FOUTEDGES_PRE(i, e2, k, { // and i->k
	  if(j<k) // Don't double-count.
	    WtIncEdge(j,k,1,wtnwp); // increment j-k.
	});
    });
}

U_CHANGESTAT_FN(u__isp_wtnet){
  GET_AUX_STORAGE(WtNetwork, wtnwp);
  int echange = (IS_OUTEDGE(tail, head) == 0) ? 1 : -1;

  // Update all h<-t->k shared partners.
  EXEC_THROUGH_FOUTEDGES_PRE(tail, e, k, {
      if(head!=k)
	WtIncEdge(head,k,echange,wtnwp);
    });
}

F_CHANGESTAT_FN(f__isp_wtnet){
  GET_AUX_STORAGE(WtNetwork, wtnwp);

  WtNetworkDestroy(wtnwp);
}

/* Construct and maintain an undirected weighted network whose (i,j)
   value is the number of undirected shared partners of i and j. */

I_CHANGESTAT_FN(i__utp_wtnet){
  ALLOC_AUX_STORAGE(1, WtNetwork, wtnwp);
  *wtnwp = WtNetworkInitialize(NULL, NULL, NULL, 0, N_NODES, FALSE, BIPARTITE, 0, 0, NULL); // Always undirected.
  EXEC_THROUGH_NET_EDGES_PRE(i, j, e1, { // Since i-j
      EXEC_THROUGH_EDGES_PRE(i, e2, k, { // and i-k
	  if(j<k)
	    WtIncEdge(j,k,1,wtnwp); // increment j-k.
	});
      EXEC_THROUGH_EDGES_PRE(j, e2, k, { // and j-k
	  if(i<k)
	    WtIncEdge(i,k,1,wtnwp); // increment i-k.
	});
    });
}

U_CHANGESTAT_FN(u__utp_wtnet){
  GET_AUX_STORAGE(WtNetwork, wtnwp);
  int echange = (IS_OUTEDGE(tail, head) == 0) ? 1 : -1;

  // Update all h-t-k shared partners.
  EXEC_THROUGH_EDGES_PRE(tail, e, k, {
      if(head!=k)
	WtIncEdge(head,k,echange,wtnwp);
    });

  // Update all t-h-k shared partners.
  EXEC_THROUGH_EDGES_PRE(head, e, k, {
      if(tail!=k)
	WtIncEdge(tail,k,echange,wtnwp);
    });

}

F_CHANGESTAT_FN(f__utp_wtnet){
  GET_AUX_STORAGE(WtNetwork, wtnwp);

  WtNetworkDestroy(wtnwp);
}

