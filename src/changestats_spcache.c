#include "ergm_changestat.h"
#include "ergm_wtedgetree.h"
#include "ergm_storage.h"


/* Construct and maintain a directed weighted network whose (i,j)
   value is the number of directed two-paths from i to j. */

I_CHANGESTAT_FN(i__otp_wtnet){
  ALLOC_AUX_STORAGE(1, WtNetwork, wtnwp);
  *wtnwp = WtNetworkInitialize(NULL, NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, 0, 0, NULL);
  EXEC_THROUGH_NET_EDGES(i, j, e1, { // Since i->j
      EXEC_THROUGH_FOUTEDGES(j, e2, k, { // and j->k
	  if(i!=k)
	    WtSetEdge(i,k,WtGetEdge(i,k,wtnwp)+1,wtnwp); // increment i->k.
	});
    });
}

U_CHANGESTAT_FN(u__otp_wtnet){
  GET_AUX_STORAGE(WtNetwork, wtnwp);
  int echange = (IS_OUTEDGE(tail, head) == 0) ? 1 : -1;

  {
    // Update all t->h->k two-paths.
    EXEC_THROUGH_FOUTEDGES(head, e, k, {
	if(tail!=k)
	  WtSetEdge(tail,k,WtGetEdge(tail,k,wtnwp)+echange,wtnwp);
      });
  }
  {
    // Update all k->t->h two-paths.
    EXEC_THROUGH_FINEDGES(tail, e, k, {
	if(k!=head)
	  WtSetEdge(k,head,WtGetEdge(k,head,wtnwp)+echange,wtnwp);
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
  EXEC_THROUGH_NET_EDGES(i, j, e1, { // Since i->j
      EXEC_THROUGH_FINEDGES(j, e2, k, { // and k->j
	  if(i<k) // Don't double-count.
	    WtSetEdge(i,k,WtGetEdge(i,k,wtnwp)+1,wtnwp); // increment i-k.
	});
    });
}

U_CHANGESTAT_FN(u__osp_wtnet){
  GET_AUX_STORAGE(WtNetwork, wtnwp);
  int echange = (IS_OUTEDGE(tail, head) == 0) ? 1 : -1;

  // Update all t->h<-k shared partners.
  EXEC_THROUGH_FINEDGES(head, e, k, {
      if(tail!=k)
	WtSetEdge(tail,k,WtGetEdge(tail,k,wtnwp)+echange,wtnwp);
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
  EXEC_THROUGH_NET_EDGES(i, j, e1, { // Since i->j
      EXEC_THROUGH_FOUTEDGES(i, e2, k, { // and i->k
	  if(j<k) // Don't double-count.
	    WtSetEdge(j,k,WtGetEdge(j,k,wtnwp)+1,wtnwp); // increment j-k.
	});
    });
}

U_CHANGESTAT_FN(u__isp_wtnet){
  GET_AUX_STORAGE(WtNetwork, wtnwp);
  int echange = (IS_OUTEDGE(tail, head) == 0) ? 1 : -1;

  // Update all h<-t->k shared partners.
  EXEC_THROUGH_FOUTEDGES(tail, e, k, {
      if(head!=k)
	WtSetEdge(head,k,WtGetEdge(head,k,wtnwp)+echange,wtnwp);
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
  EXEC_THROUGH_NET_EDGES(i, j, e1, { // Since i-j
      EXEC_THROUGH_EDGES(i, e2, k, { // and i-k
	  if(j<k)
	    WtSetEdge(j,k,WtGetEdge(j,k,wtnwp)+1,wtnwp); // increment j-k.
	});
      EXEC_THROUGH_EDGES(j, e2, k, { // and j-k
	  if(i<k)
	    WtSetEdge(i,k,WtGetEdge(i,k,wtnwp)+1,wtnwp); // increment i-k.
	});
    });
}

U_CHANGESTAT_FN(u__utp_wtnet){
  GET_AUX_STORAGE(WtNetwork, wtnwp);
  int echange = (IS_OUTEDGE(tail, head) == 0) ? 1 : -1;

  // Update all h-t-k shared partners.
  EXEC_THROUGH_EDGES(tail, e, k, {
      if(head!=k)
	WtSetEdge(head,k,WtGetEdge(head,k,wtnwp)+echange,wtnwp);
    });

  // Update all t-h-k shared partners.
  EXEC_THROUGH_EDGES(head, e, k, {
      if(tail!=k)
	WtSetEdge(tail,k,WtGetEdge(tail,k,wtnwp)+echange,wtnwp);
    });

}

F_CHANGESTAT_FN(f__utp_wtnet){
  GET_AUX_STORAGE(WtNetwork, wtnwp);

  WtNetworkDestroy(wtnwp);
}

