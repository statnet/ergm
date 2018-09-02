#include "ergm_changestat.h"
#include "ergm_storage.h"

#include "changestats_spcache.h"


// Increment the value of the edge by inc.
/* static inline void WtIncEdge (Vertex tail, Vertex head, double inc, WtNetwork *nwp)  */
/* { */
/*   ENSURE_TH_ORDER; */

/*   if(inc==0) return; */
  
/*   // Find the out-edge */
/*   Edge oe=WtEdgetreeSearch(tail,head,nwp->outedges); */
/*   if(oe){ */
/*     // Find the corresponding in-edge. */
/*     Edge ie=WtEdgetreeSearch(head,tail,nwp->inedges); */
/*     double weight = nwp->outedges[oe].weight + inc; */
/*     if(weight==0){ */
/*       // If the function is to set the edge value to 0, just delete it. */
/*       WtDeleteEdgeFromTrees(tail,head,nwp); */
/*     }else{ */
/*       nwp->inedges[ie].weight=nwp->outedges[oe].weight = weight; */
/*     } */
/*   }else{ */
/*     // Otherwise, create a new edge with that weight 0+inc==inc. */
/*     WtAddEdgeToTrees(tail,head,inc,nwp); */
/*   } */
/* } */

#define IncEdgeMapUInt(tail, head, dir, inc, wtnwp)			\
  {									\
    if(inc!=0){								\
      unsigned int _IEMUI_val;						\
      kh_getval(EdgeMapUInt, wtnwp, TH(tail,head,dir), 0, _IEMUI_val);	\
      _IEMUI_val += inc;						\
      if(_IEMUI_val==0){kh_unset(EdgeMapUInt, wtnwp, TH(tail,head,dir));} \
      else{kh_set(EdgeMapUInt, wtnwp, TH(tail,head,dir), _IEMUI_val);}	\
    }									\
  }

void PrintEdgeMapUInt(StoreEdgeMapUInt *h){
  for(khiter_t i = kh_begin(h); i!=kh_end(h); ++i){
    if(kh_exist(h, i)){
      struct TailHead k = kh_key(h, i);
      unsigned int v = kh_val(h, i);
      Rprintf("(%d,%d)->%u\n",k.tail,k.head,v);
    }
  }
}

/* Construct and maintain a directed weighted network whose (i,j)
   value is the number of directed two-paths from i to j. */

I_CHANGESTAT_FN(i__otp_wtnet){
  StoreEdgeMapUInt *wtnwp = AUX_STORAGE = kh_init(EdgeMapUInt);
  EXEC_THROUGH_NET_EDGES(i, j, e1, { // Since i->j
      EXEC_THROUGH_FOUTEDGES(j, e2, k, { // and j->k
	  if(i!=k)
	    IncEdgeMapUInt(i,k,DIRECTED,1,wtnwp); // increment i->k.
	});
    });
}

U_CHANGESTAT_FN(u__otp_wtnet){
  GET_AUX_STORAGE(StoreEdgeMapUInt, wtnwp);
  int echange = (IS_OUTEDGE(tail, head) == 0) ? 1 : -1;

  {
    // Update all t->h->k two-paths.
    EXEC_THROUGH_FOUTEDGES(head, e, k, {
	if(tail!=k)
	  IncEdgeMapUInt(tail,k,DIRECTED,echange,wtnwp);
      });
  }
  {
    // Update all k->t->h two-paths.
    EXEC_THROUGH_FINEDGES(tail, e, k, {
	if(k!=head)
	  IncEdgeMapUInt(k,head,DIRECTED,echange,wtnwp);
      });
  }
}

F_CHANGESTAT_FN(f__otp_wtnet){
  GET_AUX_STORAGE(StoreEdgeMapUInt, wtnwp);

  kh_destroy(EdgeMapUInt,wtnwp);
  AUX_STORAGE=NULL;
}


/* Construct and maintain an undirected weighted network whose (i,j)
   value is the number of outgoing shared partners of i and j. */

I_CHANGESTAT_FN(i__osp_wtnet){
  StoreEdgeMapUInt *wtnwp = AUX_STORAGE = kh_init(EdgeMapUInt);
  EXEC_THROUGH_NET_EDGES(i, j, e1, { // Since i->j
      EXEC_THROUGH_FINEDGES(j, e2, k, { // and k->j
	  if(i<k) // Don't double-count.
	    IncEdgeMapUInt(i,k,FALSE,1,wtnwp); // increment i-k.
	});
    });
}

U_CHANGESTAT_FN(u__osp_wtnet){
  GET_AUX_STORAGE(StoreEdgeMapUInt, wtnwp);
  int echange = (IS_OUTEDGE(tail, head) == 0) ? 1 : -1;

  // Update all t->h<-k shared partners.
  EXEC_THROUGH_FINEDGES(head, e, k, {
      if(tail!=k)
	IncEdgeMapUInt(tail,k,FALSE,echange,wtnwp);
    });
}

F_CHANGESTAT_FN(f__osp_wtnet){
  GET_AUX_STORAGE(StoreEdgeMapUInt, wtnwp);

  kh_destroy(EdgeMapUInt,wtnwp);
  AUX_STORAGE=NULL;
}

/* Construct and maintain an undirected weighted network whose (i,j)
   value is the number of incoming shared partners of i and j. */

I_CHANGESTAT_FN(i__isp_wtnet){
  StoreEdgeMapUInt *wtnwp = AUX_STORAGE = kh_init(EdgeMapUInt);
  EXEC_THROUGH_NET_EDGES(i, j, e1, { // Since i->j
      EXEC_THROUGH_FOUTEDGES(i, e2, k, { // and i->k
	  if(j<k) // Don't double-count.
	    IncEdgeMapUInt(j,k,FALSE,1,wtnwp); // increment j-k.
	});
    });
}

U_CHANGESTAT_FN(u__isp_wtnet){
  GET_AUX_STORAGE(StoreEdgeMapUInt, wtnwp);
  int echange = (IS_OUTEDGE(tail, head) == 0) ? 1 : -1;

  // Update all h<-t->k shared partners.
  EXEC_THROUGH_FOUTEDGES(tail, e, k, {
      if(head!=k)
	IncEdgeMapUInt(head,k,FALSE,echange,wtnwp);
    });
}

F_CHANGESTAT_FN(f__isp_wtnet){
  GET_AUX_STORAGE(StoreEdgeMapUInt, wtnwp);

  kh_destroy(EdgeMapUInt,wtnwp);
  AUX_STORAGE=NULL;
}

/* Construct and maintain an undirected weighted network whose (i,j)
   value is the number of undirected shared partners of i and j. */

I_CHANGESTAT_FN(i__utp_wtnet){
  StoreEdgeMapUInt *wtnwp = AUX_STORAGE = kh_init(EdgeMapUInt);
  EXEC_THROUGH_NET_EDGES(i, j, e1, { // Since i-j
      EXEC_THROUGH_EDGES(i, e2, k, { // and i-k
	  if(j<k)
	    IncEdgeMapUInt(j,k,DIRECTED,1,wtnwp); // increment j-k.
	});
      EXEC_THROUGH_EDGES(j, e2, k, { // and j-k
	  if(i<k)
	    IncEdgeMapUInt(i,k,DIRECTED,1,wtnwp); // increment i-k.
	});
    });
}

U_CHANGESTAT_FN(u__utp_wtnet){
  GET_AUX_STORAGE(StoreEdgeMapUInt, wtnwp);
  int echange = (IS_OUTEDGE(tail, head) == 0) ? 1 : -1;

  // Update all h-t-k shared partners.
  EXEC_THROUGH_EDGES(tail, e, k, {
      if(head!=k)
	IncEdgeMapUInt(head,k,DIRECTED,echange,wtnwp);
    });

  // Update all t-h-k shared partners.
  EXEC_THROUGH_EDGES(head, e, k, {
      if(tail!=k)
	IncEdgeMapUInt(tail,k,DIRECTED,echange,wtnwp);
    });

}

F_CHANGESTAT_FN(f__utp_wtnet){
  GET_AUX_STORAGE(StoreEdgeMapUInt, wtnwp);

  kh_destroy(EdgeMapUInt,wtnwp);
  AUX_STORAGE=NULL;
}

