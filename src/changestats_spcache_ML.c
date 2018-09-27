#include "ergm_changestat.h"
#include "ergm_storage.h"
#include "changestats_dgw_sp_ML.h"
#include "ergm_changestat_multilayer.h"
#include "ergm_dyad_hashmap.h"

/* Construct and maintain a directed weighted network whose (i,j)
   value is the number of directed two-paths from i to j. */

I_CHANGESTAT_FN(i__otp_wtnet_ML){
  StoreDyadMapUInt *spcache = AUX_STORAGE = kh_init(DyadMapUInt);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll0, 1);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll1, 2);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll2, 3);
  bool any_order=INPUT_PARAM[4];
  
  ML_EXEC_THROUGH_NET_EDGES(ll0, i, j, e1, { // Since i->j
      ML_EXEC_THROUGH_FOUTEDGES(ll0, j, e2, k, { // and j->k
	  if(i!=k && ergm_LayerLogic2Path(i,j,j,k, ll1, ll2, any_order))
	    IncDyadMapUInt(THD(i,k),1,spcache); // increment i->k.
	});
    });
}

U_CHANGESTAT_FN(u__otp_wtnet_ML){
  GET_AUX_STORAGE(StoreDyadMapUInt, spcache);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll0, 1);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll1, 2);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll2, 3);
  bool any_order=INPUT_PARAM[4];

  int echange = (IS_OUTEDGE(tail, head) == 0) ? 1 : -1;
  
  {
    // Update all t->h->k two-paths.
    ML_EXEC_THROUGH_FOUTEDGES(ll0, head, e, k, {
	if(tail!=k && ergm_LayerLogic2Path(tail,head,head,k, ll1, ll2, any_order))
	  IncDyadMapUInt(THD(tail,k),echange,spcache);
      });
  }
  {
    // Update all k->t->h two-paths.
    ML_EXEC_THROUGH_FINEDGES(ll0, tail, e, k, {
	if(k!=head && ergm_LayerLogic2Path(k,tail,tail,head, ll1, ll2, any_order))
	  IncDyadMapUInt(THD(k,head),echange,spcache);
      });
  }
}

F_CHANGESTAT_FN(f__otp_wtnet_ML){
  GET_AUX_STORAGE(StoreDyadMapUInt, spcache);

  kh_destroy(DyadMapUInt,spcache);
  AUX_STORAGE=NULL;
}


/* Construct and maintain an undirected weighted network whose (i,j)
   value is the number of outgoing shared partners of i and j. */

I_CHANGESTAT_FN(i__osp_wtnet_ML){
  StoreDyadMapUInt *spcache = AUX_STORAGE = kh_init(DyadMapUInt);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll0, 1);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll1, 2);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll2, 3);
  bool any_order=INPUT_PARAM[4];

  ML_EXEC_THROUGH_NET_EDGES(ll0, i, j, e1, { // Since i->j
      ML_EXEC_THROUGH_FINEDGES(ll0, j, e2, k, { // and k->j
	  if(i<k && ergm_LayerLogic2Path(i,j,k,j, ll1, ll2, any_order)) // Don't double-count.
	    IncDyadMapUInt(THU(i,k),1,spcache); // increment i-k.
	});
    });
}

U_CHANGESTAT_FN(u__osp_wtnet_ML){
  GET_AUX_STORAGE(StoreDyadMapUInt, spcache);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll0, 1);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll1, 2);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll2, 3);
  bool any_order=INPUT_PARAM[4];

  int echange = (IS_OUTEDGE(tail, head) == 0) ? 1 : -1;

  // Update all t->h<-k shared partners.
  ML_EXEC_THROUGH_FINEDGES(ll0, head, e, k, {
      if(tail!=k && ergm_LayerLogic2Path(tail,head,k,head, ll1, ll2, any_order))
	IncDyadMapUInt(THU(tail,k),echange,spcache);
    });
}

F_CHANGESTAT_FN(f__osp_wtnet_ML){
  GET_AUX_STORAGE(StoreDyadMapUInt, spcache);

  kh_destroy(DyadMapUInt,spcache);
  AUX_STORAGE=NULL;
}

/* Construct and maintain an undirected weighted network whose (i,j)
   value is the number of incoming shared partners of i and j. */

I_CHANGESTAT_FN(i__isp_wtnet_ML){
  StoreDyadMapUInt *spcache = AUX_STORAGE = kh_init(DyadMapUInt);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll0, 1);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll1, 2);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll2, 3);
  bool any_order=INPUT_PARAM[4];

  ML_EXEC_THROUGH_NET_EDGES(ll0, i, j, e1, { // Since i->j
      ML_EXEC_THROUGH_FOUTEDGES(ll0, i, e2, k, { // and i->k
	  if(j<k && ergm_LayerLogic2Path(i,j,i,k, ll1, ll2, any_order)) // Don't double-count.
	    IncDyadMapUInt(THU(j,k),1,spcache); // increment j-k.
	});
    });
}

U_CHANGESTAT_FN(u__isp_wtnet_ML){
  GET_AUX_STORAGE(StoreDyadMapUInt, spcache);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll0, 1);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll1, 2);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll2, 3);
  bool any_order=INPUT_PARAM[4];

  int echange = (IS_OUTEDGE(tail, head) == 0) ? 1 : -1;

  // Update all h<-t->k shared partners.
  ML_EXEC_THROUGH_FOUTEDGES(ll0, tail, e, k, {
      if(head!=k && ergm_LayerLogic2Path(tail,head,tail,k, ll1, ll2, any_order))
	IncDyadMapUInt(THU(head,k),echange,spcache);
    });
}

F_CHANGESTAT_FN(f__isp_wtnet_ML){
  GET_AUX_STORAGE(StoreDyadMapUInt, spcache);

  kh_destroy(DyadMapUInt,spcache);
  AUX_STORAGE=NULL;
}

/* Construct and maintain an undirected weighted network whose (i,j)
   value is the number of undirected shared partners of i and j. */

I_CHANGESTAT_FN(i__utp_wtnet_ML){
  StoreDyadMapUInt *spcache = AUX_STORAGE = kh_init(DyadMapUInt);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll0, 1);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll1, 2);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll2, 3);
  bool any_order=INPUT_PARAM[4];

  ML_EXEC_THROUGH_NET_EDGES(ll0, i, j, e1, { // Since i-j
      ML_EXEC_THROUGH_EDGES(ll0, i, e2, k, { // and i-k
	  if(j<k && ergm_LayerLogic2Path(i,j,i,k, ll1, ll2, any_order))
	    IncDyadMapUInt(THU(j,k),1,spcache); // increment j-k.
	});
      ML_EXEC_THROUGH_EDGES(ll0, j, e2, k, { // and j-k
	  if(i<k && ergm_LayerLogic2Path(i,j,j,k, ll1, ll2, any_order))
	    IncDyadMapUInt(THU(i,k),1,spcache); // increment i-k.
	});
    });
}

U_CHANGESTAT_FN(u__utp_wtnet_ML){
  GET_AUX_STORAGE(StoreDyadMapUInt, spcache);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll0, 1);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll1, 2);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll2, 3);
  bool any_order=INPUT_PARAM[4];

  int echange = (IS_OUTEDGE(tail, head) == 0) ? 1 : -1;

  // Update all h-t-k shared partners.
  ML_EXEC_THROUGH_EDGES(ll0, tail, e, k, {
      if(head!=k && ergm_LayerLogic2Path(head,tail,tail,k, ll1, ll2, any_order))
	IncDyadMapUInt(THU(head,k),echange,spcache);
    });

  // Update all t-h-k shared partners.
  ML_EXEC_THROUGH_EDGES(ll0, head, e, k, {
      if(tail!=k && ergm_LayerLogic2Path(tail,head,head,k, ll1, ll2, any_order))
	IncDyadMapUInt(THU(tail,k),echange,spcache);
    });

}

F_CHANGESTAT_FN(f__utp_wtnet_ML){
  GET_AUX_STORAGE(StoreDyadMapUInt, spcache);

  kh_destroy(DyadMapUInt,spcache);
  AUX_STORAGE=NULL;
}
