#include "ergm_changestat.h"
#include "ergm_storage.h"

#include "ergm_dyad_hashmap.h"


// Increment the value in the hashtable by inc.
static inline void IncDyadMapUInt(struct TailHead th, int inc, StoreDyadMapUInt *spcache){
  if(inc!=0){
    khiter_t pos = kh_get(DyadMapUInt, spcache, th);
    unsigned int val = pos==kh_end(spcache) ? 0 : kh_value(spcache, pos);
    val += inc;
    if(val==0){
      kh_del(DyadMapUInt, spcache, pos);
    }else{
      if(pos==kh_end(spcache)){
	int ret;
	pos = kh_put(DyadMapUInt, spcache, th, &ret);
      }
      kh_val(spcache, pos) = val;
    }
  }
}

void PrintDyadMapUInt(StoreDyadMapUInt *h){
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
  StoreDyadMapUInt *spcache = AUX_STORAGE = kh_init(DyadMapUInt);
  EXEC_THROUGH_NET_EDGES(i, j, e1, { // Since i->j
      EXEC_THROUGH_FOUTEDGES(j, e2, k, { // and j->k
	  if(i!=k)
	    IncDyadMapUInt(THD(i,k),1,spcache); // increment i->k.
	});
    });
}

U_CHANGESTAT_FN(u__otp_wtnet){
  GET_AUX_STORAGE(StoreDyadMapUInt, spcache);
  int echange = (IS_OUTEDGE(tail, head) == 0) ? 1 : -1;

  {
    // Update all t->h->k two-paths.
    EXEC_THROUGH_FOUTEDGES(head, e, k, {
	if(tail!=k)
	  IncDyadMapUInt(THD(tail,k),echange,spcache);
      });
  }
  {
    // Update all k->t->h two-paths.
    EXEC_THROUGH_FINEDGES(tail, e, k, {
	if(k!=head)
	  IncDyadMapUInt(THD(k,head),echange,spcache);
      });
  }
}

F_CHANGESTAT_FN(f__otp_wtnet){
  GET_AUX_STORAGE(StoreDyadMapUInt, spcache);

  kh_destroy(DyadMapUInt,spcache);
  AUX_STORAGE=NULL;
}


/* Construct and maintain an undirected weighted network whose (i,j)
   value is the number of outgoing shared partners of i and j. */

I_CHANGESTAT_FN(i__osp_wtnet){
  StoreDyadMapUInt *spcache = AUX_STORAGE = kh_init(DyadMapUInt);
  EXEC_THROUGH_NET_EDGES(i, j, e1, { // Since i->j
      EXEC_THROUGH_FINEDGES(j, e2, k, { // and k->j
	  if(i<k) // Don't double-count.
	    IncDyadMapUInt(THU(i,k),1,spcache); // increment i-k.
	});
    });
}

U_CHANGESTAT_FN(u__osp_wtnet){
  GET_AUX_STORAGE(StoreDyadMapUInt, spcache);
  int echange = (IS_OUTEDGE(tail, head) == 0) ? 1 : -1;

  // Update all t->h<-k shared partners.
  EXEC_THROUGH_FINEDGES(head, e, k, {
      if(tail!=k)
	IncDyadMapUInt(THU(tail,k),echange,spcache);
    });
}

F_CHANGESTAT_FN(f__osp_wtnet){
  GET_AUX_STORAGE(StoreDyadMapUInt, spcache);

  kh_destroy(DyadMapUInt,spcache);
  AUX_STORAGE=NULL;
}

/* Construct and maintain an undirected weighted network whose (i,j)
   value is the number of incoming shared partners of i and j. */

I_CHANGESTAT_FN(i__isp_wtnet){
  StoreDyadMapUInt *spcache = AUX_STORAGE = kh_init(DyadMapUInt);
  EXEC_THROUGH_NET_EDGES(i, j, e1, { // Since i->j
      EXEC_THROUGH_FOUTEDGES(i, e2, k, { // and i->k
	  if(j<k) // Don't double-count.
	    IncDyadMapUInt(THU(j,k),1,spcache); // increment j-k.
	});
    });
}

U_CHANGESTAT_FN(u__isp_wtnet){
  GET_AUX_STORAGE(StoreDyadMapUInt, spcache);
  int echange = (IS_OUTEDGE(tail, head) == 0) ? 1 : -1;

  // Update all h<-t->k shared partners.
  EXEC_THROUGH_FOUTEDGES(tail, e, k, {
      if(head!=k)
	IncDyadMapUInt(THU(head,k),echange,spcache);
    });
}

F_CHANGESTAT_FN(f__isp_wtnet){
  GET_AUX_STORAGE(StoreDyadMapUInt, spcache);

  kh_destroy(DyadMapUInt,spcache);
  AUX_STORAGE=NULL;
}

/* Construct and maintain an undirected weighted network whose (i,j)
   value is the number of undirected shared partners of i and j. */

I_CHANGESTAT_FN(i__utp_wtnet){
  StoreDyadMapUInt *spcache = AUX_STORAGE = kh_init(DyadMapUInt);
  EXEC_THROUGH_NET_EDGES(i, j, e1, { // Since i-j
      EXEC_THROUGH_EDGES(i, e2, k, { // and i-k
	  if(j<k)
	    IncDyadMapUInt(THU(j,k),1,spcache); // increment j-k.
	});
      EXEC_THROUGH_EDGES(j, e2, k, { // and j-k
	  if(i<k)
	    IncDyadMapUInt(THU(i,k),1,spcache); // increment i-k.
	});
    });
}

U_CHANGESTAT_FN(u__utp_wtnet){
  GET_AUX_STORAGE(StoreDyadMapUInt, spcache);
  int echange = (IS_OUTEDGE(tail, head) == 0) ? 1 : -1;

  // Update all h-t-k shared partners.
  EXEC_THROUGH_EDGES(tail, e, k, {
      if(head!=k)
	IncDyadMapUInt(THU(head,k),echange,spcache);
    });

  // Update all t-h-k shared partners.
  EXEC_THROUGH_EDGES(head, e, k, {
      if(tail!=k)
	IncDyadMapUInt(THU(tail,k),echange,spcache);
    });

}

F_CHANGESTAT_FN(f__utp_wtnet){
  GET_AUX_STORAGE(StoreDyadMapUInt, spcache);

  kh_destroy(DyadMapUInt,spcache);
  AUX_STORAGE=NULL;
}

