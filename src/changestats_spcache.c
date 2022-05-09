/*  File src/changestats_spcache.c in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2022 Statnet Commons
 */
#include "ergm_changestat.h"
#include "ergm_storage.h"

#include "ergm_dyad_hashmap.h"

/* Construct and maintain a directed weighted network whose (i,j)
   value is the number of directed two-paths from i to j. */

I_CHANGESTAT_FN(i__otp_wtnet){
  StoreDyadMapUInt *spcache = AUX_STORAGE = kh_init(DyadMapUInt); spcache->directed = TRUE;
  EXEC_THROUGH_NET_EDGES(i, j, e1, { // Since i->j
      EXEC_THROUGH_FOUTEDGES(j, e2, k, { // and j->k
	  if(i!=k)
	    IncDyadMapUInt(i,k,1,spcache); // increment i->k.
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
	  IncDyadMapUInt(tail,k,echange,spcache);
      });
  }
  {
    // Update all k->t->h two-paths.
    EXEC_THROUGH_FINEDGES(tail, e, k, {
	if(k!=head)
	  IncDyadMapUInt(k,head,echange,spcache);
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
  StoreDyadMapUInt *spcache = AUX_STORAGE = kh_init(DyadMapUInt); spcache->directed = FALSE;
  EXEC_THROUGH_NET_EDGES(i, j, e1, { // Since i->j
      EXEC_THROUGH_FINEDGES(j, e2, k, { // and k->j
	  if(i<k) // Don't double-count.
	    IncDyadMapUInt(i,k,1,spcache); // increment i-k.
	});
    });
}

U_CHANGESTAT_FN(u__osp_wtnet){
  GET_AUX_STORAGE(StoreDyadMapUInt, spcache);
  int echange = (IS_OUTEDGE(tail, head) == 0) ? 1 : -1;

  // Update all t->h<-k shared partners.
  EXEC_THROUGH_FINEDGES(head, e, k, {
      if(tail!=k)
	IncDyadMapUInt(tail,k,echange,spcache);
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
  StoreDyadMapUInt *spcache = AUX_STORAGE = kh_init(DyadMapUInt); spcache->directed = FALSE;
  EXEC_THROUGH_NET_EDGES(i, j, e1, { // Since i->j
      EXEC_THROUGH_FOUTEDGES(i, e2, k, { // and i->k
	  if(j<k) // Don't double-count.
	    IncDyadMapUInt(j,k,1,spcache); // increment j-k.
	});
    });
}

U_CHANGESTAT_FN(u__isp_wtnet){
  GET_AUX_STORAGE(StoreDyadMapUInt, spcache);
  int echange = (IS_OUTEDGE(tail, head) == 0) ? 1 : -1;

  // Update all h<-t->k shared partners.
  EXEC_THROUGH_FOUTEDGES(tail, e, k, {
      if(head!=k)
	IncDyadMapUInt(head,k,echange,spcache);
    });
}

F_CHANGESTAT_FN(f__isp_wtnet){
  GET_AUX_STORAGE(StoreDyadMapUInt, spcache);

  kh_destroy(DyadMapUInt,spcache);
  AUX_STORAGE=NULL;
}

/* Construct and maintain an undirected weighted network whose (i,j)
   value is the number of reciprocated partners of i and j. */

I_CHANGESTAT_FN(i__rtp_wtnet){
  StoreDyadMapUInt *spcache = AUX_STORAGE = kh_init(DyadMapUInt); spcache->directed = FALSE;
  EXEC_THROUGH_NET_EDGES(i, j, e1, { // Since i->j
      if(IS_OUTEDGE(j,i)) // and j->i
        EXEC_THROUGH_FOUTEDGES(i, e2, k, { // and i->k
            if(j<k&&IS_OUTEDGE(k,i)) // and k->i (and don't double-count)
              IncDyadMapUInt(j,k,1,spcache); // increment j-k.
          });
    });
}

U_CHANGESTAT_FN(u__rtp_wtnet){
  GET_AUX_STORAGE(StoreDyadMapUInt, spcache);
  if(!IS_OUTEDGE(head,tail)) return; // If no reciprocating edge, no effect.
  int echange = (IS_OUTEDGE(tail, head) == 0) ? 1 : -1;

  // Update all h?->t<->k shared partners.
  EXEC_THROUGH_FOUTEDGES(tail, e, k, {
      if(head!=k&&IS_OUTEDGE(k,tail))
	IncDyadMapUInt(head,k,echange,spcache);
    });
  // Update all k<->h?->t shared partners.
  EXEC_THROUGH_FOUTEDGES(head, e, k, {
      if(tail!=k&&IS_OUTEDGE(k,head))
	IncDyadMapUInt(tail,k,echange,spcache);
    });
}

F_CHANGESTAT_FN(f__rtp_wtnet){
  GET_AUX_STORAGE(StoreDyadMapUInt, spcache);

  kh_destroy(DyadMapUInt,spcache);
  AUX_STORAGE=NULL;
}

/* Construct and maintain an undirected weighted network whose (i,j)
   value is the number of undirected shared partners of i and j. */

I_CHANGESTAT_FN(i__utp_wtnet){
  StoreDyadMapUInt *spcache = AUX_STORAGE = kh_init(DyadMapUInt); spcache->directed = FALSE;
  EXEC_THROUGH_NET_EDGES(i, j, e1, { // Since i-j
      EXEC_THROUGH_EDGES(i, e2, k, { // and i-k
	  if(j<k)
	    IncDyadMapUInt(j,k,1,spcache); // increment j-k.
	});
      EXEC_THROUGH_EDGES(j, e2, k, { // and j-k
	  if(i<k)
	    IncDyadMapUInt(i,k,1,spcache); // increment i-k.
	});
    });
}

U_CHANGESTAT_FN(u__utp_wtnet){
  GET_AUX_STORAGE(StoreDyadMapUInt, spcache);
  int echange = (IS_OUTEDGE(tail, head) == 0) ? 1 : -1;

  // Update all h-t-k shared partners.
  EXEC_THROUGH_EDGES(tail, e, k, {
      if(head!=k)
	IncDyadMapUInt(head,k,echange,spcache);
    });

  // Update all t-h-k shared partners.
  EXEC_THROUGH_EDGES(head, e, k, {
      if(tail!=k)
	IncDyadMapUInt(tail,k,echange,spcache);
    });

}

F_CHANGESTAT_FN(f__utp_wtnet){
  GET_AUX_STORAGE(StoreDyadMapUInt, spcache);

  kh_destroy(DyadMapUInt,spcache);
  AUX_STORAGE=NULL;
}

