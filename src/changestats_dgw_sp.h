/*  File src/changestats_dgw_sp.h in package ergm, part of the Statnet suite of
 *  packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2026 Statnet Commons
 */
#ifndef _CHANGESTATS_DGW_SP_H_
#define _CHANGESTATS_DGW_SP_H_

#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_dyad_hashmap.h"

typedef enum {L2UTP, L2OTP, L2ITP, L2RTP, L2OSP, L2ISP} L2Type;

#include "cpp/ergm_network.h"

namespace ergm {
inline namespace v1 {
namespace sp {

/* The following functions calculate or obtain the number of two-paths
   of specified type from tail to head. */

inline int count_utp(ErgmCppNetwork& nw, Vertex tail, Vertex head, StoreStrictDyadMapUInt *spcache){
  if(spcache) return GETUDMUI(tail, head, spcache);
  int count = 0;
  // h - v - t
  for(auto v: nw.neighbors(head)) count += nw(v, tail);
  return count;
}

inline int count_otp(ErgmCppNetwork& nw, Vertex tail, Vertex head, StoreStrictDyadMapUInt *spcache){
  if(spcache) return GETDDMUI(tail, head, spcache);
  int count = 0;
  // t -> v -> h
  for(auto v: nw.out_neighbors(tail)) count += nw(v, head);
  return count;
}

/* itp == otp with tail and head swapped. */

inline int count_rtp(ErgmCppNetwork& nw, Vertex tail, Vertex head, Vertex exclude1, StoreStrictDyadMapUInt *spcache){
  if(spcache) return GETUDMUI(tail, head, spcache);
  int count = 0;
  // t <-> v <-> h
  for(auto v: nw.out_neighbors(tail))
    if(v != exclude1 && v != head && nw(v, tail))
      count += nw(v, head) && nw(head, v);
  return count;
}

inline int count_osp(ErgmCppNetwork& nw, Vertex tail, Vertex head, StoreStrictDyadMapUInt *spcache){
  if(spcache) return GETUDMUI(tail, head, spcache);
  int count = 0;
  // h -> v <- t
  for(auto v: nw.out_neighbors(head))
    if(v != tail) count += nw(tail, v);
  return count;
}

inline int count_isp(ErgmCppNetwork& nw, Vertex tail, Vertex head, StoreStrictDyadMapUInt *spcache){
  if(spcache) return GETUDMUI(tail, head, spcache);
  int count = 0;
  // h <- v -> t
  for(auto v: nw.in_neighbors(head))
    if(v != tail) count += nw(v, tail);
  return count;
}


template<typename UpdatePath, typename UpdateFocus>
inline void dsp_change(L2Type type, Vertex tail, Vertex head, ErgmCppNetwork& nw, StoreStrictDyadMapUInt *spcache, UpdatePath update_path, UpdateFocus){
  switch(type){
  case L2UTP:
    // h - u - #v - t: focus dyad: t - u
    for(auto u: nw.neighbors(head))
      if(u != tail)
        update_path(count_utp(nw, tail, u, spcache));

    // t - u - #v - h: focus dyad: h - u
    for(auto u: nw.neighbors(tail))
      if(u != head)
        update_path(count_utp(nw, u, head, spcache));
    break;

  case L2OTP:
  case L2ITP: // Identical without a focus dyad to give direction.
    // h -> u <- #v <- t: focus dyad t - u
    for(auto u: nw.out_neighbors(head))
      if(u != tail)
        update_path(count_otp(nw, tail, u, spcache));

    // t <- u -> #v -> h: focus dyad h - u
    for(auto u: nw.in_neighbors(tail))
      if(u != head)
        update_path(count_otp(nw, u, head, spcache));
    break;

  case L2RTP:
    if(nw(head, tail)){
      // t <-> u <-> #v <-> h: focus dyad h - u
      for(auto u: nw.out_neighbors(tail))
        if(u != head && nw(u, tail))
          update_path(count_rtp(nw, u, head, tail, spcache));

      // h <-> u <-> #v <-> t: focus dyad t - u
      for(auto u: nw.out_neighbors(head))
        if(u != tail && nw(u, head))
          update_path(count_rtp(nw, u, tail, head, spcache));
    }
    break;

  case L2OSP:
    // h <- u -> #v <- t: focus dyad t - u
    for(auto u: nw.in_neighbors(head))
      if(u != tail)
        update_path(count_osp(nw, tail, u, spcache));
    break;

  case L2ISP:
    // t -> u <- #v -> h: focus dyad h - u
    for(auto u: nw.out_neighbors(tail))
      if(u != head)
        update_path(count_isp(nw, head, u, spcache));
    break;

  default: error("In ergm shared partner helper, an unsupported type of triad: %d.", type);
  }
}


template<typename UpdatePath, typename UpdateFocus>
inline void esp_change(L2Type type, Vertex tail, Vertex head, ErgmCppNetwork& nw, StoreStrictDyadMapUInt *spcache, UpdatePath update_path, UpdateFocus update_focus){
  int L2th;
  bool htedge;

  switch(type){
  case L2UTP:
    L2th = spcache ? GETDDMUI(tail, head, spcache) : 0;

    // h - u - t: focus dyad t - h
    // h - u - t - #v - u: focus dyad t - u
    // u - #v - h - u - t: focus dyad h - u
    for(auto u: nw.neighbors(head))
      if(nw(u, tail)){
        if(!spcache) L2th++;
        update_path(count_utp(nw, tail, u, spcache));
        update_path(count_utp(nw, u, head, spcache));
      }
    break;

  case L2OTP:
    L2th = spcache ? GETDDMUI(tail, head, spcache) : 0;

    // t -> u -> h: focus dyad t -> h
    // u <- #v <- t -> u <- h: focus dyad t -> u
    for(auto u: nw.out_neighbors(tail)){
      if(!spcache && u != head && nw(u, head)) L2th++;
      if(u != head && nw(head, u))
        update_path(count_otp(nw, tail, u, spcache));
    }

    // u -> #v -> h <- u -> h: focus dyad u -> h
    for(auto u: nw.in_neighbors(head))
      if(u != tail && nw(u, tail))
        update_path(count_otp(nw, u, head, spcache));
    break;

  case L2ITP:
    L2th = spcache ? GETDDMUI(head, tail, spcache) : 0;

    // h -> u -> t: focus dyad t -> h
    // u -> #v -> h -> u -> t: focus dyad h -> u
    for(auto u: nw.out_neighbors(head))
      if(u != tail && nw(u, tail)){
        if(!spcache) L2th++;
        update_path(count_otp(nw, u, head, spcache));
      }

    // u <- #v <- t <- u <- h: focus dyad u -> t
    for(auto u: nw.in_neighbors(tail))
      if(u != head && nw(head, u))
        update_path(count_otp(nw, tail, u, spcache));
    break;

  case L2RTP:
    L2th = spcache ? GETUDMUI(tail, head, spcache) : 0;
    htedge = nw(head, tail);

    // t <-> u <-> h: focus dyad t -> h
    // u <-> h -> t <- u <-> #v <-> t: focus dyad t <- u
    for(auto u: nw.in_neighbors(tail)){
      if(u != head){
        if(!spcache) L2th += nw(tail, u) && nw(head, u) && nw(u, head);
        if(htedge && nw(head, u) && nw(u, head))
          update_path(count_rtp(nw, u, tail, tail, spcache));
      }
    }

    // u <-> h -> t -> u <-> #v <-> t: focus dyad t -> u
    for(auto u: nw.out_neighbors(tail))
      if(u != head && htedge && nw(head, u) && nw(u, head))
        update_path(count_rtp(nw, u, tail, tail, spcache));

    // u <-> t <- h <- u <-> #v <-> h: focus dyad u -> h
    for(auto u: nw.in_neighbors(head))
      if(u != tail && htedge && nw(tail, u) && nw(u, tail))
        update_path(count_rtp(nw, u, head, head, spcache));

    // u <-> t -> h <- u <-> #v <-> h: focus dyad u <- h
    for(auto u: nw.out_neighbors(head))
      if(u != tail && htedge && nw(tail, u) && nw(u, tail))
        update_path(count_rtp(nw, u, head, head, spcache));
    break;

  case L2OSP:
    L2th = spcache ? GETUDMUI(tail, head, spcache) : 0;

    // t -> u <- h: focus dyad t -> h
    // u -> #v <- t -> u -> h: focus dyad t -> u
    for(auto u: nw.out_neighbors(tail))
      if(u != head){
        if(!spcache) L2th += nw(head, u);
        if(nw(u, head))
          update_path(count_osp(nw, tail, u, spcache));
      }

    // u -> #v <- t <- u -> h: focus dyad u -> t
    for(auto u: nw.in_neighbors(tail))
      if(u != head && nw(u, head))
        update_path(count_osp(nw, tail, u, spcache));
    break;

  case L2ISP:
    L2th = spcache ? GETUDMUI(tail, head, spcache) : 0;

    // h <- u -> t: focus dyad t -> h
    // u <- #v -> h <- u <- t: focus dyad  u -> h
    for(auto u: nw.in_neighbors(head))
      if(u != tail){
        if(!spcache) L2th += nw(u, tail);
        if(nw(tail, u))
          update_path(count_isp(nw, head, u, spcache));
      }

    // u <- #v -> h -> u <- t: focus dyad  h -> u
    for(auto u: nw.out_neighbors(head))
      if(u != tail && nw(tail, u))
        update_path(count_isp(nw, head, u, spcache));
    break;

  default: error("In ergm shared partner helper, an unsupported type of triad: %d.", type);
  }

  update_focus(L2th);
}

inline int dsp_nonzero_change(L2Type type, Vertex tail, Vertex head, ErgmCppNetwork& nw, Rboolean edgestate, StoreStrictDyadMapUInt *spcache){
  int echange = edgestate ? -1 : 1;
  int delta = 0;
  dsp_change(type, tail, head, nw, spcache,
             [&](int L2){ delta += (L2 + echange != 0) - (L2 != 0); },
             [&](int){});
  return delta;
}

} // namespace sp
} // namespace v1
} // namespace ergm
#endif
