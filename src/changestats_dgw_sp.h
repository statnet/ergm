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

static ALWAYS_INLINE Vertex count_utp(ErgmCppNetwork& nw, Vertex tail, Vertex head, StoreStrictDyadMapUInt *spcache){
  if(spcache) return GETUDMUI(tail, head, spcache);
  Vertex count = 0;
  // h - v - t
  for(auto v: nw.neighbors(head)) count += nw(v, tail);
  return count;
}

static ALWAYS_INLINE std::tuple<Vertex, Vertex> count_utp2(ErgmCppNetwork& nw, Vertex tail, Vertex head1, Vertex head2, StoreStrictDyadMapUInt *spcache){
  if(spcache) return {
      GETUDMUI(tail, head1, spcache),
      GETUDMUI(tail, head2, spcache)
    };
  Vertex count1 = 0, count2 = 0;
  // h - v - t
  for(auto v: nw.neighbors(tail)) {
    count1 += nw(v, head1);
    count2 += nw(v, head2);
  }
  return {count1, count2};
}

static ALWAYS_INLINE Vertex count_otp(ErgmCppNetwork& nw, Vertex tail, Vertex head, StoreStrictDyadMapUInt *spcache){
  if(spcache) return GETDDMUI(tail, head, spcache);
  Vertex count = 0;
  // t -> v -> h
  for(auto v: nw.out_neighbors(tail)) count += nw(v, head);
  return count;
}

/* itp == otp with tail and head swapped. */

static ALWAYS_INLINE Vertex count_rtp(ErgmCppNetwork& nw, Vertex tail, Vertex head, Vertex exclude1, StoreStrictDyadMapUInt *spcache){
  if(spcache) return GETUDMUI(tail, head, spcache);
  Vertex count = 0;
  // t <-> v <-> h
  for(auto v: nw.out_neighbors(tail))
    if(v != exclude1 && v != head && nw(v, tail))
      count += nw(v, head) && nw(head, v);
  return count;
}

static ALWAYS_INLINE Vertex count_osp(ErgmCppNetwork& nw, Vertex tail, Vertex head, StoreStrictDyadMapUInt *spcache){
  if(spcache) return GETUDMUI(tail, head, spcache);
  Vertex count = 0;
  // h -> v <- t
  for(auto v: nw.out_neighbors(head))
    if(v != tail) count += nw(tail, v);
  return count;
}

static ALWAYS_INLINE Vertex count_isp(ErgmCppNetwork& nw, Vertex tail, Vertex head, StoreStrictDyadMapUInt *spcache){
  if(spcache) return GETUDMUI(tail, head, spcache);
  Vertex count = 0;
  // h <- v -> t
  for(auto v: nw.in_neighbors(head))
    if(v != tail) count += nw(v, tail);
  return count;
}


template<typename UpdatePath>
static ALWAYS_INLINE void dsp_change(L2Type type, Vertex tail, Vertex head, Network *nwp, StoreStrictDyadMapUInt *spcache, UpdatePath update_path){

  if (type == L2UTP) {
    auto nw = ErgmCppNetwork(nwp, false);

    // h - u - #v - t: focus dyad: t - u
    for(auto u: nw.neighbors(head))
      if(u != tail)
        update_path(count_utp(nw, tail, u, spcache));

    // t - u - #v - h: focus dyad: h - u
    for(auto u: nw.neighbors(tail))
      if(u != head)
        update_path(count_utp(nw, u, head, spcache));

  } else {
    auto nw = ErgmCppNetwork(nwp, true);
    switch(type){
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
}


template<typename UpdatePath, typename UpdateFocus>
static ALWAYS_INLINE void esp_change(L2Type type, Vertex tail, Vertex head, Network *nwp, StoreStrictDyadMapUInt *spcache, UpdatePath update_path, UpdateFocus update_focus){

  unsigned int L2th;

  if (type == L2UTP) {
    auto nw = ErgmCppNetwork(nwp, false);

    L2th = spcache ? GETDDMUI(tail, head, spcache) : 0;

    // h - u - t: focus dyad t - h
    // h - u - t - #v - u: focus dyad t - u
    // u - #v - h - u - t: focus dyad h - u
    for(auto u: nw.neighbors(head))
      if(nw(u, tail)){
        if(!spcache) L2th++;
        auto utp_counts = count_utp2(nw, u, tail, head, spcache);
        update_path(std::get<0>(utp_counts));
        update_path(std::get<1>(utp_counts));
      }

  } else {
    auto nw = ErgmCppNetwork(nwp, true);
    bool htedge;

    switch(type){
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
  }

  update_focus(L2th);
}

static ALWAYS_INLINE int dsp_nonzero_change(L2Type type, Vertex tail, Vertex head, Network *nwp, Rboolean edgestate, StoreStrictDyadMapUInt *spcache){
  int echange = edgestate ? -1 : 1;
  int delta = 0;
  dsp_change(type, tail, head, nwp, spcache,
             [&](int L2){ delta += (L2 + echange != 0) - (L2 != 0); });
  return delta;
}

} // namespace sp
} // namespace v1
} // namespace ergm
#endif
