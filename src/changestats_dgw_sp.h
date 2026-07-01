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

template<typename NetworkView>
inline int count_otp(NetworkView& nw, Vertex tail, Vertex head){
  int count = 0;
  for(auto u: nw.in_neighbors(head)) count += nw(tail, u);
  return count;
}

template<typename NetworkView>
inline int count_utp(NetworkView& nw, Vertex tail, Vertex head){
  int count = 0;
  for(auto u: nw.neighbors(head)) count += nw(u, tail);
  return count;
}

template<typename NetworkView>
inline int count_osp(NetworkView& nw, Vertex tail, Vertex head){
  int count = 0;
  for(auto u: nw.out_neighbors(head))
    if(u != tail) count += nw(tail, u);
  return count;
}

template<typename NetworkView>
inline int count_isp(NetworkView& nw, Vertex tail, Vertex head){
  int count = 0;
  for(auto u: nw.in_neighbors(head))
    if(u != tail) count += nw(u, tail);
  return count;
}

template<typename NetworkView>
inline int count_rtp(NetworkView& nw, Vertex tail, Vertex head, Vertex exclude1, Vertex exclude2){
  int count = 0;
  for(auto u: nw.out_neighbors(tail))
    if(u != exclude1 && u != exclude2 && nw(u, tail))
      count += nw(u, head) && nw(head, u);
  return count;
}

template<L2Type type, typename NetworkView, typename UpdatePath, typename UpdateFocus>
inline void dsp_change(Vertex tail, Vertex head, NetworkView& nw, StoreStrictDyadMapUInt *spcache, UpdatePath update_path, UpdateFocus){
  if constexpr(type == L2UTP){
    for(auto u: nw.neighbors(head))
      if(u != tail)
        update_path(spcache ? GETUDMUI(tail, u, spcache) : count_utp(nw, tail, u));

    for(auto u: nw.neighbors(tail))
      if(u != head)
        update_path(spcache ? GETUDMUI(u, head, spcache) : count_utp(nw, u, head));
  }else if constexpr(type == L2OTP){
    for(auto k: nw.out_neighbors(head))
      if(k != tail)
        update_path(spcache ? GETDDMUI(tail, k, spcache) : count_otp(nw, tail, k));

    for(auto k: nw.in_neighbors(tail))
      if(k != head)
        update_path(spcache ? GETDDMUI(k, head, spcache) : count_otp(nw, k, head));
  }else if constexpr(type == L2ITP){
    for(auto k: nw.out_neighbors(head))
      if(k != tail)
        update_path(spcache ? GETDDMUI(tail, k, spcache) : count_otp(nw, tail, k));

    for(auto k: nw.in_neighbors(tail))
      if(k != head)
        update_path(spcache ? GETDDMUI(k, head, spcache) : count_otp(nw, k, head));
  }else if constexpr(type == L2RTP){
    if(nw(head, tail)){
      for(auto k: nw.out_neighbors(tail))
        if(k != head && nw(k, tail))
          update_path(spcache ? GETUDMUI(k, head, spcache) : count_rtp(nw, k, head, tail, head));

      for(auto k: nw.out_neighbors(head))
        if(k != tail && nw(k, head))
          update_path(spcache ? GETUDMUI(k, tail, spcache) : count_rtp(nw, k, tail, head, tail));
    }
  }else if constexpr(type == L2OSP){
    for(auto k: nw.in_neighbors(head))
      if(k != tail)
        update_path(spcache ? GETUDMUI(tail, k, spcache) : count_osp(nw, tail, k));
  }else if constexpr(type == L2ISP){
    for(auto k: nw.out_neighbors(tail))
      if(k != head)
        update_path(spcache ? GETUDMUI(k, head, spcache) : count_isp(nw, head, k));
  }
}

template<typename NetworkView, typename UpdatePath, typename UpdateFocus>
inline void dsp_change(L2Type type, Vertex tail, Vertex head, NetworkView& nw, StoreStrictDyadMapUInt *spcache, UpdatePath update_path, UpdateFocus update_focus){
  switch(type){
  case L2UTP: dsp_change<L2UTP>(tail, head, nw, spcache, update_path, update_focus); break;
  case L2OTP: dsp_change<L2OTP>(tail, head, nw, spcache, update_path, update_focus); break;
  case L2ITP: dsp_change<L2ITP>(tail, head, nw, spcache, update_path, update_focus); break;
  case L2RTP: dsp_change<L2RTP>(tail, head, nw, spcache, update_path, update_focus); break;
  case L2OSP: dsp_change<L2OSP>(tail, head, nw, spcache, update_path, update_focus); break;
  case L2ISP: dsp_change<L2ISP>(tail, head, nw, spcache, update_path, update_focus); break;
  default: error("In ergm shared partner helper, an unsupported type of triad: %d.", type);
  }
}

template<L2Type type, typename NetworkView, typename UpdatePath, typename UpdateFocus>
inline void esp_change(Vertex tail, Vertex head, NetworkView& nw, StoreStrictDyadMapUInt *spcache, UpdatePath update_path, UpdateFocus update_focus){
  if constexpr(type == L2UTP){
    int L2th = spcache ? GETDDMUI(tail, head, spcache) : 0;

    for(auto u: nw.neighbors(head))
      if(nw(u, tail)){
        if(!spcache) L2th++;
        update_path(spcache ? GETUDMUI(tail, u, spcache) : count_utp(nw, tail, u));
        update_path(spcache ? GETUDMUI(u, head, spcache) : count_utp(nw, u, head));
      }

    update_focus(L2th);
  }else if constexpr(type == L2OTP){
    int L2th = spcache ? GETDDMUI(tail, head, spcache) : 0;

    for(auto k: nw.out_neighbors(tail)){
      if(!spcache && k != head && nw(k, head)) L2th++;
      if(k != head && nw(head, k))
        update_path(spcache ? GETDDMUI(tail, k, spcache) : count_otp(nw, tail, k));
    }

    for(auto k: nw.in_neighbors(head))
      if(k != tail && nw(k, tail))
        update_path(spcache ? GETDDMUI(k, head, spcache) : count_otp(nw, k, head));

    update_focus(L2th);
  }else if constexpr(type == L2ITP){
    int L2th = spcache ? GETDDMUI(head, tail, spcache) : 0;

    for(auto k: nw.out_neighbors(head))
      if(k != tail && nw(k, tail)){
        if(!spcache) L2th++;
        update_path(spcache ? GETDDMUI(k, head, spcache) : count_otp(nw, k, head));
      }

    for(auto k: nw.in_neighbors(tail))
      if(k != head && nw(head, k))
        update_path(spcache ? GETDDMUI(tail, k, spcache) : count_otp(nw, tail, k));

    update_focus(L2th);
  }else if constexpr(type == L2RTP){
    int L2th = spcache ? GETDDMUI(tail, head, spcache) : 0;
    bool htedge = nw(head, tail);

    for(auto k: nw.in_neighbors(tail)){
      if(k != head){
        if(!spcache) L2th += nw(tail, k) && nw(head, k) && nw(k, head);
        if(htedge && nw(head, k) && nw(k, head))
          update_path(spcache ? GETUDMUI(k, tail, spcache) : count_rtp(nw, k, tail, tail, 0));
      }
    }

    for(auto k: nw.out_neighbors(tail))
      if(k != head && htedge && nw(head, k) && nw(k, head))
        update_path(spcache ? GETUDMUI(tail, k, spcache) : count_rtp(nw, k, tail, tail, 0));

    for(auto k: nw.in_neighbors(head))
      if(k != tail && htedge && nw(tail, k) && nw(k, tail))
        update_path(spcache ? GETUDMUI(k, head, spcache) : count_rtp(nw, k, head, head, 0));

    for(auto k: nw.out_neighbors(head))
      if(k != tail && htedge && nw(tail, k) && nw(k, tail))
        update_path(spcache ? GETUDMUI(head, k, spcache) : count_rtp(nw, k, head, head, 0));

    update_focus(L2th);
  }else if constexpr(type == L2OSP){
    int L2th = spcache ? GETUDMUI(tail, head, spcache) : 0;

    for(auto k: nw.out_neighbors(tail))
      if(k != head){
        if(!spcache) L2th += nw(head, k);
        if(nw(k, head))
          update_path(spcache ? GETUDMUI(tail, k, spcache) : count_osp(nw, tail, k));
      }

    for(auto k: nw.in_neighbors(tail))
      if(k != head && nw(k, head))
        update_path(spcache ? GETUDMUI(k, tail, spcache) : count_osp(nw, tail, k));

    update_focus(L2th);
  }else if constexpr(type == L2ISP){
    int L2th = spcache ? GETUDMUI(tail, head, spcache) : 0;

    for(auto k: nw.in_neighbors(head))
      if(k != tail){
        if(!spcache) L2th += nw(k, tail);
        if(nw(tail, k))
          update_path(spcache ? GETUDMUI(k, head, spcache) : count_isp(nw, head, k));
      }

    for(auto k: nw.out_neighbors(head))
      if(k != tail && nw(tail, k))
        update_path(spcache ? GETUDMUI(head, k, spcache) : count_isp(nw, head, k));

    update_focus(L2th);
  }
}

template<typename NetworkView, typename UpdatePath, typename UpdateFocus>
inline void esp_change(L2Type type, Vertex tail, Vertex head, NetworkView& nw, StoreStrictDyadMapUInt *spcache, UpdatePath update_path, UpdateFocus update_focus){
  switch(type){
  case L2UTP: esp_change<L2UTP>(tail, head, nw, spcache, update_path, update_focus); break;
  case L2OTP: esp_change<L2OTP>(tail, head, nw, spcache, update_path, update_focus); break;
  case L2ITP: esp_change<L2ITP>(tail, head, nw, spcache, update_path, update_focus); break;
  case L2RTP: esp_change<L2RTP>(tail, head, nw, spcache, update_path, update_focus); break;
  case L2OSP: esp_change<L2OSP>(tail, head, nw, spcache, update_path, update_focus); break;
  case L2ISP: esp_change<L2ISP>(tail, head, nw, spcache, update_path, update_focus); break;
  default: error("In ergm shared partner helper, an unsupported type of triad: %d.", type);
  }
}

template<typename NetworkView>
inline int dsp_nonzero_change(L2Type type, Vertex tail, Vertex head, NetworkView& nw, Rboolean edgestate, StoreStrictDyadMapUInt *spcache){
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
