/*  File src/changestats_dgw_sp.cpp in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2026 Statnet Commons
 */
#include <cmath>

#include "cpp/ergm_changestat.h"

#include "changestats_dgw_sp.h"
#include "changestats.h"

using ergm::ErgmCppModelTerm;
using ergm::ErgmCppNetwork;

namespace {

inline StoreStrictDyadMapUInt* get_spcache(const ErgmCppModelTerm<>& mt){
  return mt.aux_storage.size() ? static_cast<StoreStrictDyadMapUInt*>(mt.aux_storage[0]) : nullptr;
}

inline int dsp_path_multiplier(L2Type type){
  switch(type){
  case L2OSP:
  case L2ISP:
  case L2RTP:
    return 2;
  default:
    return 1;
  }
}

inline void negate_change_stats(ErgmCppModelTerm<>& mt){
  for(double& stat : mt.stat) stat *= -1.0;
}

template<int PathMultiplier>
inline void dsp_vector_change(L2Type type, Vertex tail, Vertex head, ErgmCppModelTerm<>& mt, ErgmCppNetwork& nw, Rboolean edgestate, StoreStrictDyadMapUInt *spcache){
  int echange = edgestate ? -1 : 1;
  ergm::sp::dsp_change(type, tail, head, nw, spcache,
                       [&](int L2){
                         for(std::size_t j = 0; j < mt.stat.size(); ++j){
                           int deg = mt.iinput[j+1];
                           mt.stat[j] += (((L2 + echange) == deg) - (L2 == deg)) * PathMultiplier;
                         }
                       },
                       [&](int){});
}

inline void esp_vector_change(L2Type type, Vertex tail, Vertex head, ErgmCppModelTerm<>& mt, ErgmCppNetwork& nw, Rboolean edgestate, StoreStrictDyadMapUInt *spcache){
  int echange = edgestate ? -1 : 1;
  ergm::sp::esp_change(type, tail, head, nw, spcache,
                       [&](int L2){
                         for(std::size_t j = 0; j < mt.stat.size(); ++j){
                           int deg = mt.iinput[j+1];
                           mt.stat[j] += ((L2 + echange) == deg) - (L2 == deg);
                         }
                       },
                       [&](int L2){
                         for(std::size_t j = 0; j < mt.stat.size(); ++j){
                           int deg = mt.iinput[j+1];
                           mt.stat[j] += echange * (L2 == deg);
                         }
                       });
}

template<int PathMultiplier>
inline void dsp_dist_change(L2Type type, Vertex tail, Vertex head, ModelTerm *mtp, ErgmCppModelTerm<>& mt, ErgmCppNetwork& nw, Rboolean edgestate, StoreStrictDyadMapUInt *spcache){
  int echange = edgestate ? -1 : 1;
  int nd = static_cast<int>(mt.stat.size());
  ergm::sp::dsp_change(type, tail, head, nw, spcache,
                       [&](int L2){
                         int nL2 = L2 + echange;
                         if(nL2 > nd) cutoff_error(mtp);
                         if(L2) mt.stat[L2-1] -= PathMultiplier;
                         if(nL2) mt.stat[nL2-1] += PathMultiplier;
                       },
                       [&](int){});
}

inline void esp_dist_change(L2Type type, Vertex tail, Vertex head, ModelTerm *mtp, ErgmCppModelTerm<>& mt, ErgmCppNetwork& nw, Rboolean edgestate, StoreStrictDyadMapUInt *spcache){
  int echange = edgestate ? -1 : 1;
  int nd = static_cast<int>(mt.stat.size());
  ergm::sp::esp_change(type, tail, head, nw, spcache,
                       [&](int L2){
                         int nL2 = L2 + echange;
                         if(nL2 > nd) cutoff_error(mtp);
                         if(L2) mt.stat[L2-1]--;
                         if(nL2) mt.stat[nL2-1]++;
                       },
                       [&](int L2){
                         if(L2 > nd) cutoff_error(mtp);
                         if(L2) mt.stat[L2-1] += echange;
                       });
}

template<int PathMultiplier>
inline double dsp_gw_change(L2Type type, Vertex tail, Vertex head, ErgmCppNetwork& nw, Rboolean edgestate, StoreStrictDyadMapUInt *spcache, double alpha, double loneexpa){
  double cumchange = 0;
  ergm::sp::dsp_change(type, tail, head, nw, spcache,
                       [&](int L2){
                         cumchange += (alpha ? exp(loneexpa*(L2-edgestate)) : L2-edgestate == 0) * PathMultiplier;
                       },
                       [&](int){});
  return cumchange;
}

inline double esp_gw_change(L2Type type, Vertex tail, Vertex head, ErgmCppNetwork& nw, Rboolean edgestate, StoreStrictDyadMapUInt *spcache, double alpha, double loneexpa){
  double cumchange = 0;
  ergm::sp::esp_change(type, tail, head, nw, spcache,
                       [&](int L2){
                         cumchange += alpha ? exp(loneexpa*(L2-edgestate)) : L2-edgestate == 0;
                       },
                       [&](int L2){
                         cumchange += alpha ? exp(alpha + log1mexp(-loneexpa*L2)) : L2 != 0;
                       });
  return cumchange;
}

} // namespace

/*****************
 changestat: d_dsp
*****************/
C_CHANGESTAT_CPP(ddsp, {
    auto *spcache = get_spcache(mt);
    L2Type type = static_cast<L2Type>(mt.iinput[0]);

    if(dsp_path_multiplier(type) == 1) dsp_vector_change<1>(type, tail, head, mt, nw, edgestate, spcache);
    else dsp_vector_change<2>(type, tail, head, mt, nw, edgestate, spcache);
  })

C_CHANGESTAT_CPP(ddspdist, {
    auto *spcache = get_spcache(mt);
    L2Type type = static_cast<L2Type>(mt.iinput[0]);

    if(dsp_path_multiplier(type) == 1) dsp_dist_change<1>(type, tail, head, mtp, mt, nw, edgestate, spcache);
    else dsp_dist_change<2>(type, tail, head, mtp, mt, nw, edgestate, spcache);
  })

/*****************
 changestat: d_gwdsp
*****************/
C_CHANGESTAT_CPP(dgwdsp, {
    auto *spcache = get_spcache(mt);
    double alpha = mt.dinput[0];
    double loneexpa = log1mexp(alpha);
    L2Type type = static_cast<L2Type>(mt.iinput[0]);

    double cumchange = dsp_path_multiplier(type) == 1
      ? dsp_gw_change<1>(type, tail, head, nw, edgestate, spcache, alpha, loneexpa)
      : dsp_gw_change<2>(type, tail, head, nw, edgestate, spcache, alpha, loneexpa);

    mt.stat[0] = edgestate ? -cumchange : cumchange;
  })

/*****************
 changestat: d_esp
*****************/
C_CHANGESTAT_CPP(desp, {
    auto *spcache = get_spcache(mt);
    L2Type type = static_cast<L2Type>(mt.iinput[0]);

    esp_vector_change(type, tail, head, mt, nw, edgestate, spcache);
  })

C_CHANGESTAT_CPP(despdist, {
    auto *spcache = get_spcache(mt);
    L2Type type = static_cast<L2Type>(mt.iinput[0]);

    esp_dist_change(type, tail, head, mtp, mt, nw, edgestate, spcache);
  })

/*****************
 changestat: d_gwesp
*****************/
C_CHANGESTAT_CPP(dgwesp, {
    auto *spcache = get_spcache(mt);
    double alpha = mt.dinput[0];
    double loneexpa = log1mexp(alpha);
    L2Type type = static_cast<L2Type>(mt.iinput[0]);

    double cumchange = esp_gw_change(type, tail, head, nw, edgestate, spcache, alpha, loneexpa);
    mt.stat[0] = edgestate ? -cumchange : cumchange;
  })

/*****************
 changestat: d_nsp
*****************/
C_CHANGESTAT_CPP(dnsp, {
    auto *spcache = get_spcache(mt);
    L2Type type = static_cast<L2Type>(mt.iinput[0]);

    esp_vector_change(type, tail, head, mt, nw, edgestate, spcache);
    negate_change_stats(mt);
    if(dsp_path_multiplier(type) == 1) dsp_vector_change<1>(type, tail, head, mt, nw, edgestate, spcache);
    else dsp_vector_change<2>(type, tail, head, mt, nw, edgestate, spcache);
  })

C_CHANGESTAT_CPP(dnspdist, {
    auto *spcache = get_spcache(mt);
    L2Type type = static_cast<L2Type>(mt.iinput[0]);

    esp_dist_change(type, tail, head, mtp, mt, nw, edgestate, spcache);
    negate_change_stats(mt);
    if(dsp_path_multiplier(type) == 1) dsp_dist_change<1>(type, tail, head, mtp, mt, nw, edgestate, spcache);
    else dsp_dist_change<2>(type, tail, head, mtp, mt, nw, edgestate, spcache);
  })

/*****************
 changestat: d_gwnsp
*****************/
C_CHANGESTAT_CPP(dgwnsp, {
    auto *spcache = get_spcache(mt);
    double alpha = mt.dinput[0];
    double loneexpa = log1mexp(alpha);
    L2Type type = static_cast<L2Type>(mt.iinput[0]);

    double dspchange = dsp_path_multiplier(type) == 1
      ? dsp_gw_change<1>(type, tail, head, nw, edgestate, spcache, alpha, loneexpa)
      : dsp_gw_change<2>(type, tail, head, nw, edgestate, spcache, alpha, loneexpa);
    double cumchange = dspchange - esp_gw_change(type, tail, head, nw, edgestate, spcache, alpha, loneexpa);

    mt.stat[0] = edgestate ? -cumchange : cumchange;
  })

/*****************
 changestat: c_ddspbwrap
*****************/
C_CHANGESTAT_CPP(ddspbwrap, {
    c_ddsp(tail, head, mtp, nwp, edgestate);
    for(double &stat : mt.stat) stat /= 2.0;
  })

C_CHANGESTAT_CPP(ddspdistbwrap, {
    c_ddspdist(tail, head, mtp, nwp, edgestate);
    for(double &stat : mt.stat) stat /= 2.0;
  })

/*****************
 changestat: c_dgwdspbwrap
*****************/
C_CHANGESTAT_CPP(dgwdspbwrap, {
    c_dgwdsp(tail, head, mtp, nwp, edgestate);
    mt.stat[0] /= 2.0;
  })
