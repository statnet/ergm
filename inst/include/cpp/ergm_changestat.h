/*  File inst/include/cpp/ergm_changestat.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#pragma once

#include "FixedArray.h"
#include "ergm_modelterm_base.h"
#include "ergm_network.h"

#include "../ergm_changestat.h"
#include "../ergm_variadic_macros.h"

namespace ergm {
inline namespace v1 {

template<typename StorageType = void>
  using ErgmCppModelTerm = ErgmCppModelTermBase<ModelTerm, StorageType>;

} // namespace v1
} // namespace ergm

#define _C_CHANGESTAT_CPP_2(name, impl)                                 \
  extern "C" void c_ ## name (Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, Rboolean edgestate) { \
    ergm::ErgmCppNetwork nw(nwp);                                       \
    ergm::ErgmCppModelTerm<> mt(mtp);                                   \
    impl;                                                               \
  }

#define _C_CHANGESTAT_CPP_3(name, StorageType, impl)                    \
  extern "C" void c_ ## name (Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, Rboolean edgestate) { \
    ergm::ErgmCppNetwork nw(nwp);                                       \
    ergm::ErgmCppModelTerm<StorageType> mt(mtp);                        \
    impl;                                                               \
  }

#define C_CHANGESTAT_CPP(...)                                           \
  _GET_OVERRIDE3(__VA_ARGS__, _C_CHANGESTAT_CPP_3, _C_CHANGESTAT_CPP_2, , )(__VA_ARGS__)

#define _S_CHANGESTAT_CPP_2(name, impl)                                 \
  extern "C" void s_ ## name (ModelTerm *mtp, Network *nwp) {           \
    ergm::ErgmCppNetwork nw(nwp);                                       \
    ergm::ErgmCppModelTerm<> mt(mtp);                                   \
    impl;                                                               \
  }

#define _S_CHANGESTAT_CPP_3(name, StorageType, impl)                    \
  extern "C" void s_ ## name (ModelTerm *mtp, Network *nwp) {           \
    ergm::ErgmCppNetwork nw(nwp);                                       \
    ergm::ErgmCppModelTerm<StorageType> mt(mtp);                        \
    impl;                                                               \
  }

#define S_CHANGESTAT_CPP(...)                                           \
  _GET_OVERRIDE3(__VA_ARGS__, _S_CHANGESTAT_CPP_3, _S_CHANGESTAT_CPP_2, , )(__VA_ARGS__)

#define _D_CHANGESTAT_CPP_2(name, impl)                                 \
  extern "C" void d_ ## name (Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp) { \
    ergm::ErgmCppNetwork nw(nwp);                                       \
    ergm::ErgmCppModelTerm<> mt(mtp);                                   \
    impl;                                                               \
  }

#define _D_CHANGESTAT_CPP_3(name, StorageType, impl)                    \
  extern "C" void d_ ## name (Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp) { \
    ergm::ErgmCppNetwork nw(nwp);                                       \
    ergm::ErgmCppModelTerm<StorageType> mt(mtp);                        \
    impl;                                                               \
  }

#define D_CHANGESTAT_CPP(...)                                           \
  _GET_OVERRIDE3(__VA_ARGS__, _D_CHANGESTAT_CPP_3, _D_CHANGESTAT_CPP_2, , )(__VA_ARGS__)

#define _I_CHANGESTAT_CPP_2(name, impl)                                 \
  extern "C" void i_ ## name (ModelTerm *mtp, Network *nwp) {           \
    ergm::ErgmCppNetwork nw(nwp);                                       \
    ergm::ErgmCppModelTerm<> mt(mtp);                                   \
    impl;                                                               \
  }

#define _I_CHANGESTAT_CPP_3(name, StorageType, impl)                    \
  extern "C" void i_ ## name (ModelTerm *mtp, Network *nwp) {           \
    ergm::ErgmCppNetwork nw(nwp);                                       \
    ergm::ErgmCppModelTerm<StorageType> mt(mtp);                        \
    impl;                                                               \
  }

#define I_CHANGESTAT_CPP(...)                                           \
  _GET_OVERRIDE3(__VA_ARGS__, _I_CHANGESTAT_CPP_3, _I_CHANGESTAT_CPP_2, , )(__VA_ARGS__)

#define _U_CHANGESTAT_CPP_2(name, impl)                                 \
  extern "C" void u_ ## name (Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, Rboolean edgestate) { \
    ergm::ErgmCppNetwork nw(nwp);                                       \
    ergm::ErgmCppModelTerm<> mt(mtp);                                   \
    impl;                                                               \
  }

#define _U_CHANGESTAT_CPP_3(name, StorageType, impl)                    \
  extern "C" void u_ ## name (Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, Rboolean edgestate) { \
    ergm::ErgmCppNetwork nw(nwp);                                       \
    ergm::ErgmCppModelTerm<StorageType> mt(mtp);                        \
    impl;                                                               \
  }

#define U_CHANGESTAT_CPP(...)                                           \
  _GET_OVERRIDE3(__VA_ARGS__, _U_CHANGESTAT_CPP_3, _U_CHANGESTAT_CPP_2, , )(__VA_ARGS__)

#define _F_CHANGESTAT_CPP_2(name, impl)                                 \
  extern "C" void f_ ## name (ModelTerm *mtp, Network *nwp) {           \
    ergm::ErgmCppNetwork nw(nwp);                                       \
    ergm::ErgmCppModelTerm<> mt(mtp);                                   \
    impl;                                                               \
  }

#define _F_CHANGESTAT_CPP_3(name, StorageType, impl)                    \
  extern "C" void f_ ## name (ModelTerm *mtp, Network *nwp) {           \
    ergm::ErgmCppNetwork nw(nwp);                                       \
    ergm::ErgmCppModelTerm<StorageType> mt(mtp);                        \
    impl;                                                               \
  }

#define F_CHANGESTAT_CPP(...)                                           \
  _GET_OVERRIDE3(__VA_ARGS__, _F_CHANGESTAT_CPP_3, _F_CHANGESTAT_CPP_2, , )(__VA_ARGS__)

#define _W_CHANGESTAT_CPP_2(name, impl)                                 \
  extern "C" SEXP w_ ## name (ModelTerm *mtp, Network *nwp) {           \
    ergm::ErgmCppNetwork nw(nwp);                                       \
    ergm::ErgmCppModelTerm<> mt(mtp);                                   \
    impl;                                                               \
  }

#define _W_CHANGESTAT_CPP_3(name, StorageType, impl)                    \
  extern "C" SEXP w_ ## name (ModelTerm *mtp, Network *nwp) {           \
    ergm::ErgmCppNetwork nw(nwp);                                       \
    ergm::ErgmCppModelTerm<StorageType> mt(mtp);                        \
    impl;                                                               \
  }

#define W_CHANGESTAT_CPP(...)                                           \
  _GET_OVERRIDE3(__VA_ARGS__, _W_CHANGESTAT_CPP_3, _W_CHANGESTAT_CPP_2, , )(__VA_ARGS__)

#define _X_CHANGESTAT_CPP_2(name, impl)                                 \
  extern "C" void x_ ## name (unsigned int type, void *data, ModelTerm *mtp, Network *nwp) { \
    ergm::ErgmCppNetwork nw(nwp);                                       \
    ergm::ErgmCppModelTerm<> mt(mtp);                                   \
    impl;                                                               \
  }

#define _X_CHANGESTAT_CPP_3(name, StorageType, impl)                    \
  extern "C" void x_ ## name (unsigned int type, void *data, ModelTerm *mtp, Network *nwp) { \
    ergm::ErgmCppNetwork nw(nwp);                                       \
    ergm::ErgmCppModelTerm<StorageType> mt(mtp);                        \
    impl;                                                               \
  }

#define X_CHANGESTAT_CPP(...)                                           \
  _GET_OVERRIDE3(__VA_ARGS__, _X_CHANGESTAT_CPP_3, _X_CHANGESTAT_CPP_2, , )(__VA_ARGS__)

#define _Z_CHANGESTAT_CPP_2(name, impl)                                 \
  extern "C" void z_ ## name (ModelTerm *mtp, Network *nwp, Rboolean skip_s) { \
    ergm::ErgmCppNetwork nw(nwp);                                       \
    ergm::ErgmCppModelTerm<> mt(mtp);                                   \
    impl;                                                               \
  }

#define _Z_CHANGESTAT_CPP_3(name, StorageType, impl)                    \
  extern "C" void z_ ## name (ModelTerm *mtp, Network *nwp, Rboolean skip_s) { \
    ergm::ErgmCppNetwork nw(nwp);                                       \
    ergm::ErgmCppModelTerm<StorageType> mt(mtp);                        \
    impl;                                                               \
  }

#define Z_CHANGESTAT_CPP(...)                                           \
  _GET_OVERRIDE3(__VA_ARGS__, _Z_CHANGESTAT_CPP_3, _Z_CHANGESTAT_CPP_2, , )(__VA_ARGS__)
