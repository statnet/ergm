#pragma once

#include "FixedArray.h"
#include "ergm_modelterm_base.h"
#include "ergm_network.h"

extern "C" {
#include "../ergm_changestat.h"
}

#define C_CHANGESTAT_CPP(name, impl, ...)                               \
  extern "C" void c_ ## name (Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, int edgestate) { \
    ErgmCppNetwork nw(nwp);                                             \
    ErgmCppModelTerm<__VA_ARGS__> mt(mtp);                              \
    impl;                                                               \
  }

#define S_CHANGESTAT_CPP(name, impl, ...)                       \
  extern "C" void s_ ## name (ModelTerm *mtp, Network *nwp) {   \
    ErgmCppNetwork nw(nwp);                                     \
    ErgmCppModelTerm<__VA_ARGS__> mt(mtp);                      \
    impl;                                                       \
  }

#define D_CHANGESTAT_CPP(name, impl, ...)                               \
  extern "C" void d_ ## name (Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp) { \
    ErgmCppNetwork nw(nwp);                                             \
    ErgmCppModelTerm<__VA_ARGS__> mt(mtp);                              \
    impl;                                                               \
  }

#define I_CHANGESTAT_CPP(name, impl, ...)                       \
  extern "C" void i_ ## name (ModelTerm *mtp, Network *nwp) {   \
    ErgmCppNetwork nw(nwp);                                     \
    ErgmCppModelTerm<__VA_ARGS__> mt(mtp);                      \
    impl;                                                       \
  }

#define U_CHANGESTAT_CPP(name, impl, ...)                               \
  extern "C" void u_ ## name (Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, Rboolean edgestate) { \
    ErgmCppNetwork nw(nwp);                                             \
    ErgmCppModelTerm<__VA_ARGS__> mt(mtp);                              \
    impl;                                                               \
  }

#define F_CHANGESTAT_CPP(name, impl, ...)                       \
  extern "C" void f_ ## name (ModelTerm *mtp, Network *nwp) {   \
    ErgmCppNetwork nw(nwp);                                     \
    ErgmCppModelTerm<__VA_ARGS__> mt(mtp);                      \
    impl;                                                       \
  }

#define W_CHANGESTAT_CPP(name, impl, ...)                       \
  extern "C" SEXP w_ ## name (ModelTerm *mtp, Network *nwp) {   \
    ErgmCppNetwork nw(nwp);                                     \
    ErgmCppModelTerm<__VA_ARGS__> mt(mtp);                      \
    impl;                                                       \
  }

#define X_CHANGESTAT_CPP(name, impl, ...)                               \
  extern "C" void x_ ## name (unsigned int type, void *data, ModelTerm *mtp, Network *nwp) { \
    ErgmCppNetwork nw(nwp);                                             \
    ErgmCppModelTerm<__VA_ARGS__> mt(mtp);                              \
    impl;                                                               \
  }

#define Z_CHANGESTAT_CPP(name, impl, ...)                               \
  extern "C" void z_ ## name (ModelTerm *mtp, Network *nwp, Rboolean skip_s) { \
    ErgmCppNetwork nw(nwp);                                             \
    ErgmCppModelTerm<__VA_ARGS__> mt(mtp);                              \
    impl;                                                               \
  }

template<typename StorageType = void>
  using ErgmCppModelTerm = ErgmCppModelTermBase<ModelTerm, StorageType>;
