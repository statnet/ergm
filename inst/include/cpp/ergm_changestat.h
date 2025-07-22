#pragma once


#pragma once

#include "FixedArray.h"
#include "ergm_modelterm_base.h"

extern "C" {
#include "../ergm_changestat.h"
}

#define C_CHANGESTAT_CPP(name, impl, ...)                                    \
  extern "C" void c_ ## name (Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, int edgestate) { \
    ErgmCppNetwork nw(nwp);                                             \
    ErgmCppModelTerm<__VA_ARGS__> mt(mtp);                                           \
    impl;                                                            \
  }

#define S_CHANGESTAT_CPP(name, impl, ...)                                    \
  extern "C" void s_ ## name (ModelTerm *mtp, Network *nwp) { \
    ErgmCppNetwork nw(nwp);                                             \
    ErgmCppModelTerm<__VA_ARGS__> mt(mtp);                                           \
    impl;                                                             \
  }

template<typename StorageType = void>
using ErgmCppModelTerm = ErgmCppModelTermBase<ModelTerm, StorageType>;
