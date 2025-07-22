
#pragma once

#include "FixedArray.h"
#include "ergm_modelterm_base.h"

extern "C" {
#include "../ergm_wtchangestat.h"
}

#define WtC_CHANGESTAT_CPP(name, impl, ...)                                    \
  extern "C" void c_ ## name (Vertex tail, Vertex head, double weight, WtModelTerm *mtp, WtNetwork *nwp, double edgestate) { \
    ErgmCppWtNetwork nw(nwp);                                             \
    ErgmCppWtModelTerm<__VA_ARGS__> mt(mtp);                                           \
    impl;                                                             \
  }

#define WtS_CHANGESTAT_CPP(name, impl, ...)                                    \
  extern "C" void s_ ## name (WtModelTerm *mtp, WtNetwork *nwp) { \
    ErgmCppWtNetwork nw(nwp);                                             \
    ErgmCppWtModelTerm<__VA_ARGS__> mt(mtp);                                           \
    impl;                                                             \
  }

template<typename StorageType = void>
using ErgmCppWtModelTerm = ErgmCppModelTermBase<WtModelTerm, StorageType>;
