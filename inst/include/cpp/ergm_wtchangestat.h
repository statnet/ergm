#pragma once

#include "FixedArray.h"
#include "ergm_modelterm_base.h"
#include "ergm_wtnetwork.h"

#include "../ergm_wtchangestat.h"

namespace ergm {
inline namespace v1 {

template<typename StorageType = void>
  using ErgmCppWtModelTerm = ErgmCppModelTermBase<WtModelTerm, StorageType>;

} // namespace v1
} // namespace ergm

#define WtC_CHANGESTAT_CPP(name, impl, ...)                             \
  extern "C" void c_ ## name (Vertex tail, Vertex head, double weight, WtModelTerm *mtp, WtNetwork *nwp, double edgestate) { \
    ergm::ErgmCppWtNetwork nw(nwp);                                     \
    ergm::ErgmCppWtModelTerm<__VA_ARGS__> mt(mtp);                      \
    impl;                                                               \
  }

#define WtS_CHANGESTAT_CPP(name, impl, ...)                             \
  extern "C" void s_ ## name (WtModelTerm *mtp, WtNetwork *nwp) {       \
    ergm::ErgmCppWtNetwork nw(nwp);                                     \
    ergm::ErgmCppWtModelTerm<__VA_ARGS__> mt(mtp);                      \
    impl;                                                               \
  }

#define WtD_CHANGESTAT_CPP(name, impl, ...)                             \
  extern "C" void d_ ## name (Edge ntoggles, Vertex *tails, Vertex *heads, double *weights, WtModelTerm *mtp, WtNetwork *nwp) { \
    ergm::ErgmCppWtNetwork nw(nwp);                                     \
    ergm::ErgmCppWtModelTerm<__VA_ARGS__> mt(mtp);                      \
    impl;                                                               \
  }

#define WtI_CHANGESTAT_CPP(name, impl, ...)                             \
  extern "C" void i_ ## name (WtModelTerm *mtp, WtNetwork *nwp) {       \
    ergm::ErgmCppWtNetwork nw(nwp);                                     \
    ergm::ErgmCppWtModelTerm<__VA_ARGS__> mt(mtp);                      \
    impl;                                                               \
  }

#define WtU_CHANGESTAT_CPP(name, impl, ...)                             \
  extern "C" void u_ ## name (Vertex tail, Vertex head, double weight, WtModelTerm *mtp, WtNetwork *nwp, double edgestate) { \
    ergm::ErgmCppWtNetwork nw(nwp);                                     \
    ergm::ErgmCppWtModelTerm<__VA_ARGS__> mt(mtp);                      \
    impl;                                                               \
  }

#define WtF_CHANGESTAT_CPP(name, impl, ...)                             \
  extern "C" void f_ ## name (WtModelTerm *mtp, WtNetwork *nwp) {       \
    ergm::ErgmCppWtNetwork nw(nwp);                                     \
    ergm::ErgmCppWtModelTerm<__VA_ARGS__> mt(mtp);                      \
    impl;                                                               \
  }

#define WtW_CHANGESTAT_CPP(name, impl, ...)                             \
  extern "C" SEXP w_ ## name (WtModelTerm *mtp, WtNetwork *nwp) {       \
    ergm::ErgmCppWtNetwork nw(nwp);                                     \
    ergm::ErgmCppWtModelTerm<__VA_ARGS__> mt(mtp);                      \
    impl;                                                               \
  }

#define WtX_CHANGESTAT_CPP(name, impl, ...)                             \
  extern "C" void x_ ## name (unsigned int type, void *data, WtModelTerm *mtp, WtNetwork *nwp) { \
    ergm::ErgmCppWtNetwork nw(nwp);                                     \
    ergm::ErgmCppWtModelTerm<__VA_ARGS__> mt(mtp);                      \
    impl;                                                               \
  }

#define WtZ_CHANGESTAT_CPP(name, impl, ...)                             \
  extern "C" void z_ ## name (WtModelTerm *mtp, WtNetwork *nwp, Rboolean skip_s) { \
    ergm::ErgmCppWtNetwork nw(nwp);                                     \
    ergm::ErgmCppWtModelTerm<__VA_ARGS__> mt(mtp);                      \
    impl;                                                               \
  }
