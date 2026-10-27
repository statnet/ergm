/*  File inst/include/cpp/ergm_wtchangestat.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2026 Statnet Commons
 */
#pragma once

#include "FixedArray.h"
#include "ergm_modelterm_base.h"
#include "ergm_wtnetwork.h"

#include "../ergm_wtchangestat.h"
#include "../ergm_variadic_macros.h"

namespace ergm {
inline namespace v1 {

template<typename StorageType = void>
  using ErgmCppWtModelTerm = ErgmCppModelTermBase<WtModelTerm, StorageType>;

} // namespace v1
} // namespace ergm

#define _WtC_CHANGESTAT_CPP_2(name, impl)                               \
  extern "C" void c_ ## name (Vertex tail, Vertex head, double weight, WtModelTerm *mtp, WtNetwork *nwp, double edgestate) { \
    ergm::ErgmCppWtNetwork nw(nwp);                                     \
    ergm::ErgmCppWtModelTerm<> mt(mtp);                                 \
    impl;                                                               \
  }

#define _WtC_CHANGESTAT_CPP_3(name, StorageType, impl)                  \
  extern "C" void c_ ## name (Vertex tail, Vertex head, double weight, WtModelTerm *mtp, WtNetwork *nwp, double edgestate) { \
    ergm::ErgmCppWtNetwork nw(nwp);                                     \
    ergm::ErgmCppWtModelTerm<StorageType> mt(mtp);                      \
    impl;                                                               \
  }

#define WtC_CHANGESTAT_CPP(...)                                         \
  _GET_OVERRIDE3(__VA_ARGS__, _WtC_CHANGESTAT_CPP_3, _WtC_CHANGESTAT_CPP_2, , )(__VA_ARGS__)

#define _WtS_CHANGESTAT_CPP_2(name, impl)                               \
  extern "C" void s_ ## name (WtModelTerm *mtp, WtNetwork *nwp) {       \
    ergm::ErgmCppWtNetwork nw(nwp);                                     \
    ergm::ErgmCppWtModelTerm<> mt(mtp);                                 \
    impl;                                                               \
  }

#define _WtS_CHANGESTAT_CPP_3(name, StorageType, impl)                  \
  extern "C" void s_ ## name (WtModelTerm *mtp, WtNetwork *nwp) {       \
    ergm::ErgmCppWtNetwork nw(nwp);                                     \
    ergm::ErgmCppWtModelTerm<StorageType> mt(mtp);                      \
    impl;                                                               \
  }

#define WtS_CHANGESTAT_CPP(...)                                         \
  _GET_OVERRIDE3(__VA_ARGS__, _WtS_CHANGESTAT_CPP_3, _WtS_CHANGESTAT_CPP_2, , )(__VA_ARGS__)

#define _WtD_CHANGESTAT_CPP_2(name, impl)                               \
  extern "C" void d_ ## name (Edge ntoggles, Vertex *tails, Vertex *heads, double *weights, WtModelTerm *mtp, WtNetwork *nwp) { \
    ergm::ErgmCppWtNetwork nw(nwp);                                     \
    ergm::ErgmCppWtModelTerm<> mt(mtp);                                 \
    impl;                                                               \
  }

#define _WtD_CHANGESTAT_CPP_3(name, StorageType, impl)                  \
  extern "C" void d_ ## name (Edge ntoggles, Vertex *tails, Vertex *heads, double *weights, WtModelTerm *mtp, WtNetwork *nwp) { \
    ergm::ErgmCppWtNetwork nw(nwp);                                     \
    ergm::ErgmCppWtModelTerm<StorageType> mt(mtp);                      \
    impl;                                                               \
  }

#define WtD_CHANGESTAT_CPP(...)                                         \
  _GET_OVERRIDE3(__VA_ARGS__, _WtD_CHANGESTAT_CPP_3, _WtD_CHANGESTAT_CPP_2, , )(__VA_ARGS__)

#define _WtI_CHANGESTAT_CPP_2(name, impl)                               \
  extern "C" void i_ ## name (WtModelTerm *mtp, WtNetwork *nwp) {       \
    ergm::ErgmCppWtNetwork nw(nwp);                                     \
    ergm::ErgmCppWtModelTerm<> mt(mtp);                                 \
    impl;                                                               \
  }

#define _WtI_CHANGESTAT_CPP_3(name, StorageType, impl)                  \
  extern "C" void i_ ## name (WtModelTerm *mtp, WtNetwork *nwp) {       \
    ergm::ErgmCppWtNetwork nw(nwp);                                     \
    ergm::ErgmCppWtModelTerm<StorageType> mt(mtp);                      \
    impl;                                                               \
  }

#define WtI_CHANGESTAT_CPP(...)                                         \
  _GET_OVERRIDE3(__VA_ARGS__, _WtI_CHANGESTAT_CPP_3, _WtI_CHANGESTAT_CPP_2, , )(__VA_ARGS__)

#define _WtU_CHANGESTAT_CPP_2(name, impl)                               \
  extern "C" void u_ ## name (Vertex tail, Vertex head, double weight, WtModelTerm *mtp, WtNetwork *nwp, double edgestate) { \
    ergm::ErgmCppWtNetwork nw(nwp);                                     \
    ergm::ErgmCppWtModelTerm<> mt(mtp);                                 \
    impl;                                                               \
  }

#define _WtU_CHANGESTAT_CPP_3(name, StorageType, impl)                  \
  extern "C" void u_ ## name (Vertex tail, Vertex head, double weight, WtModelTerm *mtp, WtNetwork *nwp, double edgestate) { \
    ergm::ErgmCppWtNetwork nw(nwp);                                     \
    ergm::ErgmCppWtModelTerm<StorageType> mt(mtp);                      \
    impl;                                                               \
  }

#define WtU_CHANGESTAT_CPP(...)                                         \
  _GET_OVERRIDE3(__VA_ARGS__, _WtU_CHANGESTAT_CPP_3, _WtU_CHANGESTAT_CPP_2, , )(__VA_ARGS__)

#define _WtF_CHANGESTAT_CPP_2(name, impl)                               \
  extern "C" void f_ ## name (WtModelTerm *mtp, WtNetwork *nwp) {       \
    ergm::ErgmCppWtNetwork nw(nwp);                                     \
    ergm::ErgmCppWtModelTerm<> mt(mtp);                                 \
    impl;                                                               \
  }

#define _WtF_CHANGESTAT_CPP_3(name, StorageType, impl)                  \
  extern "C" void f_ ## name (WtModelTerm *mtp, WtNetwork *nwp) {       \
    ergm::ErgmCppWtNetwork nw(nwp);                                     \
    ergm::ErgmCppWtModelTerm<StorageType> mt(mtp);                      \
    impl;                                                               \
  }

#define WtF_CHANGESTAT_CPP(...)                                         \
  _GET_OVERRIDE3(__VA_ARGS__, _WtF_CHANGESTAT_CPP_3, _WtF_CHANGESTAT_CPP_2, , )(__VA_ARGS__)

#define _WtW_CHANGESTAT_CPP_2(name, impl)                               \
  extern "C" SEXP w_ ## name (WtModelTerm *mtp, WtNetwork *nwp) {       \
    ergm::ErgmCppWtNetwork nw(nwp);                                     \
    ergm::ErgmCppWtModelTerm<> mt(mtp);                                 \
    impl;                                                               \
  }

#define _WtW_CHANGESTAT_CPP_3(name, StorageType, impl)                  \
  extern "C" SEXP w_ ## name (WtModelTerm *mtp, WtNetwork *nwp) {       \
    ergm::ErgmCppWtNetwork nw(nwp);                                     \
    ergm::ErgmCppWtModelTerm<StorageType> mt(mtp);                      \
    impl;                                                               \
  }

#define WtW_CHANGESTAT_CPP(...)                                         \
  _GET_OVERRIDE3(__VA_ARGS__, _WtW_CHANGESTAT_CPP_3, _WtW_CHANGESTAT_CPP_2, , )(__VA_ARGS__)

#define _WtX_CHANGESTAT_CPP_2(name, impl)                               \
  extern "C" void x_ ## name (unsigned int type, void *data, WtModelTerm *mtp, WtNetwork *nwp) { \
    ergm::ErgmCppWtNetwork nw(nwp);                                     \
    ergm::ErgmCppWtModelTerm<> mt(mtp);                                 \
    impl;                                                               \
  }

#define _WtX_CHANGESTAT_CPP_3(name, StorageType, impl)                  \
  extern "C" void x_ ## name (unsigned int type, void *data, WtModelTerm *mtp, WtNetwork *nwp) { \
    ergm::ErgmCppWtNetwork nw(nwp);                                     \
    ergm::ErgmCppWtModelTerm<StorageType> mt(mtp);                      \
    impl;                                                               \
  }

#define WtX_CHANGESTAT_CPP(...)                                         \
  _GET_OVERRIDE3(__VA_ARGS__, _WtX_CHANGESTAT_CPP_3, _WtX_CHANGESTAT_CPP_2, , )(__VA_ARGS__)
#define _WtZ_CHANGESTAT_CPP_2(name, impl)                               \
  extern "C" void z_ ## name (WtModelTerm *mtp, WtNetwork *nwp, Rboolean skip_s) { \
    ergm::ErgmCppWtNetwork nw(nwp);                                     \
    ergm::ErgmCppWtModelTerm<> mt(mtp);                                 \
    impl;                                                               \
  }

#define _WtZ_CHANGESTAT_CPP_3(name, StorageType, impl)                  \
  extern "C" void z_ ## name (WtModelTerm *mtp, WtNetwork *nwp, Rboolean skip_s) { \
    ergm::ErgmCppWtNetwork nw(nwp);                                     \
    ergm::ErgmCppWtModelTerm<StorageType> mt(mtp);                      \
    impl;                                                               \
  }

#define WtZ_CHANGESTAT_CPP(...)                                         \
  _GET_OVERRIDE3(__VA_ARGS__, _WtZ_CHANGESTAT_CPP_3, _WtZ_CHANGESTAT_CPP_2, , )(__VA_ARGS__)


/* The following are variants for when the statistic is not
   ambidirected, which can significantly speed up neigbor iteration by
   removing some code paths at compile time. ErgmCppNetwork() does
   check that nwp's directedness matches, so at worst, these will crash if the wrong directedness is given. */

// Directed

#define _WtC_CHANGESTAT_D_CPP_2(name, impl)                               \
  extern "C" void c_ ## name (Vertex tail, Vertex head, double weight, WtModelTerm *mtp, WtNetwork *nwp, double edgestate) { \
    ergm::ErgmCppWtNetwork nw(nwp, true);                                     \
    ergm::ErgmCppWtModelTerm<> mt(mtp);                                 \
    impl;                                                               \
  }

#define _WtC_CHANGESTAT_D_CPP_3(name, StorageType, impl)                  \
  extern "C" void c_ ## name (Vertex tail, Vertex head, double weight, WtModelTerm *mtp, WtNetwork *nwp, double edgestate) { \
    ergm::ErgmCppWtNetwork nw(nwp, true);                                     \
    ergm::ErgmCppWtModelTerm<StorageType> mt(mtp);                      \
    impl;                                                               \
  }

#define WtC_CHANGESTAT_D_CPP(...)                                         \
  _GET_OVERRIDE3(__VA_ARGS__, _WtC_CHANGESTAT_D_CPP_3, _WtC_CHANGESTAT_D_CPP_2, , )(__VA_ARGS__)

#define _WtS_CHANGESTAT_D_CPP_2(name, impl)                               \
  extern "C" void s_ ## name (WtModelTerm *mtp, WtNetwork *nwp) {       \
    ergm::ErgmCppWtNetwork nw(nwp, true);                                     \
    ergm::ErgmCppWtModelTerm<> mt(mtp);                                 \
    impl;                                                               \
  }

#define _WtS_CHANGESTAT_D_CPP_3(name, StorageType, impl)                  \
  extern "C" void s_ ## name (WtModelTerm *mtp, WtNetwork *nwp) {       \
    ergm::ErgmCppWtNetwork nw(nwp, true);                                     \
    ergm::ErgmCppWtModelTerm<StorageType> mt(mtp);                      \
    impl;                                                               \
  }

#define WtS_CHANGESTAT_D_CPP(...)                                         \
  _GET_OVERRIDE3(__VA_ARGS__, _WtS_CHANGESTAT_D_CPP_3, _WtS_CHANGESTAT_D_CPP_2, , )(__VA_ARGS__)

#define _WtD_CHANGESTAT_D_CPP_2(name, impl)                               \
  extern "C" void d_ ## name (Edge ntoggles, Vertex *tails, Vertex *heads, double *weights, WtModelTerm *mtp, WtNetwork *nwp) { \
    ergm::ErgmCppWtNetwork nw(nwp, true);                                     \
    ergm::ErgmCppWtModelTerm<> mt(mtp);                                 \
    impl;                                                               \
  }

#define _WtD_CHANGESTAT_D_CPP_3(name, StorageType, impl)                  \
  extern "C" void d_ ## name (Edge ntoggles, Vertex *tails, Vertex *heads, double *weights, WtModelTerm *mtp, WtNetwork *nwp) { \
    ergm::ErgmCppWtNetwork nw(nwp, true);                                     \
    ergm::ErgmCppWtModelTerm<StorageType> mt(mtp);                      \
    impl;                                                               \
  }

#define WtD_CHANGESTAT_D_CPP(...)                                         \
  _GET_OVERRIDE3(__VA_ARGS__, _WtD_CHANGESTAT_D_CPP_3, _WtD_CHANGESTAT_D_CPP_2, , )(__VA_ARGS__)

#define _WtI_CHANGESTAT_D_CPP_2(name, impl)                               \
  extern "C" void i_ ## name (WtModelTerm *mtp, WtNetwork *nwp) {       \
    ergm::ErgmCppWtNetwork nw(nwp, true);                                     \
    ergm::ErgmCppWtModelTerm<> mt(mtp);                                 \
    impl;                                                               \
  }

#define _WtI_CHANGESTAT_D_CPP_3(name, StorageType, impl)                  \
  extern "C" void i_ ## name (WtModelTerm *mtp, WtNetwork *nwp) {       \
    ergm::ErgmCppWtNetwork nw(nwp, true);                                     \
    ergm::ErgmCppWtModelTerm<StorageType> mt(mtp);                      \
    impl;                                                               \
  }

#define WtI_CHANGESTAT_D_CPP(...)                                         \
  _GET_OVERRIDE3(__VA_ARGS__, _WtI_CHANGESTAT_D_CPP_3, _WtI_CHANGESTAT_D_CPP_2, , )(__VA_ARGS__)

#define _WtU_CHANGESTAT_D_CPP_2(name, impl)                               \
  extern "C" void u_ ## name (Vertex tail, Vertex head, double weight, WtModelTerm *mtp, WtNetwork *nwp, double edgestate) { \
    ergm::ErgmCppWtNetwork nw(nwp, true);                                     \
    ergm::ErgmCppWtModelTerm<> mt(mtp);                                 \
    impl;                                                               \
  }

#define _WtU_CHANGESTAT_D_CPP_3(name, StorageType, impl)                  \
  extern "C" void u_ ## name (Vertex tail, Vertex head, double weight, WtModelTerm *mtp, WtNetwork *nwp, double edgestate) { \
    ergm::ErgmCppWtNetwork nw(nwp, true);                                     \
    ergm::ErgmCppWtModelTerm<StorageType> mt(mtp);                      \
    impl;                                                               \
  }

#define WtU_CHANGESTAT_D_CPP(...)                                         \
  _GET_OVERRIDE3(__VA_ARGS__, _WtU_CHANGESTAT_D_CPP_3, _WtU_CHANGESTAT_D_CPP_2, , )(__VA_ARGS__)

#define _WtF_CHANGESTAT_D_CPP_2(name, impl)                               \
  extern "C" void f_ ## name (WtModelTerm *mtp, WtNetwork *nwp) {       \
    ergm::ErgmCppWtNetwork nw(nwp, true);                                     \
    ergm::ErgmCppWtModelTerm<> mt(mtp);                                 \
    impl;                                                               \
  }

#define _WtF_CHANGESTAT_D_CPP_3(name, StorageType, impl)                  \
  extern "C" void f_ ## name (WtModelTerm *mtp, WtNetwork *nwp) {       \
    ergm::ErgmCppWtNetwork nw(nwp, true);                                     \
    ergm::ErgmCppWtModelTerm<StorageType> mt(mtp);                      \
    impl;                                                               \
  }

#define WtF_CHANGESTAT_D_CPP(...)                                         \
  _GET_OVERRIDE3(__VA_ARGS__, _WtF_CHANGESTAT_D_CPP_3, _WtF_CHANGESTAT_D_CPP_2, , )(__VA_ARGS__)

#define _WtW_CHANGESTAT_D_CPP_2(name, impl)                               \
  extern "C" SEXP w_ ## name (WtModelTerm *mtp, WtNetwork *nwp) {       \
    ergm::ErgmCppWtNetwork nw(nwp, true);                                     \
    ergm::ErgmCppWtModelTerm<> mt(mtp);                                 \
    impl;                                                               \
  }

#define _WtW_CHANGESTAT_D_CPP_3(name, StorageType, impl)                  \
  extern "C" SEXP w_ ## name (WtModelTerm *mtp, WtNetwork *nwp) {       \
    ergm::ErgmCppWtNetwork nw(nwp, true);                                     \
    ergm::ErgmCppWtModelTerm<StorageType> mt(mtp);                      \
    impl;                                                               \
  }

#define WtW_CHANGESTAT_D_CPP(...)                                         \
  _GET_OVERRIDE3(__VA_ARGS__, _WtW_CHANGESTAT_D_CPP_3, _WtW_CHANGESTAT_D_CPP_2, , )(__VA_ARGS__)

#define _WtX_CHANGESTAT_D_CPP_2(name, impl)                               \
  extern "C" void x_ ## name (unsigned int type, void *data, WtModelTerm *mtp, WtNetwork *nwp) { \
    ergm::ErgmCppWtNetwork nw(nwp, true);                                     \
    ergm::ErgmCppWtModelTerm<> mt(mtp);                                 \
    impl;                                                               \
  }

#define _WtX_CHANGESTAT_D_CPP_3(name, StorageType, impl)                  \
  extern "C" void x_ ## name (unsigned int type, void *data, WtModelTerm *mtp, WtNetwork *nwp) { \
    ergm::ErgmCppWtNetwork nw(nwp, true);                                     \
    ergm::ErgmCppWtModelTerm<StorageType> mt(mtp);                      \
    impl;                                                               \
  }

#define WtX_CHANGESTAT_D_CPP(...)                                         \
  _GET_OVERRIDE3(__VA_ARGS__, _WtX_CHANGESTAT_D_CPP_3, _WtX_CHANGESTAT_D_CPP_2, , )(__VA_ARGS__)
#define _WtZ_CHANGESTAT_D_CPP_2(name, impl)                               \
  extern "C" void z_ ## name (WtModelTerm *mtp, WtNetwork *nwp, Rboolean skip_s) { \
    ergm::ErgmCppWtNetwork nw(nwp, true);                                     \
    ergm::ErgmCppWtModelTerm<> mt(mtp);                                 \
    impl;                                                               \
  }

#define _WtZ_CHANGESTAT_D_CPP_3(name, StorageType, impl)                  \
  extern "C" void z_ ## name (WtModelTerm *mtp, WtNetwork *nwp, Rboolean skip_s) { \
    ergm::ErgmCppWtNetwork nw(nwp, true);                                     \
    ergm::ErgmCppWtModelTerm<StorageType> mt(mtp);                      \
    impl;                                                               \
  }

#define WtZ_CHANGESTAT_D_CPP(...)                                         \
  _GET_OVERRIDE3(__VA_ARGS__, _WtZ_CHANGESTAT_D_CPP_3, _WtZ_CHANGESTAT_D_CPP_2, , )(__VA_ARGS__)

// Undirected

#define _WtC_CHANGESTAT_U_CPP_2(name, impl)                               \
  extern "C" void c_ ## name (Vertex tail, Vertex head, double weight, WtModelTerm *mtp, WtNetwork *nwp, double edgestate) { \
    ergm::ErgmCppWtNetwork nw(nwp, false);                                     \
    ergm::ErgmCppWtModelTerm<> mt(mtp);                                 \
    impl;                                                               \
  }

#define _WtC_CHANGESTAT_U_CPP_3(name, StorageType, impl)                  \
  extern "C" void c_ ## name (Vertex tail, Vertex head, double weight, WtModelTerm *mtp, WtNetwork *nwp, double edgestate) { \
    ergm::ErgmCppWtNetwork nw(nwp, false);                                     \
    ergm::ErgmCppWtModelTerm<StorageType> mt(mtp);                      \
    impl;                                                               \
  }

#define WtC_CHANGESTAT_U_CPP(...)                                         \
  _GET_OVERRIDE3(__VA_ARGS__, _WtC_CHANGESTAT_U_CPP_3, _WtC_CHANGESTAT_U_CPP_2, , )(__VA_ARGS__)

#define _WtS_CHANGESTAT_U_CPP_2(name, impl)                               \
  extern "C" void s_ ## name (WtModelTerm *mtp, WtNetwork *nwp) {       \
    ergm::ErgmCppWtNetwork nw(nwp, false);                                     \
    ergm::ErgmCppWtModelTerm<> mt(mtp);                                 \
    impl;                                                               \
  }

#define _WtS_CHANGESTAT_U_CPP_3(name, StorageType, impl)                  \
  extern "C" void s_ ## name (WtModelTerm *mtp, WtNetwork *nwp) {       \
    ergm::ErgmCppWtNetwork nw(nwp, false);                                     \
    ergm::ErgmCppWtModelTerm<StorageType> mt(mtp);                      \
    impl;                                                               \
  }

#define WtS_CHANGESTAT_U_CPP(...)                                         \
  _GET_OVERRIDE3(__VA_ARGS__, _WtS_CHANGESTAT_U_CPP_3, _WtS_CHANGESTAT_U_CPP_2, , )(__VA_ARGS__)

#define _WtD_CHANGESTAT_U_CPP_2(name, impl)                               \
  extern "C" void d_ ## name (Edge ntoggles, Vertex *tails, Vertex *heads, double *weights, WtModelTerm *mtp, WtNetwork *nwp) { \
    ergm::ErgmCppWtNetwork nw(nwp, false);                                     \
    ergm::ErgmCppWtModelTerm<> mt(mtp);                                 \
    impl;                                                               \
  }

#define _WtD_CHANGESTAT_U_CPP_3(name, StorageType, impl)                  \
  extern "C" void d_ ## name (Edge ntoggles, Vertex *tails, Vertex *heads, double *weights, WtModelTerm *mtp, WtNetwork *nwp) { \
    ergm::ErgmCppWtNetwork nw(nwp, false);                                     \
    ergm::ErgmCppWtModelTerm<StorageType> mt(mtp);                      \
    impl;                                                               \
  }

#define WtD_CHANGESTAT_U_CPP(...)                                         \
  _GET_OVERRIDE3(__VA_ARGS__, _WtD_CHANGESTAT_U_CPP_3, _WtD_CHANGESTAT_U_CPP_2, , )(__VA_ARGS__)

#define _WtI_CHANGESTAT_U_CPP_2(name, impl)                               \
  extern "C" void i_ ## name (WtModelTerm *mtp, WtNetwork *nwp) {       \
    ergm::ErgmCppWtNetwork nw(nwp, false);                                     \
    ergm::ErgmCppWtModelTerm<> mt(mtp);                                 \
    impl;                                                               \
  }

#define _WtI_CHANGESTAT_U_CPP_3(name, StorageType, impl)                  \
  extern "C" void i_ ## name (WtModelTerm *mtp, WtNetwork *nwp) {       \
    ergm::ErgmCppWtNetwork nw(nwp, false);                                     \
    ergm::ErgmCppWtModelTerm<StorageType> mt(mtp);                      \
    impl;                                                               \
  }

#define WtI_CHANGESTAT_U_CPP(...)                                         \
  _GET_OVERRIDE3(__VA_ARGS__, _WtI_CHANGESTAT_U_CPP_3, _WtI_CHANGESTAT_U_CPP_2, , )(__VA_ARGS__)

#define _WtU_CHANGESTAT_U_CPP_2(name, impl)                               \
  extern "C" void u_ ## name (Vertex tail, Vertex head, double weight, WtModelTerm *mtp, WtNetwork *nwp, double edgestate) { \
    ergm::ErgmCppWtNetwork nw(nwp, false);                                     \
    ergm::ErgmCppWtModelTerm<> mt(mtp);                                 \
    impl;                                                               \
  }

#define _WtU_CHANGESTAT_U_CPP_3(name, StorageType, impl)                  \
  extern "C" void u_ ## name (Vertex tail, Vertex head, double weight, WtModelTerm *mtp, WtNetwork *nwp, double edgestate) { \
    ergm::ErgmCppWtNetwork nw(nwp, false);                                     \
    ergm::ErgmCppWtModelTerm<StorageType> mt(mtp);                      \
    impl;                                                               \
  }

#define WtU_CHANGESTAT_U_CPP(...)                                         \
  _GET_OVERRIDE3(__VA_ARGS__, _WtU_CHANGESTAT_U_CPP_3, _WtU_CHANGESTAT_U_CPP_2, , )(__VA_ARGS__)

#define _WtF_CHANGESTAT_U_CPP_2(name, impl)                               \
  extern "C" void f_ ## name (WtModelTerm *mtp, WtNetwork *nwp) {       \
    ergm::ErgmCppWtNetwork nw(nwp, false);                                     \
    ergm::ErgmCppWtModelTerm<> mt(mtp);                                 \
    impl;                                                               \
  }

#define _WtF_CHANGESTAT_U_CPP_3(name, StorageType, impl)                  \
  extern "C" void f_ ## name (WtModelTerm *mtp, WtNetwork *nwp) {       \
    ergm::ErgmCppWtNetwork nw(nwp, false);                                     \
    ergm::ErgmCppWtModelTerm<StorageType> mt(mtp);                      \
    impl;                                                               \
  }

#define WtF_CHANGESTAT_U_CPP(...)                                         \
  _GET_OVERRIDE3(__VA_ARGS__, _WtF_CHANGESTAT_U_CPP_3, _WtF_CHANGESTAT_U_CPP_2, , )(__VA_ARGS__)

#define _WtW_CHANGESTAT_U_CPP_2(name, impl)                               \
  extern "C" SEXP w_ ## name (WtModelTerm *mtp, WtNetwork *nwp) {       \
    ergm::ErgmCppWtNetwork nw(nwp, false);                                     \
    ergm::ErgmCppWtModelTerm<> mt(mtp);                                 \
    impl;                                                               \
  }

#define _WtW_CHANGESTAT_U_CPP_3(name, StorageType, impl)                  \
  extern "C" SEXP w_ ## name (WtModelTerm *mtp, WtNetwork *nwp) {       \
    ergm::ErgmCppWtNetwork nw(nwp, false);                                     \
    ergm::ErgmCppWtModelTerm<StorageType> mt(mtp);                      \
    impl;                                                               \
  }

#define WtW_CHANGESTAT_U_CPP(...)                                         \
  _GET_OVERRIDE3(__VA_ARGS__, _WtW_CHANGESTAT_U_CPP_3, _WtW_CHANGESTAT_U_CPP_2, , )(__VA_ARGS__)

#define _WtX_CHANGESTAT_U_CPP_2(name, impl)                               \
  extern "C" void x_ ## name (unsigned int type, void *data, WtModelTerm *mtp, WtNetwork *nwp) { \
    ergm::ErgmCppWtNetwork nw(nwp, false);                                     \
    ergm::ErgmCppWtModelTerm<> mt(mtp);                                 \
    impl;                                                               \
  }

#define _WtX_CHANGESTAT_U_CPP_3(name, StorageType, impl)                  \
  extern "C" void x_ ## name (unsigned int type, void *data, WtModelTerm *mtp, WtNetwork *nwp) { \
    ergm::ErgmCppWtNetwork nw(nwp, false);                                     \
    ergm::ErgmCppWtModelTerm<StorageType> mt(mtp);                      \
    impl;                                                               \
  }

#define WtX_CHANGESTAT_U_CPP(...)                                         \
  _GET_OVERRIDE3(__VA_ARGS__, _WtX_CHANGESTAT_U_CPP_3, _WtX_CHANGESTAT_U_CPP_2, , )(__VA_ARGS__)
#define _WtZ_CHANGESTAT_U_CPP_2(name, impl)                               \
  extern "C" void z_ ## name (WtModelTerm *mtp, WtNetwork *nwp, Rboolean skip_s) { \
    ergm::ErgmCppWtNetwork nw(nwp, false);                                     \
    ergm::ErgmCppWtModelTerm<> mt(mtp);                                 \
    impl;                                                               \
  }

#define _WtZ_CHANGESTAT_U_CPP_3(name, StorageType, impl)                  \
  extern "C" void z_ ## name (WtModelTerm *mtp, WtNetwork *nwp, Rboolean skip_s) { \
    ergm::ErgmCppWtNetwork nw(nwp, false);                                     \
    ergm::ErgmCppWtModelTerm<StorageType> mt(mtp);                      \
    impl;                                                               \
  }

#define WtZ_CHANGESTAT_U_CPP(...)                                         \
  _GET_OVERRIDE3(__VA_ARGS__, _WtZ_CHANGESTAT_U_CPP_3, _WtZ_CHANGESTAT_U_CPP_2, , )(__VA_ARGS__)
