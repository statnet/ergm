/*  File inst/include/cpp/ergm_wtnetwork.h in package ergm, part of the Statnet
 *  suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#pragma once

#include "ergm_network_base.h"

#include "ergm_wtedgetree.h"

// Functor wrappers for static inline C functions for template use (weighted)

namespace ergm {
inline namespace v1 {

struct WtGetEdgeFunc {
  static double call(Vertex tail, Vertex head, WtNetwork* nwp) {
    return WtGetEdge(tail, head, nwp);
  }
};
struct WtEdgetreeMinimumFunc {
  static Edge call(WtTreeNode* edges, Vertex node) {
    return WtEdgetreeMinimum(edges, node);
  }
};
struct WtEdgetreeSuccessorFunc {
  static Edge call(WtTreeNode* edges, Edge e) {
    return WtEdgetreeSuccessor(edges, e);
  }
};

using ErgmCppWtNetwork = ErgmCppNetworkBase<
  WtNetwork,
  WtTreeNode,
  double, // operator() return type
  std::pair<Vertex, double>, // EdgeIterator::value_type
  std::tuple<Vertex, Vertex, double>, // NetworkEdgeIterator::value_type
  WtGetEdgeFunc,
  WtEdgetreeMinimumFunc,
  WtEdgetreeSuccessorFunc
>;

} // namespace v1
} // namespace ergm
