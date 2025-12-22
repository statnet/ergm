/*  File inst/include/cpp/ergm_network.h in package ergm, part of the Statnet
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

#include "ergm_edgetree.h"

// Functor wrappers for static inline C functions for template use (unweighted)

namespace ergm {
inline namespace v1 {

struct GetEdgeFunc {
  static Rboolean call(Vertex tail, Vertex head, Network* nwp) {
    return GetEdge(tail, head, nwp);
  }
};
struct EdgetreeMinimumFunc {
  static Edge call(TreeNode* edges, Vertex node) {
    return EdgetreeMinimum(edges, node);
  }
};
struct EdgetreeSuccessorFunc {
  static Edge call(TreeNode* edges, Edge e) {
    return EdgetreeSuccessor(edges, e);
  }
};

using ErgmCppNetwork = ErgmCppNetworkBase<
  Network,
  TreeNode,
  Rboolean, // operator() return type
  Vertex, // EdgeIterator::value_type
  std::tuple<Vertex, Vertex>, // NetworkEdgeIterator::value_type
  GetEdgeFunc,
  EdgetreeMinimumFunc,
  EdgetreeSuccessorFunc
>;

} // namespace v1
} // namespace ergm
