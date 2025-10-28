#pragma once

#include "ergm_network_base.h"

extern "C" {
#include "ergm_edgetree.h"
}

// Functor wrappers for static inline C functions for template use (unweighted)

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
