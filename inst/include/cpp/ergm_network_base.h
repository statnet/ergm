/*  File inst/include/cpp/ergm_network_base.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2026 Statnet Commons
 */
#pragma once

#include <iterator>
#include <tuple>
#include <type_traits>

#include "ergm_debug.h"
#include "ergm_edgetree_types.h"
#include "inc/ergm_edgetree_common.do_not_include_directly.h"

namespace ergm {
inline namespace v1 {

// Template for both weighted and unweighted network wrappers

template <
  typename NetType,
  typename TreeType,
  typename ValueType, // type for operator() return value
  typename EdgeValueType, // type for EdgeIterator::value_type
  typename NetEdgeValueType, // type for NetworkEdgeIterator::value_type
  typename GetEdgeFunc,
  typename EdgetreeSearchFunc,
  typename EdgetreeMinimumFunc,
  typename EdgetreeSuccessorFunc
>
class ErgmCppNetworkBase {
public:
  explicit ErgmCppNetworkBase(NetType* nwp)
    : dir(nwp->directed_flag != 0), n(nwp->nnodes), bip(nwp->bipartite), ptr(nwp) {}

  class NeighborIterator {
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = EdgeValueType;
    using difference_type = std::ptrdiff_t;
    using pointer = value_type*;
    using reference = value_type&;

    ALWAYS_INLINE NeighborIterator(TreeType* first_edges, TreeType* second_edges, Vertex node, bool end = false)
      : edges_(end && second_edges ? second_edges : first_edges),
        edge_(0),
        second_edges_(second_edges),
        node_(node),
        second_phase_(end && second_edges) {
      if (!end) {
        edge_ = node_ && edges_[node_].value ? EdgetreeMinimumFunc::call(edges_, node_) : 0;
        advance_to_nonempty_second();
      }
    }
    ALWAYS_INLINE value_type operator*() const {
      if constexpr (std::is_same_v<value_type, Vertex>) {
        return edges_[edge_].value;
      } else {
        return {edges_[edge_].value, edges_[edge_].weight};
      }
    }
    ALWAYS_INLINE NeighborIterator& operator++() {
      edge_ = EdgetreeSuccessorFunc::call(edges_, edge_);
      advance_to_nonempty_second();
      return *this;
    }
     ALWAYS_INLINE bool operator!=(const NeighborIterator& other) const {
      return edge_ != other.edge_ || second_phase_ != other.second_phase_;
    }

  private:
    TreeType* edges_;
    Edge edge_;
    TreeType* second_edges_;
    Vertex node_;
    bool second_phase_;

     ALWAYS_INLINE void advance_to_nonempty_second() {
      if (!edge_ && !second_phase_ && second_edges_) {
        second_phase_ = true;
        edges_ = second_edges_;
        edge_ = node_ && edges_[node_].value ? EdgetreeMinimumFunc::call(edges_, node_) : 0;
      }
    }
  };

  class NeighborRange {
  public:
    using iterator = NeighborIterator;
     ALWAYS_INLINE NeighborRange(TreeType* first_edges, Vertex node, TreeType* second_edges = nullptr)
      : first_edges_(first_edges), second_edges_(second_edges), node_(node) {}
     ALWAYS_INLINE NeighborIterator begin() const { return NeighborIterator(first_edges_, second_edges_, node_); }
     ALWAYS_INLINE NeighborIterator end() const { return NeighborIterator(first_edges_, second_edges_, node_, true); }
  private:
    TreeType* first_edges_;
    TreeType* second_edges_;
    Vertex node_;
  };

  ALWAYS_INLINE NeighborRange out_neighbors(Vertex node) {
    if(dir) {
      return NeighborRange(ptr->outedges, node);
    } else {
      return NeighborRange(ptr->inedges, node, ptr->outedges);
    }
  }
  ALWAYS_INLINE NeighborRange in_neighbors(Vertex node) {
    if(dir) {
      return NeighborRange(ptr->inedges, node);
    } else {
      return NeighborRange(ptr->inedges, node, ptr->outedges);
    }
  }
  ALWAYS_INLINE NeighborRange neighbors(Vertex node) {
    return NeighborRange(ptr->inedges, node, ptr->outedges);
  }

  class NetworkEdgeIterator {
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = NetEdgeValueType;
    using difference_type = std::ptrdiff_t;
    using pointer = value_type*;
    using reference = value_type&;

  NetworkEdgeIterator(TreeType* outedges, Vertex nnodes, Vertex tail = 1)
    : outedges_(outedges), nnodes_(nnodes), tail_(tail), range_(outedges, tail <= nnodes ? tail : 0), it_(range_.begin()), end_it_(range_.end()) {
      advance_vertex();
    }
    ALWAYS_INLINE value_type operator*() const {
      if constexpr (std::is_same_v<value_type, std::tuple<Vertex, Vertex>>) {
        return std::make_tuple(tail_, *it_);
      } else {
        auto pair = *it_;
        return std::make_tuple(tail_, pair.first, pair.second);
      }
    }
    ALWAYS_INLINE NetworkEdgeIterator& operator++() {
      ++it_;
      advance_vertex();
      return *this;
    }
    ALWAYS_INLINE bool is_end() const {
      return tail_ > nnodes_;
    }
    ALWAYS_INLINE bool operator!=(const NetworkEdgeIterator& other) const {
      // Consider both iterators at end if both are past-the-end
      return is_end() != other.is_end();
    }
  private:
    TreeType* outedges_;
    Vertex nnodes_;
    Vertex tail_;
    NeighborRange range_;
    NeighborIterator it_, end_it_;
    ALWAYS_INLINE void advance_vertex() {
      while (tail_ <= nnodes_ && !(it_ != end_it_)) {
        ++tail_;
        if (tail_ <= nnodes_) {
          range_ = NeighborRange(outedges_, tail_);
          it_ = range_.begin();
          end_it_ = range_.end();
        }
      }
    }
  };

  class NetworkEdgeRange {
  public:
    ALWAYS_INLINE NetworkEdgeRange(TreeType* outedges, Vertex nnodes)
      : outedges_(outedges), nnodes_(nnodes) {}
    ALWAYS_INLINE NetworkEdgeIterator begin() const { return NetworkEdgeIterator(outedges_, nnodes_, 1); }
    ALWAYS_INLINE NetworkEdgeIterator end() const { return NetworkEdgeIterator(outedges_, nnodes_, nnodes_ + 1); }
  private:
    TreeType* outedges_;
    Vertex nnodes_;
  };

  NetworkEdgeRange edges() const {
    return NetworkEdgeRange(ptr->outedges, ptr->nnodes);
  }

  class NodeIterator {
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = Vertex;
    using difference_type = std::ptrdiff_t;
    using pointer = Vertex*;
    using reference = Vertex&;
    NodeIterator(Vertex v) : v_(v) {}
    ALWAYS_INLINE Vertex operator*() const { return v_; }
    ALWAYS_INLINE NodeIterator& operator++() { ++v_; return *this; }
    ALWAYS_INLINE bool operator!=(const NodeIterator& other) const { return v_ != other.v_; }
  private:
    Vertex v_;
  };

  class NodeRange {
  public:
    ALWAYS_INLINE NodeRange(Vertex n) : start_(1), end_(n + 1) {}
    ALWAYS_INLINE NodeRange(Vertex start, Vertex end) : start_(start), end_(end) {}
    ALWAYS_INLINE NodeIterator begin() const { return NodeIterator(start_); }
    ALWAYS_INLINE NodeIterator end() const { return NodeIterator(end_); }
    ALWAYS_INLINE Vertex size() const { return end_ - start_; }
  private:
    Vertex start_;
    Vertex end_;
  };

  NodeRange nodes() const { return NodeRange(1, n + 1); }
  NodeRange b1() const { return NodeRange(1, bip + 1); }
  NodeRange b2() const { return NodeRange(bip + 1, n + 1); }

  ALWAYS_INLINE ValueType operator()(Vertex tail, Vertex head) const {
    return GetEdgeFunc::call(tail, head, ptr);
  }

  ALWAYS_INLINE ValueType operator()(Vertex tail, Vertex head, bool dir) const {
    if (dir) return EdgetreeSearchFunc::call(tail, head, ptr->outedges) ? TRUE : FALSE;
    else return EdgetreeSearchFunc::call(MIN(tail, head), MAX(tail, head), ptr->outedges) ? TRUE : FALSE;
  }

  // Degree access methods
  ALWAYS_INLINE Vertex out_degree(Vertex i) const {
    if (dir || bip) {
      return ptr->outdegree[i];
    } else {
      return ptr->outdegree[i] + ptr->indegree[i];
    }
  }
  ALWAYS_INLINE Vertex in_degree(Vertex i) const {
    if (dir || bip) {
      return ptr->indegree[i];
    } else {
      return ptr->indegree[i] + ptr->outdegree[i];
    }
  }
  // Total degree convenience method (always out + in)
  ALWAYS_INLINE Vertex degree(Vertex i) const {
    return ptr->outdegree[i] + ptr->indegree[i];
  }
  const bool dir;
  const Vertex n;
  const Vertex bip;
  NetType* ptr;
};

} // namespace v1
} // namespace ergm
