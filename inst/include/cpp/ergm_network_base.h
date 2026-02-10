/*  File inst/include/cpp/ergm_network_base.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#pragma once

#include <optional>
#include <tuple>
#include <type_traits>

#include "ergm_edgetree_types.h"

#include "combined_range.h"

// Template for both weighted and unweighted network wrappers
#include "combined_range.h"

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
  typename EdgetreeMinimumFunc,
  typename EdgetreeSuccessorFunc
>
class ErgmCppNetworkBase {
public:
  explicit ErgmCppNetworkBase(NetType* nwp)
    : dir(nwp->directed_flag != 0), n(nwp->nnodes), bip(nwp->bipartite), loops(nwp->loops_flag != 0), nwp_(nwp) {}

  class EdgeIterator {
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = EdgeValueType;
    using difference_type = std::ptrdiff_t;
    using pointer = value_type*;
    using reference = value_type&;

    EdgeIterator(TreeType* edges, Vertex node)
      : edges_(edges), e_(node && edges_[node].value ? EdgetreeMinimumFunc::call(edges_, node) : 0) {}
    EdgeIterator() : edges_(nullptr), e_(0) {}
    value_type operator*() const {
      if constexpr (std::is_same_v<value_type, Vertex>) {
        return edges_[e_].value;
      } else {
        return {edges_[e_].value, edges_[e_].weight};
      }
    }
    EdgeIterator& operator++() {
      e_ = EdgetreeSuccessorFunc::call(edges_, e_);
      return *this;
    }
    bool operator!=(const EdgeIterator& other) const {
      return (e_ == 0) != (other.e_ == 0);
    }
  private:
    TreeType* edges_;
    Edge e_;
  };

  class EdgeRange {
  public:
    using iterator = EdgeIterator;
    EdgeRange(TreeType* edges, Vertex node)
      : edges_(edges), node_(node) {}
    EdgeRange() : edges_(nullptr), node_(0) {}
    EdgeIterator begin() const { return EdgeIterator(edges_, node_); }
    EdgeIterator end() const { return EdgeIterator(edges_, 0); }
  private:
    TreeType* edges_;
    Vertex node_;
  };

  class FilteredEdgeIterator {
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = typename EdgeIterator::value_type;
    using difference_type = typename EdgeIterator::difference_type;
    using pointer = typename EdgeIterator::pointer;
    using reference = typename EdgeIterator::reference;

    FilteredEdgeIterator() : it_(), end_(), skip_(std::nullopt) {}
    FilteredEdgeIterator(EdgeIterator it, EdgeIterator end, std::optional<Vertex> skip)
      : it_(it), end_(end), skip_(skip) {
      advance();
    }
    value_type operator*() const { return *it_; }
    FilteredEdgeIterator& operator++() {
      ++it_;
      advance();
      return *this;
    }
    bool operator!=(const FilteredEdgeIterator& other) const {
      return it_ != other.it_;
    }
  private:
    EdgeIterator it_, end_;
    std::optional<Vertex> skip_;

    static Vertex edge_vertex(const value_type& value) {
      if constexpr (std::is_same_v<value_type, Vertex>) {
        return value;
      } else {
        return value.first;
      }
    }

    void advance() {
      if (!skip_.has_value()) {
        return;
      }
      while (it_ != end_ && edge_vertex(*it_) == *skip_) {
        ++it_;
      }
    }
  };

  class FilteredEdgeRange {
  public:
    using iterator = FilteredEdgeIterator;
    FilteredEdgeRange(EdgeRange range, std::optional<Vertex> skip = std::nullopt)
      : range_(range), skip_(skip) {}
    FilteredEdgeIterator begin() const { return FilteredEdgeIterator(range_.begin(), range_.end(), skip_); }
    FilteredEdgeIterator end() const { return FilteredEdgeIterator(range_.end(), range_.end(), skip_); }
  private:
    EdgeRange range_;
    std::optional<Vertex> skip_;
  };

  using CombinedEdgeRange = CombinedRange<FilteredEdgeRange>;

  CombinedEdgeRange out_neighbors(Vertex node) {
    if(dir) {
      return CombinedEdgeRange(FilteredEdgeRange(EdgeRange(nwp_->outedges, node)));
    } else {
      auto skip = has_loop(node) ? std::make_optional(node) : std::nullopt;
      return CombinedEdgeRange(
        FilteredEdgeRange(EdgeRange(nwp_->inedges, node)),
        std::make_optional(FilteredEdgeRange(EdgeRange(nwp_->outedges, node), skip))
      );
    }
  }
  CombinedEdgeRange in_neighbors(Vertex node) {
    if(dir) {
      return CombinedEdgeRange(FilteredEdgeRange(EdgeRange(nwp_->inedges, node)));
    } else {
      auto skip = has_loop(node) ? std::make_optional(node) : std::nullopt;
      return CombinedEdgeRange(
        FilteredEdgeRange(EdgeRange(nwp_->inedges, node)),
        std::make_optional(FilteredEdgeRange(EdgeRange(nwp_->outedges, node), skip))
      );
    }
  }
  CombinedEdgeRange neighbors(Vertex node) {
    auto skip = has_loop(node) ? std::make_optional(node) : std::nullopt;
    return CombinedEdgeRange(
      FilteredEdgeRange(EdgeRange(nwp_->inedges, node)),
      std::make_optional(FilteredEdgeRange(EdgeRange(nwp_->outedges, node), skip))
    );
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
    value_type operator*() const {
      if constexpr (std::is_same_v<value_type, std::tuple<Vertex, Vertex>>) {
        return std::make_tuple(tail_, *it_);
      } else {
        auto pair = *it_;
        return std::make_tuple(tail_, pair.first, pair.second);
      }
    }
    NetworkEdgeIterator& operator++() {
      ++it_;
      advance_vertex();
      return *this;
    }
    bool is_end() const {
      return tail_ > nnodes_;
    }
    bool operator!=(const NetworkEdgeIterator& other) const {
      // Consider both iterators at end if both are past-the-end
      return is_end() != other.is_end();
    }
  private:
    TreeType* outedges_;
    Vertex nnodes_;
    Vertex tail_;
    EdgeRange range_;
    EdgeIterator it_, end_it_;
    void advance_vertex() {
      while (tail_ <= nnodes_ && !(it_ != end_it_)) {
        ++tail_;
        if (tail_ <= nnodes_) {
          range_ = EdgeRange(outedges_, tail_);
          it_ = range_.begin();
          end_it_ = range_.end();
        }
      }
    }
  };

  class NetworkEdgeRange {
  public:
    NetworkEdgeRange(TreeType* outedges, Vertex nnodes)
      : outedges_(outedges), nnodes_(nnodes) {}
    NetworkEdgeIterator begin() const { return NetworkEdgeIterator(outedges_, nnodes_, 1); }
    NetworkEdgeIterator end() const { return NetworkEdgeIterator(outedges_, nnodes_, nnodes_ + 1); }
  private:
    TreeType* outedges_;
    Vertex nnodes_;
  };

  NetworkEdgeRange edges() const {
    return NetworkEdgeRange(nwp_->outedges, nwp_->nnodes);
  }

  class NodeIterator {
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = Vertex;
    using difference_type = std::ptrdiff_t;
    using pointer = Vertex*;
    using reference = Vertex&;
    NodeIterator(Vertex v) : v_(v) {}
    Vertex operator*() const { return v_; }
    NodeIterator& operator++() { ++v_; return *this; }
    bool operator!=(const NodeIterator& other) const { return v_ != other.v_; }
  private:
    Vertex v_;
  };

  class NodeRange {
  public:
    NodeRange(Vertex n) : start_(1), end_(n + 1) {}
    NodeRange(Vertex start, Vertex end) : start_(start), end_(end) {}
    NodeIterator begin() const { return NodeIterator(start_); }
    NodeIterator end() const { return NodeIterator(end_); }
    Vertex size() const { return end_ - start_; }
  private:
    Vertex start_;
    Vertex end_;
  };

  NodeRange nodes() const { return NodeRange(1, n + 1); }
  NodeRange b1() const { return NodeRange(1, bip + 1); }
  NodeRange b2() const { return NodeRange(bip + 1, n + 1); }

  ValueType operator()(Vertex tail, Vertex head) const {
    return GetEdgeFunc::call(tail, head, nwp_);
  }

  // Degree access methods
  Vertex out_degree(Vertex i) const {
    if (dir || bip) {
      return nwp_->outdegree[i];
    } else {
      return nwp_->outdegree[i] + nwp_->indegree[i] - loop_adjust(i);
    }
  }
  Vertex in_degree(Vertex i) const {
    if (dir || bip) {
      return nwp_->indegree[i];
    } else {
      return nwp_->indegree[i] + nwp_->outdegree[i] - loop_adjust(i);
    }
  }
  // Total degree convenience method (out + in, adjusted for loops)
  Vertex degree(Vertex i) const {
    return nwp_->outdegree[i] + nwp_->indegree[i] - loop_adjust(i);
  }
  const bool dir;
  const Vertex n;
  const Vertex bip;
  const bool loops;

private:
  NetType* nwp_;
  Vertex loop_adjust(Vertex i) const {
    return has_loop(i) ? 1 : 0;
  }
  bool has_loop(Vertex i) const {
    return loops && GetEdgeFunc::call(i, i, nwp_) != static_cast<ValueType>(0);
  }
};

} // namespace v1
} // namespace ergm
