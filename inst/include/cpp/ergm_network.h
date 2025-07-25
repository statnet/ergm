#pragma once

#include <optional>
#include <tuple>

#include "combined_range.h"

extern "C" {
#include "ergm_edgetree.h"
}

class ErgmCppNetwork {
public:
  explicit ErgmCppNetwork(Network* nwp)
    : dir(nwp->directed_flag != 0), n(nwp->nnodes), bip(nwp->bipartite), nwp_(nwp) {}

  // Unified edge iterator
  class EdgeIterator {
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = Vertex;
    using difference_type = std::ptrdiff_t;
    using pointer = Vertex*;
    using reference = Vertex&;

    EdgeIterator(TreeNode* edges, Vertex node)
      : edges_(edges), e_(node && edges_[node].value? EdgetreeMinimum(edges_, node) : 0) {}
    EdgeIterator() : edges_(nullptr), e_(0) {}
    Vertex operator*() const { return edges_[e_].value; }
    EdgeIterator& operator++() {
      e_ = EdgetreeSuccessor(edges_, e_);
      return *this;
    }
    bool operator!=(const EdgeIterator& other) const { return (e_ == 0 || edges_[e_].value == 0) != (other.e_ == 0 || other.edges_[other.e_].value == 0); }

  private:
    TreeNode* edges_;
    Edge e_;
  };

  class EdgeRange {
  public:
    using iterator = EdgeIterator;
    EdgeRange(TreeNode* edges, Vertex node)
      : edges_(edges), node_(node) {}
    EdgeRange() : edges_(nullptr), node_(0) {}
    EdgeIterator begin() const { return EdgeIterator(edges_, node_); }
    EdgeIterator end() const { return EdgeIterator(edges_, 0); }
  private:
    TreeNode* edges_;
    Vertex node_;
  };

  // Use generic CombinedRange<EdgeRange> for combined edge iteration
  using CombinedEdgeRange = CombinedRange<EdgeRange>;

  CombinedEdgeRange out_neighbors(Vertex node) {
    if(dir) {
      return CombinedEdgeRange(EdgeRange(nwp_->outedges, node));
    } else {
      return CombinedEdgeRange(EdgeRange(nwp_->inedges, node), std::make_optional(EdgeRange(nwp_->outedges, node)));
    }
  }
  CombinedEdgeRange in_neighbors(Vertex node) {
    if(dir) {
      return CombinedEdgeRange(EdgeRange(nwp_->inedges, node));
    } else {
      return CombinedEdgeRange(EdgeRange(nwp_->inedges, node), std::make_optional(EdgeRange(nwp_->outedges, node)));
    }
  }
  CombinedEdgeRange neighbors(Vertex node) {
    return CombinedEdgeRange(EdgeRange(nwp_->inedges, node), std::make_optional(EdgeRange(nwp_->outedges, node)));
  }

  class NetworkEdgeIterator {
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = std::tuple<Vertex, Vertex>;
    using difference_type = std::ptrdiff_t;
    using pointer = value_type*;
    using reference = value_type&;

    NetworkEdgeIterator(TreeNode* outedges, Vertex nnodes, Vertex tail = 1)
      : outedges_(outedges), nnodes_(nnodes), tail_(tail), range_(outedges, tail <= nnodes ? tail : 0), it_(range_.begin()), end_it_(range_.end()) {
      advance_to_next_valid();
    }

    value_type operator*() const { return std::make_tuple(tail_, *it_); }
    NetworkEdgeIterator& operator++() {
      ++it_;
      advance_to_next_valid();
      return *this;
    }
    bool operator!=(const NetworkEdgeIterator& other) const {
      return tail_ != other.tail_ || it_ != other.it_;
    }

  private:
    TreeNode* outedges_;
    Vertex nnodes_;
    Vertex tail_;
    EdgeRange range_;
    EdgeIterator it_, end_it_;

    void advance_to_next_valid() {
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
    NetworkEdgeRange(TreeNode* outedges, Vertex nnodes)
      : outedges_(outedges), nnodes_(nnodes) {}
    NetworkEdgeIterator begin() const { return NetworkEdgeIterator(outedges_, nnodes_, 1); }
    NetworkEdgeIterator end() const { return NetworkEdgeIterator(outedges_, nnodes_, nnodes_ + 1); }
  private:
    TreeNode* outedges_;
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
  private:
    Vertex start_;
    Vertex end_;
  };

  NodeRange nodes() const { return NodeRange(1, n + 1); }
  NodeRange b1() const { return NodeRange(1, bip + 1); }
  NodeRange b2() const { return NodeRange(bip + 1, n + 1); }

  bool operator()(Vertex tail, Vertex head) const {
    return GetEdge(tail, head, nwp_);
  }

  const bool dir;
  const Vertex n;
  const Vertex bip;

private:
  Network* nwp_;
};
