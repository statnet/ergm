#pragma once

extern "C" {
#include "ergm_edgetree.h"
}

class ErgmCppNetwork {
public:
  explicit ErgmCppNetwork(Network* nwp)
    : nwp_(nwp), dir(nwp->directed_flag != 0), n(nwp->nnodes) {}

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
    // Return the neighbor vertex for this edge
    Vertex operator*() const { return edges_[e_].value; }
    EdgeIterator& operator++() {
      e_ = EdgetreeSuccessor(edges_, e_);
      return *this;
    }
    bool operator!=(const EdgeIterator& other) const { return e_ != other.e_; }

  private:
    TreeNode* edges_;
    Edge e_;
  };

  // Unified edge range
  class EdgeRange {
  public:
    EdgeRange(TreeNode* edges, Vertex node)
      : edges_(edges), node_(node) {}
    EdgeIterator begin() const { return EdgeIterator(edges_, node_); }
    EdgeIterator end() const { return EdgeIterator(edges_, 0); }
  private:
    TreeNode* edges_;
    Vertex node_;
  };

  class CombinedEdgeIterator {
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = Vertex;
    using difference_type = std::ptrdiff_t;
    using pointer = Vertex*;
    using reference = Vertex&;

    CombinedEdgeIterator(EdgeRange in_range, EdgeRange out_range, bool end = false)
      : out_range_(out_range), it_(in_range.begin()), in_end_(in_range.end()),
        out_end_(out_range_.end()), in_(!end) {
      if (in_ && !(it_ != in_end_)) {
        in_ = false;
        it_ = out_range_.begin();
      } else if (end) {
        it_ = out_range_.end();
      }
    }

    Vertex operator*() const {
      return *it_;
    }

    CombinedEdgeIterator& operator++() {
      ++it_;
      if (in_ && !(it_ != in_end_)) {
        in_ = false;
        it_ = out_range_.begin();
      }
      return *this;
    }

    bool operator!=(const CombinedEdgeIterator& other) const {
      return in_ != other.in_ || it_ != other.it_;
    }

  private:
    EdgeRange out_range_;
    EdgeIterator it_, in_end_, out_end_;
    bool in_;
  };

  class CombinedEdgeRange {
  public:
    CombinedEdgeRange(TreeNode* inedges, TreeNode* outedges, Vertex node)
      : in_range_(inedges, node), out_range_(outedges, node) {}
    CombinedEdgeIterator begin() const { return CombinedEdgeIterator(in_range_, out_range_, false); }
    CombinedEdgeIterator end() const { return CombinedEdgeIterator(in_range_, out_range_, true); }
  private:
    EdgeRange in_range_;
    EdgeRange out_range_;
  };

  EdgeRange out_neighbors(Vertex node) { return EdgeRange(nwp_->outedges, node); }
  EdgeRange in_neighbors(Vertex node) { return EdgeRange(nwp_->inedges, node); }
  CombinedEdgeRange neighbors(Vertex node) {
    return CombinedEdgeRange(nwp_->inedges, nwp_->outedges, node);
  }

  struct EdgePair {
    Vertex tail;
    Vertex head;
  };

  class NetworkEdgeIterator {
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = EdgePair;
    using difference_type = std::ptrdiff_t;
    using pointer = EdgePair*;
    using reference = EdgePair&;

    NetworkEdgeIterator(TreeNode* outedges, Vertex nnodes, Vertex tail = 1)
      : outedges_(outedges), nnodes_(nnodes), tail_(tail), range_(outedges, tail), it_(range_.begin()), end_it_(range_.end()) {
      advance_to_next_valid();
    }

    EdgePair operator*() const { return {tail_, *it_}; }
    NetworkEdgeIterator& operator++() {
      ++it_;
      if (!(it_ != end_it_)) {
        ++tail_;
        if (tail_ <= nnodes_) {
          range_ = EdgeRange(outedges_, tail_);
          it_ = range_.begin();
          end_it_ = range_.end();
          advance_to_next_valid();
        }
      }
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

  // Add more wrappers for other Network functions as needed...

  bool operator()(Vertex tail, Vertex head) const {
    return GetEdge(tail, head, nwp_);
  }

  const bool dir;
  const Vertex n;

private:
  Network* nwp_;
};
