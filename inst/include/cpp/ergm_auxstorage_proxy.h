/*  File inst/include/cpp/ergm_auxstorage_proxy.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#pragma once
#include <cstddef>

namespace ergm {
inline namespace v1 {

template<typename StructType>
class AuxStorageProxy {
public:
  StructType* ptr_;
  AuxStorageProxy(StructType* ptr) : ptr_(ptr) {}
  size_t size() const { return ptr_->n_aux; }
  void* operator[](size_t i) const { return ptr_->aux_storage[ptr_->aux_slots[i]]; }
  struct iterator {
    StructType* ptr_;
    size_t idx_;
    iterator(StructType* ptr, size_t idx) : ptr_(ptr), idx_(idx) {}
    void* operator*() const { return ptr_->aux_storage[ptr_->aux_slots[idx_]]; }
    iterator& operator++() { ++idx_; return *this; }
    bool operator!=(const iterator& other) const { return idx_ != other.idx_; }
  };
  iterator begin() const { return iterator(ptr_, 0); }
  iterator end() const { return iterator(ptr_, size()); }
};

} // namespace v1
} // namespace ergm
