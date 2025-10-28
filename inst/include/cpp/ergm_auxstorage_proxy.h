#pragma once
#include <cstddef>

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
