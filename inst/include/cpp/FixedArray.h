/*  File inst/include/cpp/FixedArray.h in package ergm, part of the Statnet
 *  suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2026 Statnet Commons
 */
#pragma once

#include <cstddef>
#include <stdexcept>
#include <iterator>

namespace ergm {
inline namespace v1 {

template<typename T>
class FixedArray {
public:
  using value_type = T;
  using size_type = std::size_t;
  using reference = T&;
  using const_reference = const T&;
  using pointer = T*;
  using const_pointer = const T*;
  using iterator = T*;
  using const_iterator = const T*;

  // Wrap an existing C array
  FixedArray(pointer arr, size_type n)
    : size_(n), data_(arr) {}

  // Copy constructor: only copy the pointer and size
  FixedArray(const FixedArray& other)
    : size_(other.size_), data_(other.data_) {}

  // Assignment operator: only copy the pointer and size
  FixedArray& operator=(const FixedArray& other) {
    if(this != &other) {
      size_ = other.size_;
      data_ = other.data_;
    }
    return *this;
  }

  // Destructor: do not delete the array
  ~FixedArray() = default;

  reference operator[](size_type idx) {
    return data_[idx];
  }

  const_reference operator[](size_type idx) const {
    return data_[idx];
  }

  reference at(size_type idx) {
    if(idx >= size_) throw std::out_of_range("FixedArray::at");
    return data_[idx];
  }

  const_reference at(size_type idx) const {
    if(idx >= size_) throw std::out_of_range("FixedArray::at");
    return data_[idx];
  }

  size_type size() const noexcept {
    return size_;
  }

  iterator begin() noexcept { return data_; }
  iterator end() noexcept { return data_ + size_; }
  const_iterator begin() const noexcept { return data_; }
  const_iterator end() const noexcept { return data_ + size_; }
  const_iterator cbegin() const noexcept { return data_; }
  const_iterator cend() const noexcept { return data_ + size_; }

  pointer data() noexcept { return data_; }
  const_pointer data() const noexcept { return data_; }

private:
  size_type size_;
  pointer data_;
};

} // namespace v1
} // namespace ergm
