// Arrays that live on host or device
#pragma once

#include "cutil.h"
namespace mandelbrot {

using std::conditional_t;
using std::enable_if_t;
using std::is_const_v;
using std::is_trivially_copyable_v;
using std::is_trivially_destructible_v;

template<class T> struct Array : public Noncopyable {
  typedef remove_const_t<Undevice<T>> S;
  static_assert(is_trivially_copyable_v<S>);
  static_assert(is_trivially_destructible_v<S>);
  typedef T* Data;
protected:
  struct Unusable {};
  struct Deleter { void operator()(T* p) const {
    if (!p) return;
    if constexpr (is_device<T>) cuda_check(cudaFreeAsync(undevice(p), stream()));
    else free(p);
  }};
  int64_t size_;  // Allocated terms
  unique_ptr<T,Deleter> x;
public:

  Array() : size_(0) {}
  Array(int64_t size) : size_(relu(size)) {
    if (size_) {
      T* p;
      if constexpr (is_device<T>) cuda_check(cudaMallocAsync(&p, size_ * sizeof(S), stream()));
      else p = static_cast<T*>(malloc(size_ * sizeof(T)));
      x.reset(p);
    }
  }
  Array(Array&& y) : size_(y.size_), x(move(y.x)) { y.size_ = 0; }

  int64_t size() const { return size_; }
  T* data() const { return x.get(); }
  T& operator[](const int64_t i) const { return x.get()[i]; }

  std::span<T> span() const { return std::span<T>(x.get(), size_t(size_)); }
  operator std::span<T>() const { return span(); }
  operator std::span<conditional_t<is_const_v<T>,Unusable,const T>>() const { return span(); }

  void swap(Array& y) { x.swap(y.x); std::swap(size_, y.size_); }
  void clear() { x.reset(); size_ = 0; }
};

}  // namespace mandelbrot
