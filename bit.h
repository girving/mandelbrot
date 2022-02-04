// A few bit twiddling functions beyond <bit>
#pragma once

#include <bit>
#include <cstdint>
#include <type_traits>
namespace mandelbrot {

using std::is_same_v;
using std::is_unsigned_v;
using std::bit_ceil;
using std::countl_zero;
using std::countr_zero;
using std::has_single_bit;

template<class I> I byteswap(I n) {
  static_assert(is_unsigned_v<I>);
  if constexpr (is_same_v<I,uint32_t>) return __builtin_bswap32(n);
  else if constexpr (is_same_v<I,uint64_t>) return __builtin_bswap64(n);
}

template<class I> I bitreverse(I n) {
  static_assert(is_unsigned_v<I>);
  if constexpr (is_same_v<I,uint32_t>) return __builtin_bitreverse32(n);
  else if constexpr (is_same_v<I,uint64_t>) return __builtin_bitreverse64(n);
}

}  // namespace mandelbrot
