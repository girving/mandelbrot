// Shim if C++20 <bit> header is missing
#pragma once

#include <cstdint>
#include <type_traits>
namespace mandelbrot {
using std::is_same_v;
using std::is_unsigned_v;
}  // namespace mandelbrot

#if __has_include(<bit>)

#include <bit>
namespace mandelbrot {
using std::bit_ceil;
using std::countl_zero;
using std::countr_zero;
using std::has_single_bit;
}  // namespace mandelbrot

#else  // No <bit>,  so roll our own

namespace mandelbrot {

template<class I> I countl_zero(const I n) {
  static_assert(is_unsigned_v<I>);
  if constexpr (is_same_v<I,unsigned>) return n ? __builtin_clz(n) : 32;
  else if constexpr (is_same_v<I,unsigned long>) return n ? __builtin_clzl(n) : 64;
}

template<class I> I countr_zero(const I n) {
  static_assert(is_unsigned_v<I>);
  if constexpr (is_same_v<I,unsigned>) return n ? __builtin_ctz(n) : 32;
  else if constexpr (is_same_v<I,unsigned long>) return n ? __builtin_ctzl(n) : 64;
}

template<class I> I bit_ceil(const I n) {
  static_assert(is_unsigned_v<I>);
  return n <= 1 ? 1 : I(1) << (8*sizeof(I) - countl_zero(n - 1));
}

}  // namespace mandelbrot
#endif  // __has_include(<bit>)

// Bit twiddling functions not in <bit>
namespace mandelbrot {

template<class I> I byteswap(I n) {
  static_assert(is_unsigned_v<I>);
  if constexpr (is_same_v<I,uint32_t>) return __builtin_bswap32(n);
  else if constexpr (is_same_v<I,uint64_t>) return __builtin_bswap64(n);
}

template<class I> I bitreverse(I n) {
  static_assert(is_unsigned_v<I>);
#ifdef __clang__
  if constexpr (is_same_v<I,uint32_t>) return __builtin_bitreverse32(n);
  else if constexpr (is_same_v<I,uint64_t>) return __builtin_bitreverse64(n);
#else
  n = byteswap(n);
  const auto swap = [](I n, int s, unsigned mask) {
    const I m = mask * I(0x0101010101010101);
    return ((n & m) << s) | ((n & m<<s) >> s);
  };
  n = swap(n, 1, 0b01010101);
  n = swap(n, 2, 0b00110011);
  n = swap(n, 4, 0b00001111);
  return n;
#endif  // __clang__
}

}  // namespace mandelbrot
