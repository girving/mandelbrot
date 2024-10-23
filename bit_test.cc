// Bit tests

#include "bit.h"
#include "print.h"
#include "tests.h"
#include <random>
namespace mandelbrot {
namespace {

using std::mt19937;
using std::mt19937_64;
using std::remove_reference_t;

TEST(countl_zero) {
  const auto test = [](const auto zero) {
    typedef decltype(zero) I;
    const int bits = 8 * sizeof(zero);
    ASSERT_EQ(countl_zero(I(0)), bits);
    for (int i = 0; i < bits; i++)
      for (int j = 0; j <= i; j++)
        ASSERT_EQ(countl_zero(I(1) << i | I(1) << j), bits - 1 - i);
  };
  test(uint32_t(0));
  test(uint64_t(0));
}

TEST(countr_zero) {
  const auto test = [](const auto zero) {
    typedef decltype(zero) I;
    const int bits = 8 * sizeof(zero);
    ASSERT_EQ(countr_zero(I(0)), bits);
    for (int i = 0; i < bits; i++)
      for (int j = 0; j <= i; j++)
        ASSERT_EQ(countr_zero(I(1) << i | I(1) << j), j);
  };
  test(uint32_t(0));
  test(uint64_t(0));
}

TEST(bit_ceil) {
  const auto test = [](const auto zero) {
    typedef decltype(zero) I;
    const int bits = 8 * sizeof(zero);
    ASSERT_EQ(bit_ceil(I(0)), I(1));
    for (int i = 0; i < bits; i++) {
      ASSERT_EQ(bit_ceil(I(1) << i), I(1) << i)
          << tfm::format("i %d, 1<<i %d, bit_ceil %d", i, I(1) << i, bit_ceil(I(1) << i));
      if (i + 1 < bits) {
        for (int j = 0; j < i; j++)
          ASSERT_EQ(bit_ceil(I(1) << i | I(1) << j), I(1) << (i+1));
      }
    }
  };
  test(uint32_t(0));
  test(uint64_t(0));
}

TEST(byteswap) {
  mt19937_64 mt(7);
  const auto test = [&mt](auto zero) {
    typedef decltype(zero) I;
    const auto slow = [](I n) {
      const int bytes = sizeof(I);
      I s = 0;
      for (int i = 0; i < bytes; i++)
        s |= ((n >> 8*i) & 0xff) << 8*(bytes-1-i);
      return s;
    };
    for (int i = 0; i < 1024; i++) {
      const I n = mt();
      ASSERT_EQ(byteswap(n), slow(n));
    }
  };
  test(uint32_t(0));
  test(uint64_t(0));
}

TEST(bitreverse) {
  mt19937_64 mt(7);
  const auto test = [&mt](auto zero) {
    typedef decltype(zero) I;
    const auto slow = [](I n) {
      const int bits = 8*sizeof(I);
      I s = 0;
      for (int i = 0; i < bits; i++)
        s |= ((n >> i) & 1) << (bits-1-i);
      return s;
    };
    for (int i = 0; i < 1024; i++) {
      const I n = mt();
      ASSERT_EQ(bitreverse(n), slow(n));
    }
  };
  test(uint32_t(0));
  test(uint64_t(0));
}

}  // namespace
}  // namespace mandelbrot
