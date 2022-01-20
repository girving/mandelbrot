// Assertions
#pragma once

#include "format.h"
namespace mandelbrot {

// Print a message and abort
void die(const string& msg) __attribute__((noreturn, cold));
template<class... Args> static inline void __attribute__((noreturn, cold)) die(const Args&... args) {
  die(format(args...));
}

// Assert that still happens in optimized mode
#define slow_assert(condition, ...) \
  ((condition) ? (void)0 : mandelbrot::assertion_failed( \
      __PRETTY_FUNCTION__, __FILE__, __LINE__, #condition, format(__VA_ARGS__)))

// If slow_assert fails...
void __attribute__((noreturn, cold))
assertion_failed(const char* function, const char* file, unsigned int line, const char* condition,
                 const string& message);

}  // namespace mandelbrot
