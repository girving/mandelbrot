// Print formatted content with a newline
#pragma once

#include "format.h"
namespace mandelbrot {

static inline void print() { std::cout << std::endl; }
template<class T> static inline void print(T&& x) { std::cout << x << std::endl; }

template<class... Args> void print(const char* fmt, const Args&... args) {
  std::cout << format(fmt, args...) << std::endl;
}

}  // namespace mandelbrot
