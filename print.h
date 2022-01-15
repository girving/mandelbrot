// Print formatted content with a newline
#pragma once

#include "tinyformat.h"
namespace mandelbrot {

using std::cout;
using std::endl;
using tinyformat::format;

template<class... Args> void print(const char* fmt, const Args&... args) {
  cout << format(fmt, args...) << endl;
}

}  // namespace mandelbrot
