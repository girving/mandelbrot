// A few overloads on top of tinyformat::format
#pragma once

#include "tinyformat.h"
namespace mandelbrot {

using std::string;
using tinyformat::format;

static inline string format() { return string(); }

}  // namespace mandelbrot
