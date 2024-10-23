// Stringify and join a list
#pragma once

#include "format.h"
namespace mandelbrot {

template<class C> string join(const C& ss, const string& sep = ", ") {
  string j;
  for (const auto& s : ss) {
    if (j.size()) j += sep;
    j += tfm::format("%s", s);
  }
  return j;
}

}  // namespace mandelbrot
