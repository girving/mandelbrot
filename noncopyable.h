// Noncopyable
#pragma once

namespace mandelbrot {

struct Noncopyable {
  Noncopyable() = default;
  Noncopyable(const Noncopyable&) = delete;
  void operator=(const Noncopyable&) = delete;
};

}  // namespace mandelbrot
