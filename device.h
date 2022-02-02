// Declare a value lives on the GPU
#pragma once

namespace mandelbrot {

// Device<T> is a T that lives on the GPU
template<class T> struct Device {
private:
  T x;
};

}  // namespace mandelbrot
