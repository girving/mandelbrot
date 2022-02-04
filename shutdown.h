// Graceful shutdown on exit
#pragma once

#include <functional>
namespace mandelbrot {

using std::function;

// Register a function to be called from shutdown().
// The function may be called multiple times.
void on_shutdown(const function<void()>& f);

// Call all registered functions
void shutdown();

}  // namespace mandelbrot
