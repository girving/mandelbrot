// Preprocessor control for codelets
#pragma once

namespace mandelbrot {

#if CODELETS
static constexpr bool codelets = true;
#else
static constexpr bool codelets = false;
#endif

}  // namespace mandelbrot
