# CLAUDE.md

## What this is

Computes an upper bound on the area of the Mandelbrot set via coefficients of the Böttcher series, using expansion arithmetic (unevaluated sums of doubles) since plain double runs out of precision around 2^23 terms. See README.md for the math.

## Building and testing

Dependencies (Mac): `brew install meson pkgconf flint openssl@3` (arb is part of flint ≥ 3; libpng and gmp also required). Linux details are in README.md.

    ./setup                        # configure build/{debug,debugoptimized,release} (wipes build/)
    meson compile -C build/debug   # build
    meson test -C build/debug      # run all tests
    meson test -C build/debug series           # run one test suite
    ./build/debug/series_test div1p exp        # run specific tests within a suite by name

Builds use `-Werror` with `warning_level=2`. CUDA (sm_80 and sm_90, via clang `-x cuda`) is optional; `*_cuda_test` targets only exist when CUDA is found. `meson.build` downloads `tinyformat.h` and `argparse.hpp` with curl at build time, so first builds need network access.

## Architecture

- **Two arithmetic worlds, tested against each other.** Fast computation uses `double` or expansion arithmetic (`expansion.h/cc`, `expansion_arith.h`); slow-but-rigorous reference computation uses the arb/flint interval library. Most tests generate random inputs, compute both ways, and compare — so test tolerances are tuned to specific flint RNG draws and can need loosening when flint changes its generator.
- **Flint C++ wrappers**: `arb_cc.h`, `acb_cc.h`, `arf_cc.h`, `fmpq_cc.h`, `mag_cc.h`, `poly.h`, `rand.h` are thin RAII wrappers around flint C types.
- **Series machinery**: `series.h/cc` implements power-series operations (mul, div, inv, exp, log and their `1p`-shifted variants) on top of `fft.h/cc`; `area.cc` drives the Böttcher-series area computation, and `arb_area.cc` is its rigorous arb counterpart.
- **Build-time code generation**: the `codelets` executable (from `codelets.cc`, `exp.cc`, `nearest.cc`, `sig.cc`) generates the `gen-*.h` headers (expansion arithmetic, FFT butterflies, series bases) that `area`/`expansion`/`series` compile against.
- **Cluster**: `cluster/` holds a one-H200 Kueue benchmark job (`pvc.yaml`, `bench.yaml`, `bench.sh`) that builds with a conda-forge toolchain.
- **CPU/GPU sharing**: the same `.cc` sources compile for CUDA when available; `loops.h`, `device.h`, and `cutil.h` abstract over serial, OpenMP, and CUDA execution.
- **Escape-time area**: `orbit.h` has resumable orbits (`Orbit<T>`, `OrbitDE`); `engine.h` runs them on CPU threads or a persistent-thread CUDA kernel; `tree.h/cc` is the certified quadtree Monte Carlo driven by `escape_tree`, deterministic from integer counts so CPU runs check GPU runs.
- **Julia set areas**: `julia.h/cc` computes μ(K(c)) for |c| < 1/4 to ~47 digits via the area transfer operator (collocation, iterative refinement with double GMRES); driver `julia_area`, notes in `notes/julia-area.md`.
