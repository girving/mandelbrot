# Thoughts on a Rust rewrite

Status: shelved (September 2026).  This records the analysis in case we revisit it.

## What would get better

- **Floating-point semantics.** Expansion arithmetic depends on exact rounding of each operation.  Rust never
  contracts `a * b + c` into an FMA implicitly; FMA is only via explicit `mul_add`.  In C++ we rely on
  compiler flags (and hope) to keep contraction off, and a stray `-ffp-contract=fast` silently breaks
  error-free transformations.
- **Build and tooling.** cargo replaces meson plus the curl-at-build-time downloads of `tinyformat.h` and
  `argparse.hpp` (`clap`, `format!` are built in or standard crates).  `cargo test` with per-test filtering
  and `#[test]` beside the code replaces our hand-rolled `tests.h` macros.
- **Parallelism.** `rayon` gives work-stealing `par_iter` for the CPU paths (`escape_tree`, `hyperbolic`)
  with data-race freedom checked at compile time, replacing the ad hoc thread pools and mutexes.
- **Generics.** `Expansion<const N: usize>` works on stable const generics, and traits make the "same
  algorithm over double / Expansion / arb" pattern cleaner than templates plus ad hoc overloads.
- **Safety** around the RAII wrappers of flint types (`arb_cc.h` etc.): `Drop` plus ownership rules remove
  a class of aliasing and double-free mistakes.

## What would get worse

- **Shared CPU/GPU source.** Today the same `.cc` files compile for both host and CUDA (`loops.h`,
  `device.h`, `cutil.h`).  No Rust option gives that seamlessly on stable today; see below.
- **Const-generic arithmetic.** Expansion products want `Expansion<{N + M}>`, which needs
  `generic_const_exprs` (still unstable).  Workarounds: fixed sizes, macros, or typenum.
- **flint/arb bindings.** Only raw `-sys` style bindings (e.g. `flint3-sys`, `arb-sys`) exist; we would write
  and maintain our own safe wrappers — roughly what `arb_cc.h`, `acb_cc.h`, `arf_cc.h`, `fmpq_cc.h`,
  `mag_cc.h`, `poly.h` already are.
- **Code generation.** The `codelets` generator would become a `build.rs` or proc macro — fine, but a port.
- **SIMD.** `std::simd` is still nightly; stable needs `std::arch` intrinsics or crates.
- **Rewrite cost and risk** of a working, test-heavy numerical codebase, with the flint-RNG-tuned test
  tolerances needing re-derivation.

## CUDA in Rust (state as of September 2026)

- **cudarc**: mature, stable host-side driver/NVRTC/cuBLAS bindings.  Kernels still written in CUDA C++ or
  PTX.  The lowest-risk path.
- **NVIDIA's CUDA Rust** (announced September 2026,
  https://developer.nvidia.com/blog/introducing-cuda-rust-two-tracks-for-writing-gpu-kernels/):
  *cuda-oxide* (SIMT kernels in Rust, nightly) and *cutile-rs* (tile-based kernels, stable Rust 1.89+).
  Promising but very new.
- **Rust-CUDA** (`rustc_codegen_nvvm`): revived community project, nightly-pinned.
- **CubeCL**: portable GPU kernels via a Rust DSL (CUDA/ROCm/WGPU); good for tensor-shaped work, less
  natural for our scalar-heavy iteration loops.

"Idiomatic Rust CUDA" is plausible now but not yet as boring as CUDA C++.

## Recommendation

Don't rewrite the whole thing.  If we try Rust, prototype only the escape-time Monte Carlo (`escape.cc`,
`escape_tree.cc`): it is self-contained, plain `f64` (no expansion arithmetic or flint), embarrassingly
parallel, and the part we would run on GPUs.  A good trial: rayon on CPU plus one GPU path (cudarc with a
CUDA C++ kernel, or cutile-rs/cuda-oxide), compared against the C++ version for throughput and identical
classifications on a fixed sample set.
