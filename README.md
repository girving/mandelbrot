<p align="center">
  <img src="https://github.com/girving/mandelbrot/blob/main/logo.png?raw=true" title="Bottcher visualization">
</p>

# Mandelbrot set area via the Böttcher series

Let C be the complex plane, M the [Mandelbrot set](https://en.wikipedia.org/wiki/Mandelbrot_set), and D the closed unit disk.  There is an analytic [Böttcher map](https://en.wikipedia.org/wiki/External_ray)

<p align="center">
  <img src="https://render.githubusercontent.com/render/math?math={\phi : \mathbb{C} - D \to \mathbb{C} - M}#gh-light-mode-only">
  <img src="https://render.githubusercontent.com/render/math?math={\color{white}\phi : \mathbb{C} - D \to \mathbb{C} - M}#gh-dark-mode-only"><br/>
  <img src="https://render.githubusercontent.com/render/math?math={\phi(z) = z %2B \sum_n b_n z^{-n}}#gh-light-mode-only">
  <img src="https://render.githubusercontent.com/render/math?math={\color{white}\phi(z) = z %2B \sum_n b_n z^{-n}}#gh-dark-mode-only">
</p>

and the area of the Mandelbrot set is

<p align="center">
  <img src="https://render.githubusercontent.com/render/math?math={\mu(M) = \sum_n n b_n^2}#gh-light-mode-only">
  <img src="https://render.githubusercontent.com/render/math?math={\color{white}\mu(M) = \sum_n n b_n^2}#gh-dark-mode-only">
</p>

[Bittner et al. 2014](https://arxiv.org/abs/1410.1212) computed 5M terms of this series, resulting in the bound

<p align="center">
  <img src="https://render.githubusercontent.com/render/math?math={\mu(M) \le 1.68288}#gh-light-mode-only">
  <img src="https://render.githubusercontent.com/render/math?math={\color{white}\mu(M) \le 1.68288}#gh-dark-mode-only">
</p>

We can compute out to 2<sup>27</sup> = 134,217,728 terms in a couple hours on an A100, producing

<p align="center">
  <img src="https://render.githubusercontent.com/render/math?math={\mu(M) \le 1.651587035834859}#gh-light-mode-only">
  <img src="https://render.githubusercontent.com/render/math?math={\color{white}\mu(M) \le 1.651587035834859}#gh-dark-mode-only">
</p>

We use [expansion arithmetic](https://people.eecs.berkeley.edu/~jrs/papers/robustr.pdf), representing numbers as
unevaluated sums of double precision numbers, as computing in double runs out of precision around
2<sup>23</sup> terms.

## Building

First, install dependencies:

    # On Mac
    brew install meson arb

    # On Debian Buster
    echo deb http://apt.llvm.org/buster/ llvm-toolchain-buster-13 main | sudo tee -a /etc/apt/sources.list
    echo deb-src http://apt.llvm.org/buster/ llvm-toolchain-buster-13 main | sudo tee -a /etc/apt/sources.list
    sudo apt-get install clang-13 libc++-13-dev libc++abi-13-dev libomp-13-dev \
        python3 python3-pip python3-setuptools python3-wheel ninja-build \
        libmpfr-dev libflint-dev libflint-arb-dev libssl-dev
    pip3 install --user meson

Then build and test with

    ./setup
    cd build/release  # Or build/debug
    meson compile
    meson test

## CUDA profiling

    # https://developer.nvidia.com/nsight-systems
    # https://docs.nvidia.com/nsight-systems/UserGuide/index.html#example-single-command-lines
    nsys profile --stats=true ./build/release/area_cuda_test
