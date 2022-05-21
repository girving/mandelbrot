<p align="center">
  <img src="logo-2x.png?raw=true" width="282" height="256" title="Bottcher visualization">
</p>

# Mandelbrot set area via the Böttcher series

Let $\mathbb{C}$ be the complex plane, $M$ the [Mandelbrot set](https://en.wikipedia.org/wiki/Mandelbrot_set), and $D$ the closed unit disk.  There is an analytic [Böttcher map](https://en.wikipedia.org/wiki/External_ray)

$$\phi : \mathbb{C} - D \to \mathbb{C} - M$$

$$\phi(z) = z + \sum_n b_n z^{-n}$$

and the area of the Mandelbrot set is

$$\mu(M) = \pi \left(1 - \sum_n n b_n^2\right)$$

[Bittner et al. 2014](https://arxiv.org/abs/1410.1212) computed 5M terms of this series, resulting in the bound

$$\mu(M) \le 1.68288$$

We can compute out to $2^{27} = 134,217,728$ terms in a couple hours on an A100, producing

$$\mu(M) \le 1.651587035834859$$

We use [expansion arithmetic](https://people.eecs.berkeley.edu/~jrs/papers/robustr.pdf), representing numbers as
unevaluated sums of double precision numbers, as computing in double runs out of precision around
$2^{23}$ terms.

## Alas, the series approach isn't the fastest

The fastest of the mostly trustworthy methods I've seen for Mandelbrot area is [Fisher and Hill 1993](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.53.2337&rep=rep1&type=pdf), who use the [Koebe 1/4 theorem](https://en.wikipedia.org/wiki/Koebe_quarter_theorem) to prove (up to double precision) that quadree cells are either entirely outside or entirely inside the set.  Their bounds are

$$1.50296686 < \mu(M) < 1.57012937$$

Worse, accurate extrapolation of the series seems out of reach: out to $2^{26}$ terms the area contribution locally averages to a noisy power law fit with exponent -1.08, but as shown by [Bielefeld et al. 1988](https://archive.mpim-bonn.mpg.de/id/eprint/3259/1/preprint_1988_46.pdf) this would violate non-Hölder continuity at the boundary.  Attempts at extrapolating from our results out to infinity produce area estimates around 1.59 or 1.60, which is already contradicted by Fisher and Hill.

## Explanation, analysis, and downloads

[This colab](https://colab.research.google.com/drive/19FcWTtfXystwet4p06L2vMXOpP-r51ZH) has more explanation of the methods used, and some plots of the computed series.  For example, here is a plot of the locally averaged area contributions vs. the illusory -1.08 power law fit:

<p align="center">
  <img src="fit.svg?raw=true" title="Illusory power law fit">
</p>

For an `.npy` file with the series coefficients out to $2^k$ terms, set `$k` to 1, 2, ..., 26, or 27 and download

* `https://storage.googleapis.com/mandelbrot/numpy/f-k$k.npy`

The shape will be `[2^k, 2]` corresponding to expansion arithmetic with 2 doubles; sum across the last axis if you want a single double.

## Building

We use the [Meson](https://mesonbuild.com) build system, and the excellent high precision arithmetic library [Arb](https://arblib.org) for bootstrapping.  We also depend on clang even when compiling CUDA, to allow more recent C++ features.  To install dependencies:

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

For CUDA profiling, use

    # https://developer.nvidia.com/nsight-systems
    # https://docs.nvidia.com/nsight-systems/UserGuide/index.html#example-single-command-lines
    nsys profile --stats=true ./build/release/area_cuda_test
