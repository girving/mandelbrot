# Mandelbrot set area via the Böttcher series

Let C be the complex plane, M the [Mandelbrot set](https://en.wikipedia.org/wiki/Mandelbrot_set), and D the closed unit disk.  There is an analytic [Böttcher map](https://en.wikipedia.org/wiki/External_ray)

    ɸ : C - D → C - M
    ɸ(z) = z + sum_n b_n z^(-n)

and the area of the Mandelbrot set is

    𝜇(M) = sum_n n b_n^2

[Bittner et al. 2014](https://arxiv.org/abs/1410.1212) computed 5M terms of this series, resulting in the bound

    𝜇(M) ≤ 1.68288

Here we try to compute more terms.

## Dependencies

To install `meson` and `arb`:

    # On Mac
    brew install meson arb

    # On Linux
    sudo apt-get install clang-14 libc++-14-dev libc++abi-14-dev  # Need libc++-14 for std::format
    sudo apt-get install python3 python3-pip python3-setuptools python3-wheel ninja-build  # Meson deps
    sudo apt-get install libmpfr-dev libflint-dev libflint-arb-dev
    pip3 install --user meson
