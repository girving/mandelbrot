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

To install `bazel` and `arb`:

    # On Mac
    brew install bazel arb

    # On Linux
    # Install bazel via https://docs.bazel.build/versions/main/install-ubuntu.html
    sudo apt-get install libmpfr-dev libflint-dev libflint-arb-dev
