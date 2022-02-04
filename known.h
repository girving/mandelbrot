// Known Mandelbrot area bounds
#pragma once

namespace mandelbrot {

struct Known {
  int k;  // Estimate from 2^k terms
  double time;  // time of refinement step in seconds (negative for unknown)
  int prec;  // Arb precision
  const char* value;  // serialized arb_t value
};
extern const Known known_areas[16+1];

// 8 terms of f and g = log f:
//   f = [(1 +/- 0), (-0.5 +/- 0), (0.125 +/- 0), (-0.25 +/- 9.3345e-61),
//        (0.117 +/- 1.3224e-60), (2.33e-62 +/- 1.7859e-60), (-0.0459 +/- 3.3449e-60), (-0.0625 +/- 4.9784e-60)]
//   g = [(0 +/- 0), (-0.5 +/- 0), (0 +/- 0), (-0.229 +/- 1.1927e-60),
//        (7.78e-62 +/- 1.7891e-60), (0.0289 +/- 2.6331e-60), (-0.0625 +/- 4.8099e-60), (-0.0836 +/- 7.5565e-60)]

// Arb history:
//   15jan2022, prec 200:
//     k 11 refine, 2.06 s: mu = 1.854656777 +/- 6.3116e-18
//   16jan2022, prec 200, Newton refinement:
//     k 11 refine, 5.65 s: mu = 1.854656777 +/- 1.2799e-37
//     k 12 refine, 13.9 s: mu = 1.834655733 +/- 9.6986e-33
//   17jan2022, prec 200, log1p_exp_shift:
//     k 11 refine, 20.0 s: mu = 1.854656777 +/- 2.3030e-45
//     k 12 refine, 61.5 s: mu = 1.834655733 +/- 7.6881e-43
//     k 13 refine, 169 s:  mu = 1.806178886 +/- 4.6439e-40
//     k 14 refine, 422 s:  mu = 1.786389717 +/- 3.9121e-37
//   17jan2022, prec 200, poly_inv_refine → f(y0)/f'(y):
//     k 11 refine, 25.1 s: mu = 1.854656777 +/- 1.1026e-45
//     k 12 refine, 72.4 s: mu = 1.834655733 +/- 3.3779e-43
//     k 13 refine, 206 s:  mu = 1.806178886 +/- 1.8539e-40
//     k 14 refine, 694 s:  mu = 1.786389717 +/- 1.4147e-37
//   17jan2022, prec 200, solve in g = log f space:
//     k 11 refine, 26.7 s: mu = 1.854656777 +/- 8.9856e-46
//     k 12 refine, 84.4 s: mu = 1.834655733 +/- 2.7409e-43

// Bittner et al.'s Tables 1 and 2 (https://arxiv.org/abs/1410.1212):
//     0.5M: 1.72 (Ewing-Schober)
//     1.0M: 1.70393
//     1.5M: 1.69702
//     2.0M: 1.69388
//     2.5M: 1.69096
//     3.0M: 1.68895, 9 days
//     3.5M: 1.6874,  10.8 days
//     4.0M: 1.68633, 12.5 days
//     4.5M: 1.68447, 14.4 days
//     5.0M: 1.68288, 16.2 days

// Series history:
//   22jan2022, double, fft_mul:
//     k 14, 39.5 s: mu = 1.786389717, error = 7.11e-15
//   22jan2022, double, no bit reverse:
//     k 14, 36.1 s: mu = 1.786389717, error = 2.11e-14
//   27jan2022, double, srfft:
//     k 14, 17.9 s: mu = 1.786389717, error = 1.82e-14
//   27jan2022, double, cache twiddle factors:
//     k 14, 2.90 s: mu = 1.786389717, error = 1.82e-14
//     k 15, 6.99 s: mu = 1.766837674, error = 3.55e-14
//     k 16, 19.1 s: mu = 1.753375772, error = 9.37e-13
//     k 17, 44.9 s: mu = nan
//   29jan2022, double, Newton refine more often:
//     k 15, 13.7 s: mu = 1.766837674, error = 9.77e-15
//     k 16, 33.4 s: mu = 1.753375772, error = 4.29e-14
//     k 17, 86.2 s: mu = 1.736187979
//     k 18, 208 s:  mu = 1.726163785
//     k 19, 524 s:  mu = 1.712556954 (0.01 below Ewing-Schober)
//   29jan2022, double, track known zeros:
//     k 15, 13.1 s: mu = 1.766837674, error = 2.82e-14
//     k 16, 32.1 s: mu = 1.753375772, error = 3.55e-15
//     k 17, 84.4 s: mu = 1.736187979
//     k 18, 210 s:  mu = 1.726163785
//   29jan2022, double, compute twiddles via arb:
//     k 15, 13.1 s: mu = 1.766837674, error = 3.11e-15
//     k 16, 35.5 s: mu = 1.753375772, error = 3e-14
//     k 17, 86.9 s: mu = 1.736187979
//     k 18, 214 s:  mu = 1.726163785
//     k 19, 545 s:  mu = 1.712556954 (still 0.01 below Ewing-Schober, but same as above?)
//   31jan2022, exp2, expansion arithmetic works:
//     k 13, 26.1 s: mu = 1.80617888585652252315542830417844,  error = 1.36e-31
//     k 14, 63.8 s: mu = 1.78638971655198525707933106215495,  error < 2.07e-20
//     k 15, 150 s:  mu = 1.76683767419408303662440486811699,  error < 2.44e-20
//     k 16, 352 s:  mu = 1.753375772356926299663688533777756, error < 3.63e-20
//     k 17, 833 s:  mu = 1.73618797890571568324245734657981

// CUDA series history (on an A100)
//   4feb2022, cuda-double, first working version:
//     k 10, 11.0 s: mu = 1.895943075803316, error = 3.41e-15
//     k 11, 16.4 s: mu = 1.8546567767819579, error = 2.83e-15
//   4feb2022, cuda-double, mul/sqr base cases:
//     k 10,  8.0 s: mu = 1.8959430758033156, error = 3.85e-15
//     k 11, 12.4 s: mu = 1.8546567767819575, error = 3.27e-15
//   4feb2022, cuda-double, inv/exp base cases:
//     k 11, 10.1 s: mu = 1.8546567767819564, error = 4.38e-15
//     k 12, 15.0 s: mu = 1.834655732626151, error = 4.65e-15
//     k 13, 21.7 s: mu = 1.8061788858565135, error = 8.98e-15
//     k 14, 30.5 s: mu = 1.78638971655198, error = 5.28e-15
//     k 15, 41.7 s: mu = 1.7668376741940959, error = 1.28e-14
//     k 16, 56.1 s: mu = 1.7533757723569496, error = 2.33e-14
//     k 17, 74.0 s: mu = 1.7361879789057115, error ≈ 5e-15 (vs. exp2)
//     k 18, 96.3 s: mu = 1.7261637845417952
//     k 19, 124 s:  mu = 1.7125569540291936
//     k 20, 160 s:  mu = 1.7032798671348994, error ≲ 1e-3 (vs. Bittner 1M)
//     k 21, 213 s:  mu = 1.6933586065947914, error ≲ 1e-3 (vs. Bittner 2M)
//     k 22, 301 s:  mu = 1.6858651156374813, error ≲ 1e-3 (vs. Bittner 4M)

}  // namespace mandelbrot
