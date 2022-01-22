// Known Mandelbrot area bounds
#pragma once

const int known_ks = 16+1;
const char* known_areas[known_ks] = {
  "3.1415926535897932385 +/- 3.74e-20",  // k 0, prec 2000
  "3.1415926535897932385 +/- 3.74e-20",  // k 1, prec 2000
  "2.6998061866787285643 +/- 3.84e-21",  // k 2, prec 2000
  "2.4636540388624509718 +/- 1.11e-20",  // k 3, prec 2000
  "2.3089221574735649437 +/- 1.86e-20",  // k 4, prec 2000
  "2.1863969605631909502 +/- 3.77e-20",  // k 5, prec 2000
  "2.1117843557698974143 +/- 3.66e-20",  // k 6, prec 2000
  "2.0290979917037465762 +/- 2.62e-20",  // k 7, prec 2000
  "1.9793859246155479841 +/- 4.66e-21",  // k 8, prec 2000
  "1.9277162229003024017 +/- 3.17e-20",  // k 9, prec 2000
  "1.8959430758033194112 +/- 3.71e-20",  // k 10, prec 2000
  "1.8546567767819607804 +/- 1.97e-20",  // k 11, prec 2000
  "1.8346557326261556582 +/- 3.03e-20",  // k 12, prec 2000
  "1.8061788858565225232 +/- 4.46e-20",  // k 13, prec 2000
  "1.7863897165519852571 +/- 2.07e-20",  // k 14, prec 2000
  "1.7668376741940830366 +/- 2.45e-20",  // k 15, prec 2000
  "1.7533757723569262997 +/- 3.64e-20",  // k 16, prec 2000, repeats 2, radius 1.1345e-477
};

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

// Series history:
//   22jan2022, Series<double>, fft_mul:
//     k 14, 39.5 s: mu = 1.786389717, error = 7.11e-15
