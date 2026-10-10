# Bulb areas as a weighted finite automaton

Research scripts for the renormalization cascade (notes/julia-area.md, "Bulb areas are a rational series").  They
read bulb_areas outputs from `$BULB_DATA` and run from this directory (they `exec` each other's preambles).

- `hankel_words.py`, `hankel20.py`: prefix/suffix sets and the bulbs a Hankel block H[u, s] = F(u·s) needs.
- `large_words.py`: bulbs for large digits (prefix rows u·b, b ≤ 400, and suffixes [c], c ≤ 400).
- `wfa.py`: Hankel singular values; spectral learning F(w) = αᵀ N(a1)⋯N(a_{k-1}) β(a_k); validation on long words.
- `wfa_sum.py`, `wfa_sum2.py K`: the automaton summed exactly by a matrix-valued Mayer resolvent, used as a control
  variate on the tail of the bulb sum (wfa_sum2 with measured large-digit matrices).
- `wfa_check.py`: large-digit predictions against F([n,2]), F([2,n]), F(1/n).
- `wfa_scale2.py`: held-out error against Hankel size and K.

Data (local, 10 threads): words.out 111k bulbs (14 min), hankel.out 18k (2 min), hankel20.out 37k (6 min),
large2.out 7k (4 min; 523 retried with BULB_TOL=1e-9).

High-precision pipeline (bulb_areas with BULB_EXP=1, double-double polish):
- `wfa_hp.py`: double-double data, suffix-closure learning (N(b) from H[u, b·s]), held-out error against K.
- `basis.py`: a well-conditioned 56 × 56 basis of short prefixes/suffixes for shifted blocks F(u·b·s).
- `wfa_hp_sum.py K ...`: exact head (all q ≤ 1000) + automaton control variate; N(4..256) from shifted blocks.
- `diag.py K`: model error on 400 < q ≤ 500 by word class; `knobs.py`: Chebyshev degree and digit cutoff.
Results: analysis/results/bulb-total.log.
