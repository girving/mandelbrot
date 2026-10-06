# Analysis of the hybrid area runs

The scripts behind the result in `notes/hybrid-area.md`. Run them from the repo root. They need numpy, and the
image scripts also need matplotlib and PIL.

- `results/prod-d15/`, `results/tail-d13/`: the production and direct-tail runs. Each has one `escape_tree`
  result per shard, the merged summary (`merge.txt`), and the `cluster/idle_shards.py` run definition
  (`run.json`, built at tag `area-run-7`).
- `tree_results.py` parses the shard results into A(k), D(k), and their variances.
- `final_tail.py` fits k·D(k) and extrapolates the tail past 2^32: the first estimate, ±2.7e-11.
- `tail_backtest.py` backtests that extrapolation at earlier cuts.
- `tail_direct.py` gives the final result. It combines the backtest-corrected extrapolation with the direct
  measurement to 2^36, giving μ = 1.506591883653(24).
- `lo/` measures the biases of Hsing Lo's membership test against ours. The compile lines are in each file.
- `splash.py` and `area_image.py` make the announcement images: the quadtree splash and the logo with μ(M).
