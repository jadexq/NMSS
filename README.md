# NMSS

Network Mediation Structural Shrinkage — a Bayesian high-dimensional
mediation model for network (connectivity) mediators.

Reference: [Imaging mediation analysis for longitudinal outcomes: a case study
of childhood brain tumor survivorship](https://scholar.google.com/citations?view_op=view_citation&hl=en&user=qcu4dJQAAAAJ&sortby=pubdate&citation_for_view=qcu4dJQAAAAJ:Se3iqnhoufwC)

## Layout

```
src/       C++ Gibbs sampler (Eigen)
python/    Python driver + simulation utilities
scripts/   get_eigen.sh (fetches Eigen 3.3.9 -> extern/, gitignored)
bin/       compiled binary (gitignored)
```

## Build

Requires **Eigen 3.3** (3.4 breaks the `vec(double)` indexing used here).

```bash
bash scripts/get_eigen.sh
g++ -O2 -I extern/eigen -std=c++14 src/*.cpp -o bin/nmss
```

## Run

The binary is a pure function of a working directory: it reads
`<workdir>/config.txt` (tuning + control) and `<workdir>/data/*`, infers all
dimensions from the data, and writes posterior summaries to
`<workdir>/output/*`. Drive it from Python:

```python
import sys; sys.path.insert(0, "python")
from simulate import build_edge_neighbors, simulate_truth, simulate_data
from nmss_cpp import run_cpp, cv_cpp

nb = build_edge_neighbors(50)
a, b, med = simulate_truth(nb, "community", seed=1)
data = simulate_data(a, b, nb, n=500, sigD=1.0)

res = run_cpp(data, nb, config=dict(sigDFix=1.0))          # single fit
grid = [dict(sigDFix=s) for s in (0.25, 0.5, 1.0, 2.0)]
cv = cv_cpp(data, nb, grid, K=5, n_jobs=-1)                # tune by CV
```

`run_cpp` returns posterior means (`alpha_mean`, `beta_mean`), inclusion
probabilities (`pip_alpha`, `pip_beta`), and effects (`IE`, `DE`). Each run is
isolated in its own directory, so `cv_cpp` parallelises across folds x grid.

## Config keys

`config.txt` holds one `key value` per line. Control: `nIter`, `nBurn`,
`nRep`, `setting`, `applyCv`, `cv`. Tuning: `aIA`, `bIA`, `aIB`, `bIB`, `cI`,
`sig1Alpha`, `sig1Beta`, `sigDFix`. Dimensions are **not** listed — they are
inferred from the data files (with cross-file consistency checks). Any key
omitted keeps its compiled default.

## Notes

- Data files (see `python/simulate.write_cpp_inputs`): `it_std.txt` (exposure),
  `age.txt`/`sex.txt` (mediator confounders), `c3.txt` (outcome confounders),
  `mtC.txt` (subjects x edges), `y.txt`, `nbhi.txt`/`nbhj.txt` (triangular-clique
  neighbors), `idx_dis80_T2.txt`, `validlocC.txt`.
- `sigDFix` is the fixed σ²_δ (the manuscript treats it as a CV-tuned
  parameter). The `updateB` β-update was corrected so the cross-edge term is
  weighted by `1/σ²_δ` per Supplementary Eq. (2); the sampler is now stable and
  accurate for any σ²_δ, not only σ²_δ = 1.
- RNG (`rnorm`/`rgamma`) reseeds from `random_device` per call, so runs are not
  bit-reproducible; agreement across runs is up to Monte-Carlo error.
```
