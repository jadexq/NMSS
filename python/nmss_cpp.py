"""
Python driver for the NMSS C++ sampler.

The C++ binary is a pure function of a working directory: it reads
<workdir>/config.txt (tuning + control) and <workdir>/data/* (the input
files), infers all dimensions from the data, runs the Gibbs sampler, and
writes posterior summaries to <workdir>/output/*. This module writes those
inputs, invokes the binary in an isolated directory, and parses the outputs.

Because every run is self-contained in its own directory, K-fold CV over a
tuning grid parallelizes trivially (joblib) with no shared state.
"""
import os
import shutil
import subprocess
import tempfile

import numpy as np

from simulate import write_cpp_inputs

_HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(os.path.dirname(_HERE), 'bin', 'nmss')

# control + tuning defaults; override per call via `config=...`
DEFAULT_CONFIG = dict(
    nIter=2000, nBurn=1000, nRep=1, nPrint=500, setting=1,
    applyCv=1, cv=1.0,
    aIA=-4.0, aIA0=-4.0, bIA=0.02,
    aIB=-4.0, aIB0=-4.0, bIB=0.02, cI=1.0,
    sig1Alpha=1.0, sig1Beta=1.0, sigDFix=1.0,
)


def _write_config(path, cfg):
    with open(path, 'w') as f:
        for k, v in cfg.items():
            f.write(f"{k} {v}\n")


def run_cpp(data, nb, config=None, workdir=None, keep=False, bin=BIN):
    """Run the C++ sampler on one dataset; return posterior summaries.

    Returns a dict with alpha_mean, beta_mean, pip_alpha, pip_beta, b, gamma,
    IE, DE (and the workdir if kept). Estimates match run_nmss's shape.
    """
    cfg = {**DEFAULT_CONFIG, **(config or {})}
    tmp = workdir or tempfile.mkdtemp(prefix='nmss_')
    made_tmp = workdir is None
    try:
        os.makedirs(os.path.join(tmp, 'data'), exist_ok=True)
        os.makedirs(os.path.join(tmp, 'output'), exist_ok=True)
        write_cpp_inputs(data, nb, os.path.join(tmp, 'data'))
        _write_config(os.path.join(tmp, 'config.txt'), cfg)

        r = subprocess.run([bin, tmp], capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"nmss failed (rc={r.returncode})\n"
                               f"--- stdout ---\n{r.stdout}\n"
                               f"--- stderr ---\n{r.stderr}")

        out = os.path.join(tmp, 'output')

        def _load(fn):
            return np.atleast_1d(np.loadtxt(os.path.join(out, fn)).ravel())

        alpha = _load('rep_alpha.txt')
        beta = _load('rep_beta.txt')
        pipA = _load('rep_zetaA.txt')
        pipB = _load('rep_zetaB.txt')
        bvec = _load('rep_b.txt')
        gamma = float(_load('rep_gamma.txt')[0])
        cV = float(cfg['cv'])
        res = dict(alpha_mean=alpha, beta_mean=beta, pip_alpha=pipA,
                   pip_beta=pipB, b=bvec, gamma=gamma,
                   IE=cV * float(np.sum(alpha * beta)), DE=gamma,
                   stdout=r.stdout)
        if keep:
            res['workdir'] = tmp
        return res
    finally:
        if made_tmp and not keep:
            shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------------------
# cross-validation helpers
# --------------------------------------------------------------------------
def subset_data(data, idx):
    """Restrict a dataset to a subset of subjects (columns of M)."""
    d = dict(data)
    d['x'] = data['x'][idx]
    d['M'] = data['M'][:, idx]
    d['y'] = data['y'][idx]
    d['c1'] = data['c1'][idx]
    if data.get('c12') is not None:
        d['c12'] = data['c12'][idx]
    d['C'] = data['C'][idx]
    return d


def predict_y(res, data):
    """Predict y from fitted (beta, gamma, b) and observed mediators."""
    cV = data['cV']
    return (cV * (res['beta_mean'] @ data['M'])
            + res['gamma'] * data['x'] + data['C'] @ res['b'])


def cv_cpp(data, nb, grid, base_config=None, K=5, n_jobs=1, seed=0):
    """K-fold predictive-MSE cross-validation over a tuning grid.

    `grid` is a list of config-override dicts (one per tuning combination).
    Returns a list of dicts {config, cv_mse, fold_mse} sorted by cv_mse, so
    result[0] is the selected tuning.
    """
    from joblib import Parallel, delayed

    n = data['x'].shape[0]
    rng = np.random.default_rng(seed)
    folds = np.array_split(rng.permutation(n), K)

    jobs = []
    meta = []
    for gi, override in enumerate(grid):
        cfg = {**(base_config or {}), **override}
        for k in range(K):
            test = folds[k]
            train = np.concatenate([folds[j] for j in range(K) if j != k])
            jobs.append(delayed(_cv_one)(subset_data(data, train),
                                         subset_data(data, test), nb, cfg))
            meta.append((gi, k))

    mses = Parallel(n_jobs=n_jobs)(jobs)

    per_combo = [[] for _ in grid]
    for (gi, _k), mse in zip(meta, mses):
        per_combo[gi].append(mse)

    results = []
    for gi, override in enumerate(grid):
        fm = per_combo[gi]
        results.append(dict(config=override,
                            cv_mse=float(np.mean(fm)),
                            fold_mse=[float(v) for v in fm]))
    results.sort(key=lambda r: r['cv_mse'])
    return results


def _cv_one(train, test, nb, cfg):
    res = run_cpp(train, nb, config=cfg)
    yhat = predict_y(res, test)
    return float(np.mean((test['y'] - yhat) ** 2))


# --------------------------------------------------------------------------
# selection metrics
# --------------------------------------------------------------------------
def selection_metrics(pip, truth_mask, thresh=0.5):
    """TPR/FPR of a PIP-thresholded selection against a boolean truth mask."""
    sel = pip >= thresh
    t = np.asarray(truth_mask, dtype=bool)
    tp = int(np.sum(sel & t))
    fp = int(np.sum(sel & ~t))
    fn = int(np.sum(~sel & t))
    tn = int(np.sum(~sel & ~t))
    tpr = tp / (tp + fn) if (tp + fn) else 0.0
    fpr = fp / (fp + tn) if (fp + tn) else 0.0
    return dict(tpr=tpr, fpr=fpr, tp=tp, fp=fp, fn=fn, tn=tn)
