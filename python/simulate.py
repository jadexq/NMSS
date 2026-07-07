"""
Simulate network-mediator data for NMSS and write it in the file layout the
C++ sampler (Mcmc.cpp::readData, setting=1) expects.

Model (matching the C++ / manuscript):
  mediator : m_iv = alpha_v x_i + a1_v c11_i + a2_v c12_i + eps_iv,  eps ~ N(0, sigE)
  outcome  : y_i  = cV * sum_v beta_v m_iv + gamma x_i + C_i . b + delta_i,  delta ~ N(0, sigD)

Mediators are the E = V(V-1)/2 edges of a network on V nodes. The active
(true) mediators form a structured subnetwork (community / scale-free /
small-world). The triangular-clique neighbor structure used by the NMSS
smoothing prior is built by build_edge_neighbors.
"""
import os
import numpy as np


# --------------------------------------------------------------------------
# neighbor structure: triangular cliques over the edge set
# --------------------------------------------------------------------------
def build_edge_neighbors(V):
    """Edge indexing + triangular-clique neighbors.

    For edge e=(u,w), each third node t gives a triangle whose other two
    edges are (u,t) and (w,t); these are stored as nbhi[e,k], nbhj[e,k].
    Every edge has exactly V-2 such triangles, so mNnbh = V-2 (no padding).
    """
    pair2edge = -np.ones((V, V), dtype=np.int64)
    edge2pair = []
    e = 0
    for u in range(V):
        for w in range(u + 1, V):
            pair2edge[u, w] = e
            pair2edge[w, u] = e
            edge2pair.append((u, w))
            e += 1
    E = e
    edge2pair = np.asarray(edge2pair, dtype=np.int64)
    mNnbh = V - 2
    nbhi = np.empty((E, mNnbh), dtype=np.int64)
    nbhj = np.empty((E, mNnbh), dtype=np.int64)
    nodes = np.arange(V)
    for e in range(E):
        u, w = edge2pair[e]
        ts = nodes[(nodes != u) & (nodes != w)]
        nbhi[e] = pair2edge[u, ts]
        nbhj[e] = pair2edge[w, ts]
    return dict(E=E, V=V, mNnbh=mNnbh, edge2pair=edge2pair,
                pair2edge=pair2edge, nbhi=nbhi, nbhj=nbhj)


# --------------------------------------------------------------------------
# active-edge sets for the three network structures
# --------------------------------------------------------------------------
def _edges_community(V, pair2edge, rng, comm_size=10, n_comm=1):
    """Planted communities: all within-community edges are active."""
    active = []
    nodes = rng.permutation(V)
    for c in range(n_comm):
        block = nodes[c * comm_size:(c + 1) * comm_size]
        for a in range(len(block)):
            for b in range(a + 1, len(block)):
                active.append(pair2edge[block[a], block[b]])
    return np.unique(np.asarray(active, dtype=np.int64))


def _edges_scale_free(V, pair2edge, rng, m=2, m0=3):
    """Barabasi-Albert preferential attachment; active edges = its edges."""
    deg = np.zeros(V, dtype=np.int64)
    active = []
    for a in range(m0):                       # initial clique
        for b in range(a + 1, m0):
            active.append(pair2edge[a, b]); deg[a] += 1; deg[b] += 1
    for new in range(m0, V):
        cand = np.arange(new)
        w = deg[cand].astype(float) + 1e-9
        w /= w.sum()
        targets = rng.choice(cand, size=min(m, new), replace=False, p=w)
        for t in targets:
            active.append(pair2edge[new, t]); deg[new] += 1; deg[t] += 1
    return np.unique(np.asarray(active, dtype=np.int64))


def _edges_small_world(V, pair2edge, rng, k=4, p=0.1):
    """Watts-Strogatz ring lattice with rewiring; active edges = its edges."""
    active = set()
    for u in range(V):
        for j in range(1, k // 2 + 1):
            w = (u + j) % V
            if rng.random() < p:              # rewire
                w = rng.integers(V)
                while w == u:
                    w = rng.integers(V)
            if u != w:
                active.add(pair2edge[u, w])
    return np.unique(np.asarray(sorted(active), dtype=np.int64))


_STRUCTS = {'community': _edges_community,
            'scale_free': _edges_scale_free,
            'small_world': _edges_small_world}


def simulate_truth(nb, structure='community', seed=1, amp=(0.4, 0.8),
                   opposite_sign=True, **struct_kw):
    """Plant a structured sparse (alpha, beta) truth on the edge set.

    Active edges get alpha>0 and beta (opposite sign by default, so the
    per-edge mediation contribution alpha*beta is consistently signed).
    """
    rng = np.random.default_rng(seed)
    V, E = nb['V'], nb['E']
    active = _STRUCTS[structure](V, nb['pair2edge'], rng, **struct_kw)
    alpha = np.zeros(E)
    beta = np.zeros(E)
    lo, hi = amp
    alpha[active] = rng.uniform(lo, hi, size=active.size)
    bsign = -1.0 if opposite_sign else 1.0
    beta[active] = bsign * rng.uniform(lo, hi, size=active.size)
    med = (alpha != 0) & (beta != 0)
    return alpha, beta, med


def simulate_data(alpha, beta, nb, n, seed=7, cV=1.0, sigE=0.5, sigD=1.0,
                  gamma=-0.5, py=2, a1=0.0, a2=0.0, b=None):
    """Draw one dataset from the mediation model given a planted truth."""
    rng = np.random.default_rng(seed)
    E = alpha.shape[0]
    if b is None:
        b = np.zeros(py)
    x = rng.standard_normal(n)
    c11 = rng.standard_normal(n)                     # age (mediator confounder)
    c12 = (rng.random(n) < 0.5).astype(float)        # sex (mediator confounder)
    C = rng.standard_normal((n, py))                 # outcome confounders
    eps = rng.standard_normal((E, n)) * np.sqrt(sigE)
    M = np.outer(alpha, x) + a1 * np.outer(np.ones(E), c11) \
        + a2 * np.outer(np.ones(E), c12) + eps
    ybar = cV * (beta @ M) + gamma * x + C @ np.asarray(b)
    y = ybar + rng.standard_normal(n) * np.sqrt(sigD)
    return dict(x=x, M=M, y=y, c1=c11, c12=c12, C=C,
                gamma=gamma, cV=cV, sigE=sigE, sigD=sigD, b=np.asarray(b))


def ie_true(alpha, beta, cV=1.0):
    return cV * float(np.sum(alpha * beta))


# --------------------------------------------------------------------------
# write the C++ input layout
# --------------------------------------------------------------------------
def write_cpp_inputs(data, nb, datadir):
    """Write the 10 files Mcmc.cpp::readData reads under <datadir>."""
    os.makedirs(datadir, exist_ok=True)
    x, M, y = data['x'], data['M'], data['y']
    c11 = data['c1']
    c12 = data.get('c12')
    C = data['C']
    n, E, py = x.shape[0], M.shape[0], C.shape[1]

    def wcol(fn, arr):
        np.savetxt(os.path.join(datadir, fn),
                   np.asarray(arr).reshape(-1, 1), fmt='%.10g')

    wcol('it_std.txt', x)                                   # exposure
    wcol('age.txt', c11)                                    # c11
    wcol('sex.txt', c12 if c12 is not None else np.zeros(n))  # c12
    np.savetxt(os.path.join(datadir, 'c3.txt'), C, fmt='%.10g')       # c2 (n x py)
    np.savetxt(os.path.join(datadir, 'mtC.txt'), M.T, fmt='%.10g')    # subj x edges
    wcol('y.txt', y)                                        # outcome
    np.savetxt(os.path.join(datadir, 'nbhi.txt'), nb['nbhi'], fmt='%d')
    np.savetxt(os.path.join(datadir, 'nbhj.txt'), nb['nbhj'], fmt='%d')
    wcol('idx_dis80_T2.txt', np.zeros(E, dtype=int))        # 0 = important edge
    wcol('validlocC.txt', np.arange(E, dtype=int))          # update all edges
