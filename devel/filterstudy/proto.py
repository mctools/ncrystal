
################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2026 NCrystal developers                                   ##
##                                                                            ##
##  Licensed under the Apache License, Version 2.0 (the "License");           ##
##  you may not use this file except in compliance with the License.          ##
##  You may obtain a copy of the License at                                   ##
##                                                                            ##
##      http://www.apache.org/licenses/LICENSE-2.0                            ##
##                                                                            ##
##  Unless required by applicable law or agreed to in writing, software       ##
##  distributed under the License is distributed on an "AS IS" BASIS,         ##
##  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  ##
##  See the License for the specific language governing permissions and       ##
##  limitations under the License.                                            ##
##                                                                            ##
################################################################################

"""Python prototypes of algorithms for piecewise-linear cross section tables.

A table is a set of knots (x_i, y_i) in some (X, Y) space. The cross section
at wavelength wl is obtained by linear interpolation in (X, Y) space, followed
by the inverse transform of Y. The quality measure is always the relative
error of sigma.
"""
import numpy as np

WL2EKIN = 0.081804209605330899  # E[eV] = WL2EKIN / wl[Aa]^2

# ---------------------------------------------------------------------------
# Spaces: functions (wl -> x), (x -> wl), (sigma -> y), (y -> sigma)
# ---------------------------------------------------------------------------
SPACES = {
    'wl_lin': (lambda wl: wl, lambda x: x, lambda s: s, lambda y: y),
    'ekin_lin': (lambda wl: WL2EKIN / wl**2, lambda x: np.sqrt(WL2EKIN / x),
                 lambda s: s, lambda y: y),
    'logwl_lin': (np.log, np.exp, lambda s: s, lambda y: y),
    'logwl_log': (np.log, np.exp, np.log, np.exp),
}


def interp_sigma(table, wl, space):
    """Evaluate a table (xk, yk; sorted in increasing x) at wavelengths wl."""
    fx, _, _, fyinv = SPACES[space]
    xk, yk = table
    x = fx(np.asarray(wl, dtype=float))
    # np.interp requires increasing xk and clamps outside (we only evaluate
    # inside the table range here).
    return fyinv(np.interp(x, xk, yk))


def _to_space(wl, s, space):
    fx, _, fy, _ = SPACES[space]
    x = fx(wl)
    y = fy(s)
    o = np.argsort(x, kind='stable')
    return x[o], y[o], s[o]


def _chord_relerr(x, y, s, i, j, fyinv):
    """Max relative sigma error of the chord (i,j) at the interior points."""
    if j - i < 2:
        return 0.0, -1
    xi, xj = x[i], x[j]
    t = (x[i + 1:j] - xi) / (xj - xi)
    yc = y[i] + t * (y[j] - y[i])
    sc = fyinv(yc)
    se = s[i + 1:j]
    with np.errstate(divide='ignore', invalid='ignore'):
        err = np.abs(sc - se) / np.where(se > 0, se, np.inf)
        err = np.where(se > 0, err, np.abs(sc - se) * 1e300)  # sigma=0: any deviation is infinite
    k = int(np.argmax(err))
    return float(err[k]), i + 1 + k


def reduce_dp(wl, s, space, tol):
    """Douglas-Peucker style reduction of a dense sample: recursively split
    at the point of largest relative error, until all points are within tol.
    Knots are a subset of the sample points."""
    x, y, ss = _to_space(wl, s, space)
    fyinv = SPACES[space][3]
    keep = np.zeros(len(x), dtype=bool)
    keep[0] = keep[-1] = True
    stack = [(0, len(x) - 1)]
    while stack:
        i, j = stack.pop()
        err, k = _chord_relerr(x, y, ss, i, j, fyinv)
        if err > tol:
            keep[k] = True
            stack.append((i, k))
            stack.append((k, j))
    return x[keep], y[keep]


def reduce_greedy(wl, s, space, tol):
    """Greedy forward reduction: from the current knot, extend the segment as
    far as possible (exponential + binary search on the end index, assuming
    approximately monotonic feasibility), then verify."""
    x, y, ss = _to_space(wl, s, space)
    fyinv = SPACES[space][3]
    n = len(x)
    knots = [0]
    i = 0
    while i < n - 1:
        # exponential search for an infeasible end
        step = 1
        good = i + 1
        while True:
            j = min(i + step, n - 1)
            err, _ = _chord_relerr(x, y, ss, i, j, fyinv)
            if err <= tol:
                good = j
                if j == n - 1:
                    break
                step *= 2
            else:
                break
        if good < n - 1 and step > 1:
            lo, hi = good, min(i + step, n - 1)
            while hi - lo > 1:
                mid = (lo + hi) // 2
                err, _ = _chord_relerr(x, y, ss, i, mid, fyinv)
                if err <= tol:
                    lo = mid
                else:
                    hi = mid
            good = lo
        knots.append(good)
        i = good
    knots = np.asarray(knots)
    return x[knots], y[knots]


def adaptive_bisection(mat, wlmin, wlmax, space, tol, nseed=100,
                       edges=None, edge_pad=1e-9, maxiter=60, ncheck=1):
    """Adaptive refinement with direct cross section evaluations (as done by
    e.g. NJOY when linearising cross sections): start from a log-spaced seed
    grid plus points just below and above each Bragg edge, and split each
    interval (in x space) at its midpoint until the table reproduces sigma
    within tol at the midpoint (and optionally at ncheck>1 interior points)."""
    fx, fxinv, fy, fyinv = SPACES[space]
    wl = np.geomspace(wlmin, wlmax, nseed)
    if edges is not None and len(edges):
        e = edges[(edges > wlmin) & (edges < wlmax)]
        wl = np.concatenate([wl, e * (1 - edge_pad), e * (1 + edge_pad)])
    x = np.unique(fx(wl))
    s = mat.sigma(fxinv(x))
    nevals = len(x)
    for _ in range(maxiter):
        # candidate check points in each interval
        fr = np.arange(1, ncheck + 1) / (ncheck + 1)
        xa, xb = x[:-1], x[1:]
        # do not split the tiny intervals across edges
        splittable = (xb - xa) > 1e-7 * np.maximum(np.abs(xa), np.abs(xb))
        xm = xa[:, None] + (xb - xa)[:, None] * fr[None, :]
        sm = mat.sigma(fxinv(xm.ravel())).reshape(xm.shape)
        nevals += sm.size
        ya, yb = fy(s[:-1]), fy(s[1:])
        sc = fyinv(ya[:, None] + (yb - ya)[:, None] * fr[None, :])
        with np.errstate(divide='ignore', invalid='ignore'):
            err = np.where(sm > 0, np.abs(sc - sm) / sm, np.where(np.abs(sc - sm) > 0, np.inf, 0.0))
        bad = (err.max(axis=1) > tol) & splittable
        if not bad.any():
            break
        # add all check points of bad intervals (keeps evaluations useful)
        newx = xm[bad].ravel()
        news = sm[bad].ravel()
        x = np.concatenate([x, newx])
        s = np.concatenate([s, news])
        o = np.argsort(x, kind='stable')
        x, s = x[o], s[o]
    return (x, fy(s)), nevals


def validate(mat, table, space, wlmin, wlmax, nrand=300000, seed=123, edges=None):
    """Max and 99.9% quantile of the relative error on an independent sample:
    random log-uniform wavelengths plus points near each Bragg edge (which are
    the places where tables are most likely to be wrong)."""
    rng = np.random.default_rng(seed)
    wl = np.exp(rng.uniform(np.log(wlmin), np.log(wlmax), nrand))
    if edges is not None and len(edges):
        e = edges[(edges > wlmin * (1 + 1e-6)) & (edges < wlmax * (1 - 1e-6))]
        for f in (1 - 1e-6, 1 + 1e-6, 1 - 1e-4, 1 + 1e-4):
            wl = np.concatenate([wl, e * f])
    s = mat.sigma(wl)
    st = interp_sigma(table, wl, space)
    ok = s > 0
    rel = np.abs(st[ok] - s[ok]) / s[ok]
    return float(rel.max()), float(np.quantile(rel, 0.999)), float(np.mean(rel))
