
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

"""Validate the C++ filtertable query against exact cross sections (task 4)."""
import json
import time

import numpy as np
import xsmat
from NCrystalDev.misc import evaluate_query


def table_eval(wl_t, xs_t, wl):
    """Linear interpolation, where a pair of identical wavelengths (a Bragg
    edge) uses the second value (above the edge) for wl >= edge."""
    wl_t = np.asarray(wl_t)
    xs_t = np.asarray(xs_t)
    i = np.searchsorted(wl_t, wl, side='right') - 1
    i = np.clip(i, 0, len(wl_t) - 2)
    x0, x1 = wl_t[i], wl_t[i + 1]
    y0, y1 = xs_t[i], xs_t[i + 1]
    with np.errstate(divide='ignore', invalid='ignore'):
        t = np.where(x1 > x0, (wl - x0) / (x1 - x0), 0.0)
    return y0 + t * (y1 - y0)


def validation_points(m, wlmin, wlmax, nrand=300000, seed=4321):
    rng = np.random.default_rng(seed)
    wl = np.exp(rng.uniform(np.log(wlmin), np.log(wlmax), nrand))
    e = m.edges(wlmin * 1.001, wlmax / 1.001)
    if len(e):
        wl = np.concatenate([wl] + [e * f for f in (1 - 1e-8, 1 + 1e-8, 1 - 1e-5, 1 + 1e-5)])
        # exclude the tiny zones within 1e-9 of any edge (inside merged edge clusters)
        j = np.searchsorted(e, wl)
        jl = np.clip(j - 1, 0, len(e) - 1)
        jr = np.clip(j, 0, len(e) - 1)
        d = np.minimum(np.abs(wl - e[jl]), np.abs(wl - e[jr])) / wl
        wl = wl[d > 1e-9]
    return wl


def run_one(cfg, opts, wlmin=0.01, wlmax=100.0, m=None, wlv=None, sv=None):
    q = ['filtertable', cfg, f'wlmin={wlmin:g}', f'wlmax={wlmax:g}'] + opts
    evaluate_query(q)  # warm up (material loading and caches)
    t0 = time.time()
    r = evaluate_query(q)
    t = time.time() - t0
    st = table_eval(r['wl'], r['xs'], wlv)
    ok = sv > 0
    rel = np.abs(st[ok] - sv[ok]) / sv[ok]
    zero_ok = bool(np.all(st[~ok] == 0.0))
    return {'npts': r['npts'], 'nedges': r['nedges'], 't': t, 'nrefine': r['nrefine'], 'ndense_total': r['ndense_total'],
                'maxerr': float(rel.max()) if len(rel) else 0.0,
                'q9999': float(np.quantile(rel, 0.9999)) if len(rel) else 0.0,
                'zero_ok': zero_ok}


def main():
    variants = [('greedy_1e-2', ['tol=1e-2']),
                ('greedy_1e-3', ['tol=1e-3']),
                ('greedy_1e-4', ['tol=1e-4']),
                ('dp_1e-3', ['tol=1e-3', 'algo=dp']),
                ('noedges_1e-3', ['tol=1e-3', 'edges=0']),
                ('greedy_1e-3_nd5k', ['tol=1e-3', 'ndense=5000'])]
    res = []
    for cfg in xsmat.all_cfgs():
        m = xsmat.Mat(cfg)
        wlv = validation_points(m, 0.01, 100.0)
        sv = m.sigma(wlv)
        r = {'cfg': cfg}
        for name, opts in variants:
            r[name] = run_one(cfg, opts, m=m, wlv=wlv, sv=sv)
        res.append(r)
        print(f"{cfg.replace('stdlib::', '')[:38]:38s} "
              + ' '.join(f"{n}:{r[n]['npts']}/{r[n]['maxerr']:.1e}/{1e3*r[n]['t']:.0f}ms" for n, _ in variants),
              flush=True)
    with open('validate_cpp.json', 'w') as f:
        json.dump(res, f, indent=1)


if __name__ == '__main__':
    main()
