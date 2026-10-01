
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

"""Validation of the filter tables with the default range [0,500] Aa, for all
configurations: success (or the exception), number of points, time, the leaf
physics processes, and the worst relative error vs. the exact cross section at
random wavelengths (log-uniform in [1e-10,500] Aa, and uniform in [0,1e-5],
[0,1e-7] and [0,1e-9] Aa) and next to the Bragg edges. Also checks that materials with
physics processes not validated for filter tables give a warning (and a valid
table), and that oriented materials are rejected.

Usage: ./run.sh validate_full.py [TOL]
"""
import json
import sys
import time

import numpy as np
from NCrystalDev.misc import evaluate_query as ncquery
from xsmat import Mat, all_cfgs


def table_eval(wl_t, xs_t, wl):
    # Linear interpolation (as evalTable in C++, within the table range),
    # where a pair of identical wavelengths uses the second value for wl >= it:
    i = np.clip(np.searchsorted(wl_t, wl, side='right') - 1, 0, len(wl_t) - 2)
    x0, x1, y0, y1 = wl_t[i], wl_t[i + 1], xs_t[i], xs_t[i + 1]
    t = np.where(x1 > x0, (wl - x0) / np.where(x1 > x0, x1 - x0, 1.0), 0.0)
    return y0 + t * (y1 - y0)


tol = float(sys.argv[1]) if len(sys.argv) > 1 else 1e-3
rng = np.random.default_rng(12345)
wl_rand = np.concatenate([10**rng.uniform(-10, np.log10(500.0), 250000),
                          rng.uniform(0.0, 1e-5, 20000),
                          rng.uniform(0.0, 1e-7, 15000),
                          rng.uniform(0.0, 1e-9, 15000)])

res = {}
for cfg in all_cfgs():
    m = Mat(cfg)
    t0 = time.perf_counter()
    try:
        r = ncquery(['filtertable', cfg, f'tol={tol}'])
    except Exception as e:  # noqa: BLE001 (report all failures)
        res[cfg] = {'error': str(e)}
        print('FAILED', cfg, e, flush=True)
        continue
    dt = time.perf_counter() - t0
    wl, xs = np.asarray(r['wl']), np.asarray(r['xs'])
    e = m.edges(1e-9, 500.0)
    wl_test = np.concatenate([wl_rand, e * (1 - 1e-8), e * (1 + 1e-8),
                              e * (1 - 1e-5), e * (1 + 1e-5), [0.0, 500.0]])
    exact = m.sigma(wl_test)
    exact[wl_test == 0.0] = xs[0]  # (the exact value at 0 is the limit)
    approx = table_eval(wl, xs, wl_test)
    # The tolerance is relative to max(xs, 1e-12 barn):
    relerr = np.abs(approx - exact) / np.maximum(exact, 1e-12)
    res[cfg] = {'npts': len(wl), 'time': dt, 'worst': float(relerr.max() / tol),
                'processes': r['processes'], 'nrefine': r['nrefine'],
                'ndense': r['ndense_total'], 'xs0': float(xs[0])}

with open(f'validate_full_{tol:g}.json', 'w') as f:
    json.dump(res, f, indent=1)

ok = {c: r for c, r in res.items() if 'error' not in r}
print(f'tol={tol:g}: {len(ok)} of {len(res)} configurations OK')
for key, fmt in (('npts', '{:.0f}'), ('time', '{:.3f}'), ('worst', '{:.4f}')):
    v = np.array([r[key] for r in ok.values()])
    worst = max(ok, key=lambda c: ok[c][key])
    print(f'  {key}: median {fmt.format(np.median(v))}, 90% {fmt.format(np.percentile(v, 90))},'
          f' max {fmt.format(v.max())} ({worst})')
procs = sorted({p for r in ok.values() for p in r['processes']})
print('  leaf processes:', ', '.join(procs))

# Materials with physics processes not validated for filter tables (which give
# a warning, but valid tables):
for cfg in ['stdlib::Polyethylene_CH2.ncmat;ucnmode=remove',
            'stdlib::Polyethylene_CH2.ncmat;ucnmode=only']:
    r = ncquery(['filtertable', cfg, f'tol={tol}'])
    wl, xs = np.asarray(r['wl']), np.asarray(r['xs'])
    exact = Mat(cfg).sigma(wl_rand)
    worst = (np.abs(table_eval(wl, xs, wl_rand) - exact)
             / np.maximum(exact, 1e-12)).max() / tol
    print(f'  {cfg}: npts={len(wl)} worst={worst:.4f} (processes: {", ".join(r["processes"])})')

# Oriented materials must be rejected:
cfg = ('stdlib::Al_sg225.ncmat;dir1=@crys_hkl:0,0,1@lab:0,0,1;'
       'dir2=@crys_hkl:0,1,0@lab:0,1,0;mos=1deg')
try:
    ncquery(['filtertable', cfg])
    print('  NOT REJECTED:', cfg)
except Exception as e:  # noqa: BLE001 (any exception is a rejection)
    print('  rejected as expected:', cfg, '->', str(e)[:150])
