
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

"""Study of the table end points (for all configurations):

A) The limit of sigma for wavelength -> 0: convergence at very short
   wavelengths, and the difference to sigma at 0.01 Aa.

B) Extending the table from 100 Aa to 300 Aa: extra points and time, and the
   error of linear extrapolation beyond the table end (to 1000 Aa).
"""
import json
import time

import numpy as np
from NCrystalDev.misc import evaluate_query as ncquery
from xsmat import Mat, all_cfgs


def table(cfg, wlmax, tol=1e-3):
    t0 = time.perf_counter()
    r = ncquery(['filtertable', cfg, f'tol={tol}', f'wlmax={wlmax}'])
    return np.asarray(r['wl']), np.asarray(r['xs']), time.perf_counter() - t0


def extrap(wl_t, xs_t, wl):
    # Linear continuation of the last segment, clamped at >= 0:
    x0, x1, y0, y1 = wl_t[-2], wl_t[-1], xs_t[-2], xs_t[-1]
    return np.maximum(0.0, y1 + (wl - x1) * (y1 - y0) / (x1 - x0))


def relerr(a, b):
    return abs(a - b) / b if b > 0 else abs(a - b)


res = {}
for cfg in all_cfgs():
    try:
        m = Mat(cfg)
    except ValueError:
        continue
    r = {}
    # A) Short wavelengths:
    wls = np.array([1e-2, 1e-3, 1e-4, 1e-5, 1e-6])
    s = m.sigma(wls)
    r['short_sigma'] = s.tolist()
    r['short_abs'] = m.sigma_abs(wls).tolist()
    lim = s[-1]
    r['conv_1e-5_vs_1e-6'] = relerr(s[-2], lim)
    r['diff_0.01_vs_lim'] = relerr(s[0], lim)
    # B) Long wavelengths:
    wl100, xs100, t100 = table(cfg, 100.0)
    wl300, xs300, t300 = table(cfg, 300.0)
    r['npts'] = [len(wl100), len(wl300)]
    r['time'] = [t100, t300]
    for wl in (300.0, 1000.0):
        exact = float(m.sigma(wl))
        r[f'extrap100_{wl:g}'] = relerr(float(extrap(wl100, xs100, wl)), exact)
        if wl > 300.0:
            r[f'extrap300_{wl:g}'] = relerr(float(extrap(wl300, xs300, wl)), exact)
    # log-log slope of sigma at 100 and 300 Aa:
    for wl in (100.0, 300.0):
        a, b = m.sigma(wl), m.sigma(wl * 1.01)
        r[f'slope_{wl:g}'] = float(np.log(b / a) / np.log(1.01)) if a > 0 and b > 0 else None
    res[cfg] = r

with open('limits.json', 'w') as f:
    json.dump(res, f, indent=1)


def summ(key, idx=None):
    v = np.array([r[key] if idx is None else r[key][idx] for r in res.values()])
    worst = max(res, key=lambda c: res[c][key] if idx is None else res[c][key][idx])
    return f'median {np.median(v):.3g}, 90% {np.percentile(v, 90):.3g}, max {v.max():.3g} ({worst})'


print(f'{len(res)} configurations')
print('A) sigma(1e-5 Aa) vs sigma(1e-6 Aa):', summ('conv_1e-5_vs_1e-6'))
print('A) sigma(0.01 Aa) vs sigma(1e-6 Aa):', summ('diff_0.01_vs_lim'))
print('B) points, wlmax=100:', summ('npts', 0))
print('B) points, wlmax=300:', summ('npts', 1))
d = np.array([r['npts'][1] - r['npts'][0] for r in res.values()])
print(f'B) extra points for 300: median {np.median(d):g}, max {d.max()}')
print('B) time wlmax=300 [s]:', summ('time', 1))
print('B) extrap from 100 to 300:', summ('extrap100_300'))
print('B) extrap from 100 to 1000:', summ('extrap100_1000'))
print('B) extrap from 300 to 1000:', summ('extrap300_1000'))
for key in ('slope_100', 'slope_300'):
    v = [(r[key], c) for c, r in res.items() if r[key] is not None]
    lo = sorted(v)[:5]
    print(f'B) {key}: lowest', ', '.join(f'{c.replace("stdlib::", "")} {s:.2f}' for s, c in lo))
