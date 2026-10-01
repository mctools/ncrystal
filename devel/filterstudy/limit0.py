
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

"""Study of a table point at wavelength 0 (for all configurations).

The limit for wavelength -> 0 is estimated from the total cross section
(scattering + absorption, both evaluated with the NCrystal process objects),
by linear extrapolation to 0 from very short wavelengths, as done in the C++
code. Checks the convergence of this estimate, and the error of linear
interpolation between (0, limit) and the point at 0.01 Aa, compared to the
exact values in between.
"""
import numpy as np
from xsmat import Mat, all_cfgs

WL0 = 0.01
wl_between = np.geomspace(1e-5, WL0, 200)[:-1]

rows = []
for cfg in all_cfgs():
    try:
        m = Mat(cfg)
    except ValueError:
        continue
    def extrap0(h, m=m):
        a, b = m.sigma([h, 2 * h])
        return 2 * a - b
    lim_a, lim = extrap0(1e-7), extrap0(1e-8)
    conv = abs(lim_a - lim) / lim if lim > 0 else abs(lim_a - lim)
    y0 = float(m.sigma(WL0))
    interp = lim + (y0 - lim) * wl_between / WL0
    exact = m.sigma(wl_between)
    err = np.max(np.abs(interp - exact) / np.where(exact > 0, exact, 1.0))
    abs_frac = float(m.sigma_abs(WL0)) / y0 if y0 > 0 else 0.0
    rows.append((cfg, lim, conv, err, abs_frac, float(m.sigma_scat(WL0))))

rows.sort(key=lambda r: -r[3])
conv = np.array([r[2] for r in rows])
err = np.array([r[3] for r in rows])
print(f'{len(rows)} configurations')
print(f'limit convergence (estimates from 1e-7 and 1e-8 Aa): median {np.median(conv):.2g},'
      f' max {conv.max():.2g}')
print(f'linear interpolation error on [1e-5, 0.01] Aa: median {np.median(err):.2g},'
      f' 90% {np.percentile(err, 90):.2g}, max {err.max():.2g}')
print(f'configurations with error > 1e-3: {np.sum(err > 1e-3)}, > 1e-4: {np.sum(err > 1e-4)}')
print('worst 10 (cfg, limit [barn], interp error, absorption fraction at 0.01,'
      ' scattering at 0.01 [barn]):')
for cfg, lim, c, e, a, s in rows[:10]:
    print(f'  {cfg.replace("stdlib::", "")}: {lim:.4g}, {e:.2g}, {a:.2g}, {s:.4g}')
