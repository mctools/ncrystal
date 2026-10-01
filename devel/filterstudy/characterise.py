
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

"""Characterise total cross section curves of many materials (task 1)."""
import json
import sys
import time

import numpy as np
import xsmat

WLMIN, WLMAX = 0.01, 100.0


def loglog_slope(mat, wl):
    a, b = wl * (1 - 1e-3), wl * (1 + 1e-3)
    sa, sb = mat.sigma([a, b])
    return float(np.log(sb / sa) / np.log(b / a))


def main():
    cfgs = xsmat.all_cfgs()
    out = []
    for cfg in cfgs:
        try:
            t0 = time.time()
            m = xsmat.Mat(cfg)
            tinit = time.time() - t0
        except (ValueError, xsmat.NC.NCException) as ex:
            print(f'SKIP {cfg}: {ex}', file=sys.stderr)
            continue
        t0 = time.time()
        wl, s = xsmat.dense_reference(m, WLMIN, WLMAX, n=100000)
        tdense = time.time() - t0
        e = m.edges(WLMIN, WLMAX)
        jumps = m.edge_jumps(e) if len(e) else np.zeros(0)
        # relative jump: relative to sigma just above the edge
        rel = jumps / m.sigma(e * (1 + 1e-9)) if len(e) else np.zeros(0)
        # smoothness of the non-edge parts: max |second log-derivative|
        # on the regular log grid, away from edges
        g = np.geomspace(WLMIN, WLMAX, 20001)
        sg = m.sigma(g)
        d2 = np.abs(np.diff(np.log(sg), 2))
        r = {'cfg': cfg, 'tinit': tinit, 'tdense': tdense, 'npts_dense': len(wl),
                 'smin': float(s.min()), 'smax': float(s.max()),
                 's_at': {str(x): float(m.sigma([x])[0]) for x in (0.01, 0.1, 1.0, 10.0, 100.0)},
                 'slope_lo': loglog_slope(m, 0.012), 'slope_hi': loglog_slope(m, 90.0),
                 'nedges': len(e),
                 'nedges_rel_gt': {str(t): int((rel > t).sum()) for t in (1e-5, 1e-4, 1e-3, 1e-2, 1e-1)},
                 'max_rel_jump': float(rel.max()) if len(rel) else 0.0,
                 'max_edge': float(e.max()) if len(e) else 0.0,
                 'frac_abs_at_1AA': float(m.sigma_abs([1.0])[0] / m.sigma([1.0])[0]),
                 'max_d2_loggrid': float(d2.max())}
        out.append(r)
        print(f"{cfg[:60]:60s} edges={r['nedges']:6d} >1e-3:{r['nedges_rel_gt']['0.001']:5d} "
              f"s=[{r['smin']:.3g},{r['smax']:.3g}] slopes={r['slope_lo']:+.2f}/{r['slope_hi']:+.2f} "
              f"tinit={tinit:.2f}s", flush=True)
    with open('characterise.json', 'w') as f:
        json.dump(out, f, indent=1)


if __name__ == '__main__':
    main()
