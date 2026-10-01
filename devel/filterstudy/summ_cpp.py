
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

import json

import numpy as np

with open('validate_cpp.json') as f:
    res=json.load(f)
names=[k for k in res[0] if k!='cfg']
tols={'greedy_1e-2':1e-2,'greedy_1e-3':1e-3,'greedy_1e-4':1e-4,'dp_1e-3':1e-3,'noedges_1e-3':1e-3,'greedy_1e-3_nd5k':1e-3}
print(f"{'variant':18s} {'npts med':>8s} {'npts max':>8s} {'sum npts':>8s} {'#viol(>1.01tol)':>15s} {'worst err':>9s} {'t med[ms]':>9s} {'t max[ms]':>9s}")
for n in names:
    npts=np.array([r[n]['npts'] for r in res]); err=np.array([r[n]['maxerr'] for r in res]); t=np.array([r[n]['t'] for r in res])*1e3
    viol=(err>1.01*tols[n]).sum()
    print(f"{n:18s} {np.median(npts):8.0f} {npts.max():8d} {npts.sum():8d} {viol:15d} {err.max():9.2e} {np.median(t):9.1f} {t.max():9.1f}")
    if viol:
        for r in sorted(res,key=lambda r:-r[n]['maxerr'])[:4]: print('     ',r['cfg'][:60],f"{r[n]['maxerr']:.2e}")
print('zero_ok all:', all(r[n]['zero_ok'] for r in res for n in names))
print()
for n in names:
    nr=np.array([r[n]['nrefine'] for r in res]); nd=np.array([r[n]['ndense_total'] for r in res])
    print(f"{n:18s} nrefine: median={np.median(nr):.0f} max={nr.max()}  ndense_total max={nd.max()}")
