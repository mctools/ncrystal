
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

import numpy as np
import validate_cpp as V
import xsmat

variants=[('nd2',['ndense=2']),('nd100',['ndense=100']),('nd1000',['ndense=1000']),('auto',[])]
acc={n:[] for n,_ in variants}
for cfg in xsmat.all_cfgs():
    m=xsmat.Mat(cfg); wlv=V.validation_points(m,0.01,100.0); sv=m.sigma(wlv)
    for n,o in variants: acc[n].append((cfg,V.run_one(cfg,o,m=m,wlv=wlv,sv=sv)))
for n,_ in variants:
    err=np.array([r['maxerr'] for _,r in acc[n]]); npts=np.array([r['npts'] for _,r in acc[n]]); t=np.array([r['t'] for _,r in acc[n]])*1e3
    nd=np.array([r['ndense_total'] for _,r in acc[n]])
    print(f"{n:6s} #>1.01tol={int((err>1.01e-3).sum()):3d} worst={err.max()/1e-3:8.3f}*tol npts sum={npts.sum()} evals med={np.median(nd):.0f} t med={np.median(t):.1f}ms max={t.max():.1f}ms")
    for c,r in sorted(acc[n],key=lambda cr:-cr[1]['maxerr'])[:3]:
        if r['maxerr']>1.01e-3: print('     ',c[:55],f"{r['maxerr']:.2e}")
