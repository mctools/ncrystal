
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
import validate_cpp as V
import xsmat

tols=[1e-2,1e-3,1e-4]
res=[]
for cfg in xsmat.all_cfgs():
    m=xsmat.Mat(cfg); wlv=V.validation_points(m,0.01,100.0); sv=m.sigma(wlv)
    r={'cfg': cfg}
    for tol in tols: r[str(tol)]=V.run_one(cfg,[f'tol={tol:g}'],m=m,wlv=wlv,sv=sv)
    res.append(r)
with open('final_val.json','w') as f:
    json.dump(res,f,indent=1)
for tol in tols:
    k=str(tol); err=np.array([r[k]['maxerr'] for r in res]); npts=np.array([r[k]['npts'] for r in res]); t=np.array([r[k]['t'] for r in res])*1e3
    print(f"tol={tol:g}: worst={err.max()/tol:.3f}*tol, #>1.01tol={int((err>1.01*tol).sum())}, npts median={np.median(npts):.0f} 90%={np.quantile(npts,0.9):.0f} max={npts.max()}, time median={np.median(t):.1f}ms max={t.max():.1f}ms")
