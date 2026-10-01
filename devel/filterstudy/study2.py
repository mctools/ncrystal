
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

"""npts vs tolerance and dense-sample size, all materials, greedy in (wl, sigma)."""
import json
import sys
import time

import proto
import xsmat

WLMIN,WLMAX=0.01,100.0
res=[]
for cfg in xsmat.all_cfgs():
    try: m=xsmat.Mat(cfg)
    except (ValueError, xsmat.NC.NCException) as ex:
        print('SKIP',cfg,ex,file=sys.stderr)
        continue
    if m.sigma([1.0])[0] <= 0: print('SKIP zero xs',cfg,file=sys.stderr); continue
    e=m.edges(WLMIN,WLMAX)
    r={'cfg': cfg,'nedges': len(e)}
    for n in (5000,20000,200000):
        wl,s=xsmat.dense_reference(m,WLMIN,WLMAX,n=n)
        for tol in (1e-2,1e-3,1e-4):
            if n!=200000 and tol!=1e-3: continue
            t0=time.time(); tab=proto.reduce_greedy(wl,s,'wl_lin',tol); t=time.time()-t0
            mx,q,_=proto.validate(m,tab,'wl_lin',WLMIN,WLMAX,edges=e,nrand=200000)
            r[f'n{n}_tol{tol:g}']={'npts': len(tab[0]),'maxerr': mx,'q999': q,'t': t}
    res.append(r)
    a=r['n200000_tol0.001']; b=r['n5000_tol0.001']; c=r['n20000_tol0.001']
    print(f"{cfg.replace('stdlib::','')[:40]:40s} e={len(e):6d} "
          + ' '.join(f"{k.split('_')[1]}:{v['npts']:5d}/{v['maxerr']:.1e}" for k,v in r.items() if k.startswith('n200000'))
          + f" | n5k:{b['npts']:5d}/{b['maxerr']:.1e} n20k:{c['npts']:5d}/{c['maxerr']:.1e}", flush=True)
with open('study2.json','w') as f:
    json.dump(res,f,indent=1)
