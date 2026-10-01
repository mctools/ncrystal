
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
import proto
import xsmat

A,B=0.1,20.0
cfgs=[c for c in xsmat.all_cfgs() if 'void' not in c]
def extrap(tab,wl,mode,side):
    x,y=tab
    if side=='hi': x0,x1,y0,y1=x[-2],x[-1],y[-2],y[-1]
    else: x0,x1,y0,y1=x[0],x[1],y[0],y[1]
    if mode=='clamp': return np.full_like(wl, y1 if side=='hi' else y0)
    if mode=='linear': return y0+(wl-x0)*(y1-y0)/(x1-x0)
    if mode=='power':
        k=np.log(y1/y0)/np.log(x1/x0); xr,yr=(x1,y1) if side=='hi' else (x0,y0)
        return yr*(wl/xr)**k
rows=[]
for cfg in cfgs:
    m=xsmat.Mat(cfg)
    wl,s=xsmat.dense_reference(m,A,B,n=20000)
    tab=proto.reduce_greedy(wl,s,'wl_lin',1e-3)
    r=[cfg]
    for side,pts in (('hi',[40.,100.]),('lo',[0.05,0.01])):
        p=np.array(pts); sv=m.sigma(p)
        for mode in ('clamp','linear','power'):
            with np.errstate(all='ignore'):
                r.append(np.abs(extrap(tab,p,mode,side)/sv-1))
    rows.append(r)
def summ(idx,label):
    a=np.array([r[idx] for r in rows]);
    print(f"{label:22s} median={np.median(a,axis=0)}  90%={np.quantile(a,0.9,axis=0)}  max={a.max(axis=0)}")
names=['hi clamp','hi linear','hi power','lo clamp','lo linear','lo power']
print('relative error of extrapolated sigma at (2x,5x) beyond range [hi: 40,100 Aa; lo: 0.05, 0.01 Aa]')
for i,n in enumerate(names): summ(i+1,n)
print('\nworst cases for hi linear @100Aa and lo linear @0.01Aa:')
for idx,k in ((2,1),(5,1)):
    w=sorted(rows,key=lambda r:-r[idx][k])[:6]
    for r in w: print(f"  {names[idx-1]:10s} {r[0][:55]:55s} {r[idx][k]:.3f}")
