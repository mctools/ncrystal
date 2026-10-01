
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

m=xsmat.Mat('stdlib::Y_sg194.ncmat'); e=m.edges(0.01,100)
wl,s=xsmat.dense_reference(m,0.01,100,n=200000)
tab=proto.reduce_greedy(wl,s,'wl_lin',1e-4)
rng=np.random.default_rng(123)
w=np.exp(rng.uniform(np.log(0.01),np.log(100),200000))
ee=e[(e>0.01*(1+1e-6))&(e<100*(1-1e-6))]
for f in (1-1e-6,1+1e-6,1-1e-4,1+1e-4): w=np.concatenate([w,ee*f])
sv=m.sigma(w); st=proto.interp_sigma(tab,w,'wl_lin'); rel=np.abs(st-sv)/sv
i=np.argsort(rel)[-5:]
for k in i:
    x=w[k]; j=np.searchsorted(e,x); near=e[max(j-1,0):j+1]
    print(f"wl={x:.9f} rel={rel[k]:.2e} nearest edges={near} relfromedge={(x-near)/near}")
    # dense ref points around
    jj=np.searchsorted(wl,x); print('   dense nbrs', wl[jj-2:jj+2], s[jj-2:jj+2], 'sigma(x)=',sv[k], 'tab=',st[k])
