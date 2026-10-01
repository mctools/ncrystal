
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

cfgs=['stdlib::Al_sg225.ncmat','stdlib::Be_sg194.ncmat;temp=80K','stdlib::Polyethylene_CH2.ncmat','stdlib::Y2SiO5_sg15_YSO.ncmat',
      'stdlib::CaSiO3_sg2_Wollastonite.ncmat','stdlib::C_sg194_pyrolytic_graphite.ncmat','gasmix::air','stdlib::Fe_sg229_Iron-alpha.ncmat']
print(f"{'material':34s} {'[0.01,100]':>10s} {'[0.1,20]':>9s} {'[0.5,10]':>9s} {'[1,6]':>6s} | edge-knots share @[0.01,100]")
for cfg in cfgs:
    m=xsmat.Mat(cfg); out=[]
    for a,b in ((0.01,100),(0.1,20),(0.5,10),(1,6)):
        wl,s=xsmat.dense_reference(m,a,b,n=20000)
        tab=proto.reduce_greedy(wl,s,'wl_lin',1e-3); out.append(len(tab[0]))
        if a==0.01:
            e=m.edges(a,b); x=tab[0]
            j=np.searchsorted(e,x); j=np.clip(j,1,len(e)-1) if len(e)>1 else np.zeros_like(j)
            near = np.minimum(np.abs(x-e[j]),np.abs(x-e[j-1]))/x < 1e-8 if len(e)>1 else np.zeros(len(x),bool)
            share=near.mean()
    print(f"{cfg.replace('stdlib::','')[:34]:34s} {out[0]:10d} {out[1]:9d} {out[2]:9d} {out[3]:6d} | {share:.2f}")
