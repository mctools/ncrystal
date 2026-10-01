
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

import time

import proto
import xsmat

cfgs = ['stdlib::Al_sg225.ncmat','stdlib::Be_sg194.ncmat;temp=80K','stdlib::Polyethylene_CH2.ncmat',
        'stdlib::LiquidWaterH2O_T293.6K.ncmat','stdlib::Y2SiO5_sg15_YSO.ncmat','stdlib::B4C_sg166_BoronCarbide.ncmat',
        'gasmix::air','stdlib::C_sg194_pyrolytic_graphite.ncmat','stdlib::CaSiO3_sg2_Wollastonite.ncmat']
WLMIN,WLMAX,TOL=0.01,100.0,1e-3
print(f"{'material':34s} {'method':22s} {'npts':>6s} {'maxerr':>9s} {'q999':>9s} {'nevals':>8s} {'t[s]':>6s}")
for cfg in cfgs:
    m = xsmat.Mat(cfg); e = m.edges(WLMIN,WLMAX)
    wl,s = xsmat.dense_reference(m,WLMIN,WLMAX,n=200000)
    for meth,space in [('dp','wl_lin'),('dp','ekin_lin'),('dp','logwl_lin'),('dp','logwl_log'),('greedy','wl_lin'),('bisect','wl_lin'),('bisect','ekin_lin'),('bisect','logwl_log')]:
        t0=time.time()
        if meth=='dp': tab=proto.reduce_dp(wl,s,space,TOL); nev=len(wl)
        elif meth=='greedy': tab=proto.reduce_greedy(wl,s,space,TOL); nev=len(wl)
        else: tab,nev=proto.adaptive_bisection(m,WLMIN,WLMAX,space,TOL,edges=e)
        t=time.time()-t0
        mx,q,_=proto.validate(m,tab,space,WLMIN,WLMAX,edges=e)
        print(f"{cfg.replace('stdlib::','')[:34]:34s} {meth+':'+space:22s} {len(tab[0]):6d} {mx:9.2e} {q:9.2e} {nev:8d} {t:6.2f}",flush=True)
