
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

from NCrystalDev.exceptions import NCException
from NCrystalDev.misc import evaluate_query as q


def show(qq):
    t0=time.time()
    try:
        r=q(qq); t=time.time()-t0
        print(f"OK  {qq[1:]}: npts={r['npts']} nedges={r['nedges']} nrefine={r['nrefine']} ndense_total={r['ndense_total']} t={t*1e3:.0f}ms xs[0]={r['xs'][0]:.4g} xs[-1]={r['xs'][-1]:.4g}")
    except NCException as e:
        print(f"ERR {qq[1:]}: {str(e)[:170]}")
A='stdlib::Al_sg225.ncmat'
show(['filtertable','stdlib::Al_sg225.ncmat;mos=0.3deg;dir1=@crys_hkl:0,0,1@lab:0,0,1;dir2=@crys_hkl:0,1,0@lab:0,1,0'])
show(['filtertable','stdlib::void.ncmat'])
show(['filtertable',A,'wlmin=1e-4','wlmax=1e4'])
show(['filtertable',A,'wlmin=1','wlmax=1.001'])
show(['filtertable',A,'wlmin=4.6','wlmax=4.7'])   # contains the Bragg cutoff 4.676
show(['filtertable',A,'wlmin=10','wlmax=20'])
show(['filtertable',A,'tol=1e-6'])
show(['filtertable','stdlib::CaSiO3_sg2_Wollastonite.ncmat','tol=1e-5'])
show(['filtertable',A,'tol=0'])
show(['filtertable',A,'wlmin=0'])
show(['filtertable',A,'wlmin=5','wlmax=1'])
show(['filtertable',A,'ndense=1'])
show(['filtertable',A,'ndense=2'])
show(['filtertable',A,'algo=foo'])
show(['filtertable',A,'tol'])
show(['filtertable'])
show(['filtertable','stdlib::Al_sg225.ncmat;temp=10K','tol=1e-3'])
show(['filtertable','phases<0.5*stdlib::Al_sg225.ncmat&0.5*stdlib::void.ncmat>'])
