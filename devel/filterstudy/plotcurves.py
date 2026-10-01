
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

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import xsmat

cfgs = ['stdlib::Al_sg225.ncmat','stdlib::Be_sg194.ncmat;temp=80K','stdlib::Polyethylene_CH2.ncmat',
        'stdlib::LiquidWaterH2O_T293.6K.ncmat','stdlib::Y2SiO5_sg15_YSO.ncmat','stdlib::B4C_sg166_BoronCarbide.ncmat',
        'gasmix::air','stdlib::V_sg229.ncmat','stdlib::C_sg194_pyrolytic_graphite.ncmat']
fig, axs = plt.subplots(3,3,figsize=(15,12))
for ax,cfg in zip(axs.flat,cfgs):
    m = xsmat.Mat(cfg)
    wl,s = xsmat.dense_reference(m,0.01,100,n=50000)
    ax.loglog(wl,s,lw=0.7,label='total')
    ax.loglog(wl,m.sigma_abs(wl),lw=0.7,label='abs')
    ax.set_title(cfg.replace('stdlib::',''),fontsize=9); ax.grid(alpha=0.3)
    ax.set_xlabel('wavelength [Aa]'); ax.set_ylabel('barn/atom')
axs.flat[0].legend()
fig.tight_layout(); fig.savefig('curves.png',dpi=70)
