#!/usr/bin/env python3

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

# NEEDS: numpy

# Dedicated test for various plotting code, to increase test coverage.
import NCrystalDev.core as nccore
from NCrystalDev.ncmat import NCMATComposer
from NCrystalDev._numpy import _np_linspace
from NCTestUtils.env import ncsetenv

ncmatfile = 'Al_sg225.ncmat'

def extract_alvdosdata():
    info = nccore.createInfo(f'stdlib::{ncmatfile}')
    return info.dyninfos[0].vdosData()

def main(do_plot):
    if not do_plot:
        ncsetenv('FAKEPYPLOT','1')
    matbase = NCMATComposer()
    matbase.set_density(2.7,'g/cm3')
    matbase.set_dyninfo_vdos('Al',*extract_alvdosdata())

    wls = _np_linspace( 0.0, 6.0, 1000 if do_plot else 20 )
    if do_plot:
        import matplotlib.pyplot as plt

    # str(Gnmax) -> (wl,xs_min,xs_max)
    gnmax_2_expectedxs = {
        '1' : ( 1.0, 0.41, 0.43 ),
        '2' : ( 0.88, 0.63, 0.655 ),
        '3' : ( 0.8, 0.765, 0.790 ),
        '6' : ( 0.58, 1.035, 1.055 ),
        '8' : ( 0.52, 1.115, 1.13 ),
        '16' : ( 0.37, 1.258, 1.2628 ),
        'ignore1+2' : ( 1.0, 0.055, 0.07 ),
        'None' : ( 0.37, 1.264, 1.267 ),
    }
    def expected( Gnmax ):
        return gnmax_2_expectedxs[str(Gnmax)]

    failures = []
    for Gnmax in [1,2,3,6,8,16,None,'ignore1+2']:
        mat = matbase.clone()
        if Gnmax == 'ignore1+2':
            mat._unofficial_vdos2sab_ignore( order_low=1,order_high=2 )
        elif Gnmax is not None:
            mat._unofficial_vdos2sab_ignore( order_low=Gnmax+1,
                                             order_high=9999 )
        ref_wl, exp_xsmin, exp_xsmax = expected(Gnmax)
        sc = mat.loadScatter('vdoslux=2;comp=inelas')
        xs = sc.xsect(wl=wls)
        xs_ref = sc.xsect(wl=ref_wl)
        label = ( 'Gn: ignore G1+G2' if Gnmax=='ignore1+2'
                  else ( 'Gn: up to %i'%Gnmax
                         if Gnmax is not None else 'Gn: all') )
        if not ( exp_xsmin <= xs_ref <= exp_xsmax ):
            failures.append(f'FAIL at "{label}": xs({ref_wl}\u00C5)={xs_ref}'
                            f' (required to be in [{exp_xsmin},{exp_xsmax}])')

        if do_plot:
            _ = plt.plot(wls,xs,label=label )
            color = _[0].get_color()
            _ = plt.plot([ref_wl,ref_wl],[exp_xsmin, exp_xsmax],
                         color=color,lw=4,alpha=0.5,label=None)
    if do_plot:
        plt.grid()
        plt.legend()
        plt.xlim(0)
        plt.ylim(0)
        plt.ylabel(f'Inelastic XS (barn/atom) for {ncmatfile}')
        plt.xlabel('Neutron wavelength [\u00C5]')
        plt.show()
    for f in failures:
        print(f)
    if failures:
        raise SystemExit('Failures detected (NB: this test supports --plot)')

if __name__ == '__main__':
    import sys
    main( do_plot = '--plot' in sys.argv[1:] )
