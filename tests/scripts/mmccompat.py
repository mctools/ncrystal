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

# Check the support for the now deprecated MiniMC interface which was in the
# first round of notebooks and presentations. It makes the migration easier if
# we keep them around for a while. The tested and supported functionality is the
# one actually used in our notebooks, so any changes to these notebooks (or
# other peoples code) should migrate to the modern API instead.

#The obsolete functions are available in both _mmc.py and mmc.py:
import NCTestUtils.enable_fpe # noqa F401
import NCTestUtils.reprint_escaped_warnings # noqa F401
import NCrystalDev.mmc as ncmmc
import NCrystalDev._mmc as ncmmc_alt
from NCTestUtils.env import ncsetenv
from NCTestUtils.common import ensure_error
from NCrystalDev.exceptions import NCException

def main(do_plot):

    assert ( ncmmc_alt.runsim_diffraction_pattern
             is ncmmc.runsim_diffraction_pattern )

    mmc_qdp = ncmmc.quick_diffraction_pattern
    assert ncmmc_alt.quick_diffraction_pattern is mmc_qdp

    if not do_plot:
        ncsetenv('FAKEPYPLOT','1')

    #NB: Using "comp=", nthreads=None, nstat=1e3 for reproducible reflogs
    pattern = mmc_qdp('AgBr_sg225_SilverBromide.ncmat;temp=200K;comp=',
                      neutron_energy = '1.8Aa',
                      material_thickness = '1cm',
                      nstat=1e3,
                      nthreads=1)
    pattern.plot_breakdown(rebin_factor=10,logy=True)

    pattern = mmc_qdp('AgBr_sg225_SilverBromide.ncmat;temp=200K;comp=',
                      neutron_energy = '1 eV',
                      material_thickness = '1 cm',
                      nstat=1e3,
                      nthreads=1)
    pattern.plot_breakdown(rebin_factor=50)

    with ensure_error(NCException,
                      '".whatever" is not available on results from the'
                      ' obsolete quick_diffraction_pattern function. Only'
                      ' .plot_breakdown(..) is. Please migrate your code to'
                      ' use the NCrystal.minimc.run(..) function instead'
                      ' (more info at '
                      'https://github.com/mctools/ncrystal/wiki/minimc).'):
        pattern.whatever()

    with ensure_error(NCException,
                      'The runsim_diffraction_pattern(..) function is'
                      ' obsolete. Please migrate your code to use the'
                      ' NCrystal.minimc.run(..) function instead (more info'
                      ' at https://github.com/mctools/ncrystal/wiki/minimc).'):
        ncmmc.runsim_diffraction_pattern('anything','here',it='does not matter')

    #NB: Using "comp=", nthreads=None, nstat=1e3 for reproducible reflogs
    pattern2 = mmc_qdp('void.ncmat',
                      neutron_energy = '1.8Aa',
                      material_thickness = '1cm',
                      nstat=None,
                      nthreads=1)
    nstat = pattern2._real_mmcresult().setup['src']['decoded']['n']
    assert 10000 <= nstat <= 10000000
    print("All ok")

if __name__ == '__main__':
    import sys
    main(do_plot = '--plot' in sys.argv[1:])
