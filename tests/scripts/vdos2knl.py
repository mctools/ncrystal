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
import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev.vdos as nc_vdos
import NCrystalDev.exceptions as nc_exceptions
from NCTestUtils.common import ensure_error

def test( vdos, m, T, do_plot, vdoslux, target_emax = None ):
    if do_plot:
        import NCrystalDev.plot as nc_plot
        nc_plot.plot_vdos(vdos)
    print(f"Calling nc_vdos.extractKnl( m={m:g}, T={T:g}, "
          f"vdoslux={vdoslux}, target_emax={target_emax or 0.0:g} )")
    knl = nc_vdos.extractKnl( vdos,
                              mass_amu = m,
                              temperature = T,
                              target_emax = target_emax,
                              vdoslux = vdoslux,
                              plot = do_plot )
    for k,v in knl.items():
        if hasattr(v,'shape'):
            v = 'NumpyArray( shape=%s )'%(v.shape)
        print(f" Got {k} : {v}")
    print()

def main( do_plot ):
    def t( *a, **kw ):
        test( *a, **kw, do_plot = do_plot )
    vdos = nc_vdos.createVDOSDebye(400.0)
    t( vdos, m = 27.0, T = 0.5, vdoslux = 1 )
    expected_error = ( "VDOS expansion too slow - can not reach E=5000eV"
                       " after 10000 phonon convolutions (likely causes:"
                       " either the target energy value is too high, vdoslux"
                       " too low, the temperature too high, or the VDOS is"
                       " very unusual).")
    t( vdos, m = 27.0, T = 300.0, vdoslux = 1 )
    t( vdos, m = 27.0, T = 300.0, vdoslux = 1, target_emax = 50 )
    t( vdos, m = 27.0, T = 300.0, vdoslux = 4 )
    t( vdos, m = 27.0, T = 300.0, vdoslux = 4, target_emax = 0.5 )
    with ensure_error(nc_exceptions.NCCalcError,expected_error):
        t( vdos, m = 27.0, T = 0.5, vdoslux = 1, target_emax = 5000 )


if __name__ == '__main__':
    import sys
    main( do_plot = '--plot' in sys.argv[1:] )
